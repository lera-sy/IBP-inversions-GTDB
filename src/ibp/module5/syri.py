"""Run minimap2 and SyRI for a single chromosome comparison."""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

from ibp.module5.results import write_standardized

CALL_COLUMNS = [
    "ref_chr", "ref_start", "ref_end",
    "qry_chr", "qry_start", "qry_end",
    "length", "query_span_bp", "syri_id",
    "short_circular_arc_bp",
]


def fasta_length(path):
    with path.open(encoding="ascii") as handle:
        return sum(
            len(line.strip())
            for line in handle
            if not line.startswith(">")
        )


def parse_syri(path, ref_length, circular):
    """Extract parent INV records using native SyRI coordinates."""
    calls = []
    records = 0

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            records += 1

            if fields[9] != "-" or fields[10] != "INV":
                continue

            start, end = int(fields[1]), int(fields[2])
            qstart, qend = int(fields[6]), int(fields[7])
            length = abs(end - start) + 1

            calls.append({
                "ref_chr": fields[0],
                "ref_start": start,
                "ref_end": end,
                "qry_chr": fields[5],
                "qry_start": qstart,
                "qry_end": qend,
                "length": length,
                "query_span_bp": abs(qend - qstart) + 1,
                "syri_id": fields[8],
                "short_circular_arc_bp": (
                    min(length, ref_length - length)
                    if circular else ""
                ),
            })

    return calls, records


def write_calls(path, calls):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=CALL_COLUMNS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(calls)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--query-id", required=True)
    parser.add_argument(
        "--reference-topology",
        choices=("circular", "linear", "unknown"),
        required=True,
    )
    parser.add_argument(
        "--query-topology",
        choices=("circular", "linear", "unknown"),
        required=True,
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--preset", default="asm5")
    parser.add_argument("--bandwidth", default="1k,20k")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument(
        "--syri-env",
        default=None,
        help="Conda environment for SyRI; omit if SyRI is on PATH",
    )

    args = parser.parse_args()

    if args.threads < 1:
        parser.error("--threads must be positive")

    args.outdir.mkdir(parents=True, exist_ok=True)

    calls_path = args.outdir / "calls.tsv"
    standard_path = args.outdir / "inversions.tsv"
    status_path = args.outdir / "status.json"
    log_path = args.outdir / "run.log"
    paf_path = args.outdir / "alignment.paf"
    syri_path = args.outdir / "syri_syri.out"

    for path in (calls_path, status_path, syri_path, standard_path):
        path.unlink(missing_ok=True)

    status = {
        "caller": "syri",
        "reference_assembly": args.reference_id,
        "query_assembly": args.query_id,
        "reference_topology": args.reference_topology,
        "query_topology": args.query_topology,
        "status": "failed",
        "call_count": None,
        "coordinate_system": "SyRI native: 1-based inclusive",
        "parameters": {
            "preset": args.preset,
            "bandwidth": args.bandwidth,
            "threads": args.threads,
        },
    }

    try:
        if args.reference_id == args.query_id:
            raise ValueError("Reference and query must differ")

        if not args.reference.is_file() or not args.query.is_file():
            raise FileNotFoundError("Reference or query FASTA missing")

        ref_length = fasta_length(args.reference)
        status["reference_length_bp"] = ref_length
        status["query_length_bp"] = fasta_length(args.query)

        with log_path.open("w") as log:
            log.write(
                f"Reference: {args.reference}\n"
                f"Query: {args.query}\n"
            )
            log.flush()

            mm_command = [
                args.minimap2,
                "-t", str(args.threads),
                "-x", args.preset,
                "-c", "--eqx",
            ]

            if args.bandwidth:
                mm_command += ["-r", args.bandwidth]

            mm_command += [
                str(args.reference), str(args.query)
            ]

            start = time.perf_counter()

            with paf_path.open("w") as paf:
                subprocess.run(
                    mm_command,
                    stdout=paf,
                    stderr=log,
                    check=True,
                )

            status["minimap_seconds"] = round(
                time.perf_counter() - start, 3
            )

            if paf_path.stat().st_size == 0:
                status.update(
                    status="unresolved",
                    call_count=0,
                    reason="minimap2 produced no alignments",
                )
                write_calls(calls_path, [])
                write_standardized(standard_path, status, [])
                print(f"SyRI: {args.query_id}: unresolved")
                return

            # VSC stores the named SyRI environment outside MAMBA_ROOT_PREFIX.
            # Avoid `micromamba run` per comparison: concurrent launches may
            # serialize on ~/.cache/mamba/proc for minutes without using CPU.
            syri_process_env = None
            if args.syri_env:
                vsc_data = os.environ.get("VSC_DATA")
                direct_syri = (
                    Path(vsc_data) / "conda_envs" / args.syri_env / "bin" / "syri"
                    if vsc_data else None
                )
                if direct_syri is not None and direct_syri.is_file() and os.access(direct_syri, os.X_OK):
                    syri_command = [str(direct_syri)]
                    syri_process_env = dict(os.environ)
                    syri_process_env["PATH"] = (
                        str(direct_syri.parent) + os.pathsep
                        + syri_process_env.get("PATH", "")
                    )
                    log.write(f"SyRI launcher: direct {direct_syri}\n")
                else:
                    syri_command = [
                        "conda" if shutil.which("conda") else "micromamba",
                        "run", "-n", args.syri_env,
                        *(["--no-capture-output"] if shutil.which("conda") else []),
                        "syri",
                    ]
                    log.write("SyRI launcher: conda/micromamba fallback\n")
            else:
                syri_command = ["syri"]
                log.write("SyRI launcher: PATH\n")
            log.flush()

            syri_command += [
                "-c", str(paf_path),
                "-r", str(args.reference),
                "-q", str(args.query),
                "-F", "P",
                "--dir", str(args.outdir),
                "--prefix", "syri_",
            ]

            start = time.perf_counter()

            subprocess.run(
                syri_command,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                env=syri_process_env,
            )

            status["syri_seconds"] = round(
                time.perf_counter() - start, 3
            )

        if not syri_path.is_file():
            raise FileNotFoundError(
                "SyRI did not produce syri_syri.out"
            )

        calls, records = parse_syri(
            syri_path,
            ref_length,
            circular=(
                args.reference_topology == "circular"
                and args.query_topology == "circular"
            ),
        )

        write_calls(calls_path, calls)

        resolved = records > 0

        status.update(
            status="success" if resolved else "unresolved",
            call_count=len(calls),
            syri_records=records,
        )

        write_standardized(standard_path, status, calls)

        print(
            f"SyRI: {args.query_id}: "
            f"{status['status']}, {len(calls)} calls"
        )

    # Worker boundary: record unexpected failures per comparison.
    except Exception as exc:  # noqa: BLE001
        status.update(
            status="failed",
            call_count=None,
            error=f"{type(exc).__name__}: {exc}",
        )
        with log_path.open("a") as log:
            log.write("\nERROR:\n" + traceback.format_exc())

        write_calls(calls_path, [])
        write_standardized(standard_path, status, [])
        print(
            f"SyRI: {args.query_id}: FAILED - {status['error']}",
            file=sys.stderr,
        )

    finally:
        status_path.write_text(
            json.dumps(status, indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
