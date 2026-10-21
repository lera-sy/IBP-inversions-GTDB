"""Run the existing Leen algorithm on one chromosome pair."""

import argparse
import csv
import json
import sys
import traceback
from pathlib import Path

from ibp.module5.leen_core import (
    call_inversions,
    map_windows,
    read_fasta,
    read_hits,
    write_windows,
)
from ibp.module5.results import write_standardized

CALL_COLUMNS = [
    "ref_start",
    "ref_end",
    "length",
    "fraction_of_chromosome",
    "ref_intervals",
    "n_intervals",
    "spans_reference_origin",
    "call_type",
]


def write_calls(path, calls):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=CALL_COLUMNS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(calls)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--query-id", required=True)
    parser.add_argument("--reference-topology", required=True)
    parser.add_argument("--query-topology", required=True)
    parser.add_argument("--outdir", type=Path, required=True)

    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--preset", default="asm20")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--window", type=int, default=2000)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-size", type=int, default=6000)
    parser.add_argument("--max-jump", type=int, default=20000)

    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    calls_path = args.outdir / "calls.tsv"
    standard_path = args.outdir / "inversions.tsv"
    status_path = args.outdir / "status.json"
    log_path = args.outdir / "run.log"

    # Avoid retaining stale successful outputs after a failed rerun.
    calls_path.unlink(missing_ok=True)
    standard_path.unlink(missing_ok=True)
    status_path.unlink(missing_ok=True)

    status = {
        "caller": "leen",
        "reference_assembly": args.reference_id,
        "query_assembly": args.query_id,
        "status": "failed",
        "call_count": None,
        "parameters": {
            "preset": args.preset,
            "window": args.window,
            "min_mapq": args.min_mapq,
            "min_size": args.min_size,
            "max_jump": args.max_jump,
            "threads": args.threads,
        },
    }

    try:
        if args.reference_id == args.query_id:
            raise ValueError("Reference and query must differ")

        if not args.reference.is_file() or not args.query.is_file():
            raise FileNotFoundError("Reference or query FASTA missing")

        if (
            args.reference_topology != "circular"
            or args.query_topology != "circular"
        ):
            status.update(
                status="unsupported",
                call_count=0,
                reason="Current Leen core requires circular chromosomes",
            )
            write_calls(calls_path, [])
            write_standardized(standard_path, status, [])
            log_path.write_text(status["reason"] + "\n")
            print("UNSUPPORTED: chromosome topology")
            return

        work = args.outdir / "work"
        work.mkdir(exist_ok=True)

        windows = work / "windows.fna"
        paf = work / "alignment.paf"

        # The inherited mapper skips an existing nonempty PAF.
        # Delete it so this run always uses current inputs.
        paf.unlink(missing_ok=True)

        query_length = write_windows(
            str(args.query), str(windows), args.window
        )

        map_windows(
            str(args.reference),
            str(windows),
            str(paf),
            minimap2=args.minimap2,
            preset=args.preset,
            threads=args.threads,
        )

        hits, reference_length = read_hits(
            str(paf), args.window, args.min_mapq
        )

        if not reference_length:
            reference_length = len(
                next(read_fasta(str(args.reference)))[1]
            )

        calls, info = call_inversions(
            hits=hits,
            q_len=query_length,
            r_len=reference_length,
            window=args.window,
            min_size=args.min_size,
            max_jump=args.max_jump,
        )

        status.update(
            status="success" if info["resolved"] else "unresolved",
            call_count=len(calls),
            reference_length_bp=reference_length,
            query_length_bp=query_length,
            algorithm_info=info,
        )

        write_calls(calls_path, calls)
        write_standardized(standard_path, status, calls)

        log_path.write_text(
            f"Reference: {args.reference_id}\n"
            f"Query: {args.query_id}\n"
            f"Accepted windows: {len(hits)}\n"
            f"Status: {status['status']}\n"
            f"Inversion calls: {len(calls)}\n"
        )

        print(
            f"Leen: {args.query_id}: "
            f"{status['status']}, {len(calls)} calls"
        )

    # Worker boundary: record unexpected failures per comparison.
    except Exception as exc:  # noqa: BLE001
        status.update(
            status="failed",
            call_count=None,
            error=f"{type(exc).__name__}: {exc}",
        )
        log_path.write_text(traceback.format_exc())
        write_calls(calls_path, [])
        write_standardized(standard_path, status, [])
        print(
            f"Leen: {args.query_id}: FAILED - {status['error']}",
            file=sys.stderr,
        )

    finally:
        status_path.write_text(
            json.dumps(status, indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
