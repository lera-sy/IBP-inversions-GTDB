"""
Resumable NCBI downloads for large runs (min-20: ~7,800 genomes).

Why: `datasets rehydrate` fails as a whole when NCBI no longer serves a few
files (404 "File unavailable": withdrawn or replaced assemblies), and Snakemake
then deletes the half-downloaded output folder. Here:
  - everything is downloaded into a work folder that is NOT a Snakemake output,
    so progress survives a failed attempt;
  - each round retries only the files still missing;
  - files NCBI cannot serve are listed in a TSV and skipped; genome QC then
    marks those assemblies 'missing_or_multiple_fasta_files' and drops them;
  - the step fails only if more than --max-missing of the accessions are missing;
  - on success the work folder is moved into place as the rule output.

  genomes  dehydrated package -> rehydrate in rounds
  reports  sequence reports in batches (per accession for a batch that fails)
"""

import argparse
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path


def nonempty(path):
    return path.is_file() and path.stat().st_size > 0


def read_accessions(path):
    return [a.strip() for a in Path(path).read_text().split() if a.strip()]


def finish(work, out, accessions, have, missing_path, max_missing, what):
    missing = [a for a in accessions if not have(a)]
    with open(missing_path, "w") as fh:
        fh.write("accession\treason\n")
        for a in missing:
            fh.write(f"{a}\t{what}_unavailable_from_ncbi\n")
    frac = len(missing) / max(len(accessions), 1)
    print(f"{what}: {len(accessions) - len(missing)}/{len(accessions)} present, {len(missing)} missing "
          f"({frac:.1%}); list in {missing_path}", flush=True)
    if frac > max_missing:
        print(f"More than {max_missing:.0%} missing: stopping. Progress is kept in {work}; "
              "resubmit to continue.", file=sys.stderr)
        sys.exit(1)
    if out.exists():
        shutil.rmtree(out)
    work.rename(out)


# ---------------------------------------------------------------- genomes

def genomes(a):
    work, out = Path(a.work), Path(a.out)
    work.mkdir(parents=True, exist_ok=True)
    nd = work / "ncbi_dataset"
    all_fetch = nd / "fetch_all.txt"
    if not all_fetch.is_file():
        with zipfile.ZipFile(a.package) as z:
            z.extractall(work)
        shutil.copy(nd / "fetch.txt", all_fetch)
    lines = [ln for ln in all_fetch.read_text().splitlines() if ln.strip()]

    def target(line):
        rel = line.split("\t")[-1].strip()
        return work / rel if rel.startswith("ncbi_dataset/") else nd / rel

    for rnd in range(1, a.rounds + 1):
        todo = [ln for ln in lines if not nonempty(target(ln))]
        print(f"round {rnd}: {len(lines) - len(todo)}/{len(lines)} files present, {len(todo)} to fetch", flush=True)
        if not todo:
            break
        (nd / "fetch.txt").write_text("\n".join(todo) + "\n")
        r = subprocess.run(["datasets", "rehydrate", "--directory", str(work), "--max-workers", str(a.workers)])
        print(f"round {rnd}: rehydrate exit code {r.returncode}", flush=True)
        if rnd < a.rounds:
            time.sleep(30 * rnd)

    accessions = read_accessions(a.accessions)
    have = lambda acc: any(nonempty(p) for p in (nd / "data" / acc).glob("*_genomic.fna"))  # noqa: E731
    finish(work, out, accessions, have, a.missing, a.max_missing, "genome")


# ---------------------------------------------------------------- sequence reports

def fetch_reports(accs, work, tag, tries=3):
    lst = work / f"{tag}.txt"
    lst.write_text("\n".join(accs) + "\n")
    zp = work / f"{tag}.zip"
    for t in range(tries):
        r = subprocess.run(["datasets", "download", "genome", "accession", "--inputfile", str(lst),
                            "--include", "seq-report", "--filename", str(zp)])
        if r.returncode == 0 and zp.is_file():
            try:
                with zipfile.ZipFile(zp) as z:
                    z.extractall(work)
                zp.unlink()
                return True
            except zipfile.BadZipFile:
                pass
        time.sleep(20 * (t + 1))
    return False


def reports(a):
    work, out = Path(a.work), Path(a.out)
    work.mkdir(parents=True, exist_ok=True)
    accessions = read_accessions(a.accessions)
    have = lambda acc: nonempty(work / "ncbi_dataset/data" / acc / "sequence_report.jsonl")  # noqa: E731
    batches = [accessions[i:i + a.batch] for i in range(0, len(accessions), a.batch)]
    for i, batch in enumerate(batches, 1):
        todo = [acc for acc in batch if not have(acc)]
        if not todo:
            continue
        ok = fetch_reports(todo, work, f"batch_{i:03d}")
        print(f"batch {i}/{len(batches)}: {'ok' if ok else 'failed, trying one by one'}", flush=True)
        if not ok:
            for acc in [x for x in todo if not have(x)]:
                fetch_reports([acc], work, "single", tries=2)
    finish(work, out, accessions, have, a.missing, a.max_missing, "report")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="what", required=True)
    g = sub.add_parser("genomes")
    g.add_argument("--package", required=True)
    g.add_argument("--rounds", type=int, default=4)
    g.add_argument("--workers", type=int, default=2)
    r = sub.add_parser("reports")
    r.add_argument("--batch", type=int, default=500)
    for p in (g, r):
        p.add_argument("--accessions", required=True)
        p.add_argument("--work", required=True)
        p.add_argument("--out", required=True)
        p.add_argument("--missing", required=True)
        p.add_argument("--max-missing", type=float, default=0.05)
    a = ap.parse_args()
    genomes(a) if a.what == "genomes" else reports(a)


if __name__ == "__main__":
    main()
