#!/usr/bin/env python3
"""Do inversion breakpoints sit between inverted repeats? Evidence independent of the callers.

Most bacterial inversions arise by recombination between two copies of a repeat
(IS element, rRNA operon) lying in opposite orientation. A real inversion should
then show the same sequence at both breakpoints, reverse-complemented.

For every inversion family in events.tsv (output of final_inversion_model.py):
  1. take its two breakpoints on the species' reference chromosome;
  2. cut a window around each breakpoint (default +-10 kb, never past halfway
     to the other breakpoint, so the two windows never overlap);
  3. align window 1 against window 2 with minimap2 (asm20); a match counts only
     when both copies lie at their breakpoint (--max-offset, default 3 kb), and:
       inverted  - a match >= --min-repeat bp, >= --min-identity, opposite strand
       direct    - only a same-strand match (points to deletion/duplication)
       none      - no match;
  4. optionally label the repeat as an rRNA operon (barrnap) or IS element
     (ISEScan GFF), when those annotations are available.
A background control repeats step 2-3 for random breakpoint pairs of the same
length on the same chromosome, so the rate expected by chance in a repeat-rich
genome is known.

Usage (from the repository root, ibp-pipeline environment):
  python scripts/steps/analysis/breakpoint_repeats.py \
      --events results/final_model_all_species/events.tsv \
      --pilot-root results/pilot \
      --out results/breakpoint_repeats --threads 4
Optional: --barrnap (runs barrnap if installed), --isescan-dir DIR (reads
<reference>.fna.gff from an ISEScan run).
"""

import argparse
import gzip
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        name, parts = None, []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    break  # first record = the selected chromosome
                name = line[1:].split()[0]
            else:
                parts.append(line.strip())
    return name, "".join(parts).upper()


def reference_for(species, pilot_root, repo_root):
    """Path of the Mash-selected reference FASTA of one species."""
    txt = pilot_root / species / "reference_selection" / "selected_reference.txt"
    if not txt.is_file():
        raise FileNotFoundError(f"No reference selection for {species}: {txt}")
    info = dict(line.split("\t", 1) for line in txt.read_text().splitlines() if "\t" in line)
    path = Path(info["path"])
    candidates = [path] if path.is_absolute() else [repo_root / path, pilot_root.parent.parent / path]
    for candidate in candidates:
        if candidate.is_file():
            return info["genome"], candidate
    # fall back: search the species folder for the genome file
    hits = sorted((pilot_root / species).rglob(f"{info['genome']}.fna"))
    if hits:
        return info["genome"], hits[0]
    raise FileNotFoundError(f"Reference FASTA {info['path']} not found for {species}")


def circular_slice(seq, start, end):
    """seq[start:end] on a circle; start/end may be negative or past the end."""
    n = len(seq)
    start, end = start % n, end % n
    return seq[start:end] if start < end else seq[start:] + seq[:end]


def windows(p1, p2, n, width):
    """Windows around two breakpoints that never overlap.

    Arc A runs p1 -> p2 forward, arc B p2 -> p1 forward. Each breakpoint's window
    reaches into each arc by at most half that arc's length.
    """
    a = (p2 - p1) % n
    b = n - a
    in_a, in_b = min(width, a // 2), min(width, b // 2)
    w1 = (p1 - in_b, p1 + in_a)       # left of p1 lies in arc B, right in arc A
    w2 = (p2 - in_a, p2 + in_b)
    return w1, w2, (in_b, in_a)        # breakpoint position inside window 1 and window 2


def run_minimap2(minimap2, threads, flanks1, flanks2, workdir):
    q, t, paf = workdir / "q.fa", workdir / "t.fa", workdir / "hits.paf"
    q.write_text("".join(f">{k}\n{s}\n" for k, s in flanks1.items() if s))
    t.write_text("".join(f">{k}\n{s}\n" for k, s in flanks2.items() if s))
    with paf.open("w") as out:
        subprocess.run([minimap2, "-x", "asm20", "-c", "-P", "-t", str(threads), str(t), str(q)],
                       stdout=out, stderr=subprocess.DEVNULL, check=True)
    hits = []
    with paf.open() as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[0] != f[5]:
                continue  # only window 1 vs window 2 of the same breakpoint pair
            hits.append({"key": f[0], "strand": f[4], "q_start": int(f[2]), "q_end": int(f[3]),
                         "t_start": int(f[7]), "t_end": int(f[8]),
                         "matches": int(f[9]), "aln_len": int(f[10])})
    return pd.DataFrame(hits, columns=["key", "strand", "q_start", "q_end", "t_start", "t_end",
                                       "matches", "aln_len"])


def distance_to(pos, start, end):
    """0 if pos lies inside [start, end), else the distance to the nearest edge."""
    return 0 if start <= pos < end else min(abs(pos - start), abs(pos - end))


def classify(hits, min_repeat, min_identity, bp_offsets, max_offset):
    """Per key: inverted / direct / none, with the longest qualifying match.

    A match only counts when BOTH copies lie at their breakpoint (within
    max_offset bp): in repeat-rich genomes a window almost always holds some
    repeat copy, so a copy elsewhere in the window is no evidence.
    bp_offsets[key] = position of the breakpoint inside window 1 and window 2.
    """
    out = {}
    if hits.empty:
        return out
    hits = hits.assign(identity=hits["matches"] / hits["aln_len"])
    near = [distance_to(bp_offsets[k][0], qs, qe) <= max_offset and distance_to(bp_offsets[k][1], ts, te) <= max_offset
            for k, qs, qe, ts, te in zip(hits["key"], hits["q_start"], hits["q_end"], hits["t_start"], hits["t_end"])]
    hits = hits[near]
    good = hits[(hits["aln_len"] >= min_repeat) & (hits["identity"] >= min_identity)]
    for key, h in good.groupby("key"):
        inv = h[h["strand"] == "-"]
        best = (inv if not inv.empty else h).sort_values("aln_len").iloc[-1]
        out[key] = {"repeat_class": "inverted" if not inv.empty else "direct",
                    "repeat_length_bp": int(best["aln_len"]),
                    "repeat_identity": round(float(best["identity"]), 4),
                    "_q": (int(best["q_start"]), int(best["q_end"])),
                    "_t": (int(best["t_start"]), int(best["t_end"]))}
    return out


def read_gff_intervals(path, feature_filter=None):
    intervals = []
    if path is None or not Path(path).is_file():
        return intervals
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 9 or (feature_filter and not feature_filter(f)):
                continue
            intervals.append((int(f[3]) - 1, int(f[4])))
    return intervals


def run_barrnap(barrnap, fasta, threads, workdir):
    gff = workdir / "rrna.gff"
    with gff.open("w") as out:
        result = subprocess.run([barrnap, "--threads", str(threads), "--quiet", str(fasta)],
                                stdout=out, stderr=subprocess.DEVNULL, check=False)
    if result.returncode != 0:
        return []
    return read_gff_intervals(gff, lambda f: "rRNA" in f[2])


def overlaps(interval_on_ref, features, n):
    """Does a reference interval (start, end; circular) overlap any feature?"""
    s, e = interval_on_ref
    s, e = s % n, e % n
    parts = [(s, e)] if s < e else [(s, n), (0, e)]
    return any(a < fe and fs < b for a, b in parts for fs, fe in features)


def label_repeat(rec, w1, w2, n, rrna, isel):
    """rRNA / IS / unknown, from where the repeat copy in window 1 lies."""
    if "_q" not in rec:
        return ""
    copy1 = (w1[0] + rec["_q"][0], w1[0] + rec["_q"][1])
    copy2 = (w2[0] + rec["_t"][0], w2[0] + rec["_t"][1])
    if rrna and (overlaps(copy1, rrna, n) or overlaps(copy2, rrna, n)):
        return "rRNA"
    if isel and (overlaps(copy1, isel, n) or overlaps(copy2, isel, n)):
        return "IS"
    return "unannotated" if (rrna or isel) else ""


def evidence_group(row):
    if row["tier"] == "medium":
        return "medium_syri_only" if pd.notna(row.get("syri_family")) else "medium_leen_only"
    return row["tier"]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--events", type=Path, required=True)
    p.add_argument("--pilot-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--repo-root", type=Path, default=Path.cwd())
    p.add_argument("--window", type=int, default=10000, help="bp on each side of a breakpoint")
    p.add_argument("--min-repeat", type=int, default=500)
    p.add_argument("--max-offset", type=int, default=3000,
                   help="max distance (bp) between a repeat copy and its breakpoint")
    p.add_argument("--min-identity", type=float, default=0.90)
    p.add_argument("--controls", type=int, default=2, help="random breakpoint pairs per event")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--barrnap", action="store_true", help="annotate rRNA operons with barrnap")
    p.add_argument("--isescan-dir", type=Path, help="folder with ISEScan output (<genome>.fna.gff)")
    p.add_argument("--only", action="append", default=[], help="species to run (repeatable)")
    args = p.parse_args(argv)

    if shutil.which(args.minimap2) is None:
        sys.exit(f"minimap2 not found ({args.minimap2}); activate ibp-pipeline")
    barrnap = shutil.which("barrnap") if args.barrnap else None
    if args.barrnap and barrnap is None:
        print("barrnap not found: rRNA labels skipped", flush=True)

    events = pd.read_csv(args.events, sep="\t")
    if args.only:
        events = events[events["species"].isin(args.only)]
    events = events.reset_index(drop=True)
    events["evidence_group"] = events.apply(evidence_group, axis=1)
    rng = random.Random(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)

    rows, controls = [], []
    for n_sp, (species, ev) in enumerate(events.groupby("species"), 1):
        try:
            genome, fasta = reference_for(species, args.pilot_root, args.repo_root)
        except FileNotFoundError as exc:
            print(f"[{n_sp}] {species}: SKIPPED ({exc})", flush=True)
            rows += [{**r._asdict(), "repeat_class": "no_reference"} for r in ev.itertuples(index=False)]
            continue
        _, seq = read_fasta(fasta)
        n = len(seq)
        if abs(n - int(ev["reference_length_bp"].iloc[0])) > 1000:
            print(f"[{n_sp}] {species}: WARNING reference length {n} != {ev['reference_length_bp'].iloc[0]}",
                  flush=True)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            rrna = run_barrnap(barrnap, fasta, args.threads, tmp) if barrnap else []
            isel = []
            if args.isescan_dir:  # ISEScan writes <genome file name>.gff somewhere below its output folder
                gffs = sorted(args.isescan_dir.rglob(f"{fasta.name}.gff"))
                isel = read_gff_intervals(gffs[0]) if gffs else []

            flanks1, flanks2, win, offsets = {}, {}, {}, {}
            for i, r in ev.iterrows():
                w1, w2, off = windows(int(r["ref_start"]), int(r["ref_end"]), n, args.window)
                key = f"e{i}"
                flanks1[key], flanks2[key], win[key] = circular_slice(seq, *w1), circular_slice(seq, *w2), (w1, w2)
                offsets[key] = off
                for c in range(args.controls):  # same breakpoint distance, random position
                    shift = rng.randrange(n)
                    p1, p2 = (int(r["ref_start"]) + shift) % n, (int(r["ref_end"]) + shift) % n
                    c1, c2, coff = windows(p1, p2, n, args.window)
                    ckey = f"c{i}_{c}"
                    flanks1[ckey], flanks2[ckey], win[ckey] = circular_slice(seq, *c1), circular_slice(seq, *c2), (c1, c2)
                    offsets[ckey] = coff
            found = classify(run_minimap2(args.minimap2, args.threads, flanks1, flanks2, tmp),
                             args.min_repeat, args.min_identity, offsets, args.max_offset)

        for i, r in ev.iterrows():
            rec = found.get(f"e{i}", {"repeat_class": "none"})
            w1, w2 = win[f"e{i}"]
            rows.append({**r.to_dict(), "reference": genome,
                         "repeat_class": rec["repeat_class"],
                         "repeat_length_bp": rec.get("repeat_length_bp"),
                         "repeat_identity": rec.get("repeat_identity"),
                         "repeat_type": label_repeat(rec, w1, w2, n, rrna, isel),
                         "rrna_annotated": bool(rrna), "is_annotated": bool(isel)})
            for c in range(args.controls):
                controls.append({"species": species, "evidence_group": r["evidence_group"],
                                 "repeat_class": found.get(f"c{i}_{c}", {"repeat_class": "none"})["repeat_class"]})
        got = sum(found.get(f"e{i}", {}).get("repeat_class") == "inverted" for i in ev.index)
        print(f"[{n_sp}] {species}: {len(ev)} families, {got} with inverted repeats", flush=True)

    table = pd.DataFrame(rows)
    table.to_csv(args.out / "breakpoint_repeats.tsv", sep="\t", index=False)
    ctrl = pd.DataFrame(controls)

    def shares(df, by):
        tab = pd.crosstab(df[by], df["repeat_class"], normalize="index").mul(100).round(1)
        tab.insert(0, "n", df.groupby(by).size())
        return tab

    order = ["high", "medium_syri_only", "medium_leen_only", "low"]
    real = shares(table[table["repeat_class"] != "no_reference"], "evidence_group").reindex(order).dropna(how="all")
    real.to_csv(args.out / "summary_by_group.tsv", sep="\t")
    print("\nEvents, % per class (inverted repeat = support for a real inversion):")
    print(real.to_string())
    if not ctrl.empty:
        background = shares(ctrl, "evidence_group").reindex(order).dropna(how="all")
        background.to_csv(args.out / "background_by_group.tsv", sep="\t")
        print("\nBackground: random breakpoint pairs of the same lengths, same genomes:")
        print(background.to_string())
    inverted = table[(table["repeat_class"] == "inverted") & table["repeat_type"].fillna("").astype(bool)]
    if not inverted.empty:
        types = pd.crosstab(inverted["evidence_group"], inverted["repeat_type"]).reindex(order).dropna(how="all")
        types.to_csv(args.out / "repeat_types_by_group.tsv", sep="\t")
        print("\nRepeat type at inverted-repeat breakpoints:")
        print(types.to_string())
    print(f"\nWritten to {args.out}")


if __name__ == "__main__":
    main()
