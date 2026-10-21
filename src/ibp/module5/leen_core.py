"""Existing Leen inversion algorithm.

Extracted from feature/genome-parallelism at c55a96a24bba727f1f07af0b12c9239b7b40a706.
Currently implements circular-chromosome logic.
"""

import gzip
import math
import os
import subprocess
import time


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_windows(chromosome, out_path, window):
    """Cut the chromosome into consecutive windows named by their index."""
    _, sequence = next(read_fasta(chromosome))
    with open(out_path, "w") as out:
        out.writelines(f">{i // window}\n{sequence[i:i + window]}\n"
                       for i in range(0, len(sequence), window) if len(sequence) - i >= window // 2)
    return len(sequence)


def map_windows(reference, windows, paf, minimap2="minimap2", preset="asm20", threads=4):
    if os.path.exists(paf) and os.path.getsize(paf) > 0:
        return
    with open(paf, "w") as out:
        subprocess.run([minimap2, "-x", preset, "-t", str(threads), reference, windows],
                       stdout=out, stderr=subprocess.DEVNULL, check=True)


def read_hits(paf, window, min_mapq):
    """Best unique placement of each window: list of dicts sorted along the genome."""
    best, r_len = {}, 0
    with open(paf) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if "tp:A:S" in f[12:] or int(f[11]) < min_mapq:
                continue  # repeat copies and ambiguous placements
            index, span = int(f[0]), int(f[3]) - int(f[2])
            if span < int(f[1]) / 2:
                continue  # less than half of the window placed
            r_len = int(f[6])
            if index not in best or span > best[index]["span"]:
                best[index] = {"index": index, "q_mid": index * window + (int(f[2]) + int(f[3])) / 2,
                               "q_start": index * window + int(f[2]), "q_end": index * window + int(f[3]),
                               "r_start": int(f[7]), "r_end": int(f[8]), "strand": f[4], "span": span}
    return [best[i] for i in sorted(best)], r_len


def colinear(a, b, q_len, r_len, max_jump):
    """True if window b continues window a: same strand, same distance in both genomes."""
    if a["strand"] != b["strand"]:
        return False
    sign = 1 if a["strand"] == "+" else -1
    q_step = (b["q_mid"] - a["q_mid"]) % q_len
    r_step = (((b["r_start"] + b["r_end"]) - (a["r_start"] + a["r_end"])) / 2 * sign) % r_len
    difference = (r_step - q_step + r_len / 2) % r_len - r_len / 2
    return abs(difference) <= max_jump


def build_blocks(hits, q_len, r_len, max_jump, min_size, window):
    """Synteny blocks in genome order; blocks smaller than min_size are removed."""
    while True:
        groups = []
        for hit in hits:
            if groups and colinear(groups[-1][-1], hit, q_len, r_len, max_jump):
                groups[-1].append(hit)
            else:
                groups.append([hit])
        # the chromosome is circular: the last block may continue into the first
        if len(groups) > 1 and colinear(groups[-1][-1], groups[0][0], q_len, r_len, max_jump):
            groups[0] = groups.pop() + groups[0]
        kept = [g for g in groups if len(g) * window >= min_size]
        if len(kept) == len(groups):
            break
        hits = sorted((h for g in kept for h in g), key=lambda h: h["index"])
    blocks = []
    for group in groups:
        first, last = group[0], group[-1]
        forward = first["strand"] == "+"
        blocks.append({"strand": first["strand"], "q_start": first["q_start"], "q_end": last["q_end"],
                       "ref_start": first["r_start"] if forward else last["r_start"],
                       "ref_end": last["r_end"] if forward else first["r_end"],
                       "length": len(group) * window})
    return blocks


def signed_order(blocks, r_len):
    """
    Number the blocks 1..n in reference order (block 1 = the one holding the reference
    origin, or the first after it) and return their signed order around the genome,
    starting at block 1 on the + strand (the whole order is flipped if needed).
    """
    wrapping = [i for i, b in enumerate(blocks) if b["ref_start"] > b["ref_end"]]
    anchor = wrapping[0] if wrapping else min(range(len(blocks)), key=lambda i: blocks[i]["ref_start"])
    offset = blocks[anchor]["ref_start"]
    by_reference = sorted(range(len(blocks)), key=lambda i: (blocks[i]["ref_start"] - offset) % r_len)
    for rank, i in enumerate(by_reference, 1):
        blocks[i]["rank"] = rank
    order = [b["rank"] if b["strand"] == "+" else -b["rank"] for b in blocks]
    if -1 in order:
        order = [-x for x in reversed(order)]
    start = order.index(1)
    return tuple(order[start:] + order[:start])


def breakpoints(remainder, n):
    framed = (1,) + remainder + (n + 1,)
    return [k for k in range(len(framed) - 1) if framed[k + 1] - framed[k] != 1]


def sort_by_inversions(order, time_limit=15, max_solutions=50):
    """
    Smallest sets of inversions that sort `order` into 1, 2, ..., n.
    Block 1 is held fixed (on a circle, inverting a segment containing it is the same
    as inverting the rest). Returns (solutions, complete): each solution is a tuple of
    inversions, each inversion a tuple of the block numbers it contains. `complete` is
    False if the search stopped early (time or number of solutions), in which case
    there may be more solutions of the same size. Returns (None, False) if no solution
    was found within the time limit.
    """
    n, start = len(order), order[1:]
    deadline = time.monotonic() + time_limit
    solutions, failed, stopped = set(), {}, False

    def search(state, path, depth_left):
        """True if at least one solution lies below this state."""
        nonlocal stopped
        cuts = breakpoints(state, n)
        if not cuts:
            solutions.add(tuple(sorted(path)))
            return True
        if depth_left < math.ceil(len(cuts) / 2) or failed.get(state, -1) >= depth_left:
            return False
        if stopped or len(solutions) >= max_solutions or time.monotonic() > deadline:
            stopped = True
            return True  # unknown: do not record this state as a dead end
        found = False
        for a, i in enumerate(cuts):
            for j in cuts[a + 1:]:  # invert state[i:j], cutting only where order is already broken
                segment = state[i:j]
                inverted = state[:i] + tuple(-x for x in reversed(segment)) + state[j:]
                if search(inverted, path + (tuple(sorted(abs(x) for x in segment)),), depth_left - 1):
                    found = True
        if not found:
            failed[state] = depth_left
        return found

    lower = math.ceil(len(breakpoints(start, n)) / 2)
    for depth in range(lower, lower + 4):
        search(start, (), depth)
        if solutions:
            return sorted(solutions), not stopped
        if stopped:
            break
    return None, False


def clean_reversals(order):
    """
    Peel off inversions that restore order on BOTH sides when flipped back.

    `order` is the signed circular order starting with block 1 (see signed_order).
    A segment is a clean inversion if, after reversing it, its left and right
    neighbours both become consecutive (the two breakpoints disappear). This is the
    signature of a simple inversion; nested inversions are found one level at a
    time (the smallest clean segment first). Moved blocks (translocations) never
    satisfy it. Returns (inversions, remaining order); each inversion is a tuple of
    the block numbers it contains.
    """
    n = len(order)
    framed = [1] + list(order[1:]) + [n + 1]
    inversions = []
    while True:
        options = []
        for i in range(1, len(framed) - 1):
            for j in range(i, len(framed) - 1):
                left, right = framed[i - 1], framed[j + 1]
                if -framed[j] - left == 1 and right + framed[i] == 1:
                    options.append((j - i, i, j))
        if not options:
            return inversions, tuple(framed[:-1])
        _, i, j = min(options)
        segment = framed[i:j + 1]
        inversions.append(tuple(sorted(abs(x) for x in segment)))
        framed[i:j + 1] = [-x for x in reversed(segment)]


def _strips(order):
    """Maximal runs of consecutive signed blocks (k, k+1, ... or -k, -(k-1), ...)."""
    strips = [[order[0]]]
    for x in order[1:]:
        if x - strips[-1][-1] == 1:
            strips[-1].append(x)
        else:
            strips.append([x])
    return strips


def place_remaining(order, block_length):
    """
    Split what is left after clean_reversals into blocks that are in place and
    blocks that moved, like a synteny caller: the largest (by length) set of runs
    that appear in reference order stays in place; flipped runs within that set are
    inversions in place, every other run moved (translocation, possibly flipped).
    All circular relabellings are tried, so block 1 itself may be the one that moved.
    Returns (inversions, moved_blocks, ambiguous).
    """
    n = len(order)
    best = []
    for shift in range(n):
        relabel = {r: (r - 1 - shift) % n + 1 for r in range(1, n + 1)}
        seq = [relabel[abs(x)] * (1 if x > 0 else -1) for x in order]
        if 1 not in seq:
            seq = [-x for x in reversed(seq)]
        start = seq.index(1)
        seq = seq[start:] + seq[:start]
        back = {v: k for k, v in relabel.items()}
        strips = _strips(seq)
        info = []
        for s in strips:
            w = sum(block_length[back[abs(x)]] for x in s)
            # score: first the collinear (unflipped) length kept in place, then the total,
            # so a flipped reading of most of the circle never wins over the backbone
            info.append((min(abs(x) for x in s), max(abs(x) for x in s),
                         (w if s[0] > 0 else 0, w), s))
        # weighted longest chain of runs with increasing, non-overlapping ranks
        score, prev = [], []
        for k, (lo, _hi, w, _s) in enumerate(info):
            options = [((score[m][0] + w[0], score[m][1] + w[1]), m)
                       for m in range(k) if info[m][1] < lo]
            top = max(options, default=(w, -1))
            score.append(top[0])
            prev.append(top[1])
        end = max(range(len(info)), key=lambda k: score[k])
        chain, k = set(), end
        while k != -1:
            chain.add(k)
            k = prev[k]
        inversions = tuple(sorted(tuple(sorted(back[abs(x)] for x in info[k][3]))
                                  for k in chain if info[k][3][0] < 0))
        moved = tuple(sorted(back[abs(x)] for k in range(len(info)) if k not in chain
                             for x in info[k][3]))
        best.append((score[end], inversions, moved))
    top = max(b[0] for b in best)
    winners = {(inv, mv) for sc, inv, mv in best if sc == top}
    inversions, moved = min(winners)
    return list(inversions), list(moved), len(winners) > 1


def describe_inversion(content, blocks, r_len):
    """Length and reference position of one inversion, given the blocks it contains."""
    by_rank = {b["rank"]: b for b in blocks}
    total = sum(b["length"] for b in blocks)
    length = sum(by_rank[r]["length"] for r in content)
    ranks = set(content)
    if length > total / 2:  # report the shorter side of the circle
        ranks, length = set(by_rank) - ranks, total - length
    n = len(by_rank)
    # group consecutive blocks (n and 1 are neighbours on the circle) into reference intervals
    starts = sorted(r for r in ranks if (r - 2) % n + 1 not in ranks) or [min(ranks)]
    intervals = []
    for first in starts:
        last = first
        while last % n + 1 in ranks and last % n + 1 != first:
            last = last % n + 1
        intervals.append((by_rank[first]["ref_start"], by_rank[last]["ref_end"]))
    return {"ref_start": intervals[0][0], "ref_end": intervals[-1][1], "length": length,
            "fraction_of_chromosome": round(length / r_len, 4),
            "ref_intervals": ";".join(f"{a}-{b}" for a, b in intervals),
            "n_intervals": len(intervals),
            "spans_reference_origin": any(a > b for a, b in intervals)}


def call_inversions(hits, q_len, r_len, window, min_size, max_jump, max_blocks=20, time_limit=15):
    """
    Inversions of one genome relative to the reference. Returns (calls, info).

    Only rearrangements that look like inversions are called: clean inversions
    (flipping them back restores order on both sides, nested ones included) and
    flipped runs that are still in their reference position. Blocks that moved are
    counted as translocations in `info` and are NOT reported as inversions.
    max_blocks and time_limit are kept for compatibility; this classification is
    polynomial, so no genome is left unresolved for being too rearranged.
    """
    n_windows = math.ceil(q_len / window)
    info = {"mapped_fraction": round(len(hits) / n_windows, 3), "n_blocks": 0, "resolved": True,
            "ambiguous": False, "n_clean_inversions": 0, "n_in_place_inversions": 0,
            "n_moved_blocks": 0, "moved_bp": 0, "rearrangement_class": "collinear"}
    if not hits:
        info["resolved"] = False
        info["rearrangement_class"] = "unmapped"
        return [], info
    blocks = build_blocks(hits, q_len, r_len, max_jump, min_size, window)
    info["n_blocks"] = len(blocks)
    if len(blocks) < 2:
        return [], info
    order = signed_order(blocks, r_len)
    length = {b["rank"]: b["length"] for b in blocks}
    clean, remaining = clean_reversals(order)
    in_place, moved, ambiguous = place_remaining(remaining, length)
    shared = []
    if moved and all(x < 0 for x in remaining if abs(x) in moved):
        # Every leftover moved segment is flipped. Two inversions sharing a breakpoint
        # leave exactly this pattern (indistinguishable from an inverted translocation),
        # so explain the leftovers with the fewest inversions and label them.
        solutions, complete = sort_by_inversions(remaining, time_limit)
        if solutions:
            options = [[describe_inversion(c, blocks, r_len) for c in sol] for sol in solutions]
            shared = min(options, key=lambda sol: (sum(e["length"] for e in sol),
                                                   [e["ref_intervals"] for e in sol]))
            ambiguous = ambiguous or len(solutions) > 1 or not complete
            in_place, moved = [], []
    calls = ([{**describe_inversion(c, blocks, r_len), "call_type": "clean"} for c in clean]
             + [{**describe_inversion(c, blocks, r_len), "call_type": "in_place"} for c in in_place]
             + [{**e, "call_type": "shared_breakpoint"} for e in shared])
    info.update({
        "ambiguous": ambiguous,
        "n_clean_inversions": len(clean),
        "n_in_place_inversions": len(in_place),
        "n_moved_blocks": len(moved),
        "moved_bp": sum(length[r] for r in moved),
        "rearrangement_class": ("collinear" if not calls and not moved else
                                "inversions" if not moved else
                                "translocations" if not calls else
                                "inversions_and_translocations"),
    })
    return sorted(calls, key=lambda e: e["ref_start"]), info
