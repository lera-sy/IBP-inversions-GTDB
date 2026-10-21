"""Tree-free, caller-separated inversion family catalogue.

Reuses legacy complete-link breakpoint grouping and explicit 0/1/? states.
Cross-caller associations are *candidates*, never confirmed consensus.
Input coordinates are Module 5's 0-based, half-open reference coordinates.
No phylogenetic event inference is performed here.
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

from ibp.module5.results import COMPARISON_FIELDS, FIELDS
from ibp.module5.validate_input import load_package

FAMILY_FIELDS = [
    "species_id", "cohort_id", "caller", "family_id", "reference_assembly",
    "ref_start_0", "ref_end_0", "median_reported_length_bp", "observed_carriers",
    "member_calls", "duplicate_calls_within_genomes", "geometry_status",
    "interpretation",
]
MAP_FIELDS = [
    "caller", "family_id", "query_assembly", "call_id", "ref_intervals_0",
    "reported_length_bp", "geometry_status",
]
STATE_FIELDS = ["caller", "family_id", "assembly_id", "state"]
ASSOCIATION_FIELDS = [
    "query_assembly", "leen_call_id", "syri_call_id", "leen_family_id",
    "syri_family_id", "left_breakpoint_distance_bp", "right_breakpoint_distance_bp",
    "association_status", "interpretation",
]
SUMMARY_FIELDS = [
    "species_id", "cohort_id", "caller", "total_genomes", "callable_genomes",
    "unknown_genomes", "candidate_families", "families_with_multiple_carriers",
    "calls", "complex_geometry_calls", "candidate_cross_caller_pairs",
    "unique_cross_caller_pairs", "ambiguous_cross_caller_pairs", "interpretation",
]
SIZE_FIELDS = [
    "species_id", "cohort_id", "assembly_id", "chromosome_length_bp", "genome_length_bp",
    "is_reference", "leen_callable", "syri_callable",
]


def read_table(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_table(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def distance(a, b, chromosome_length, circular):
    delta = abs(a - b)
    return min(delta, chromosome_length - delta) if circular else delta


def endpoints(raw, chromosome_length, circular):
    """Return an unambiguous breakpoint pair, or None for complex geometry.

    Module 5 may emit multiple Leen segments. Only a single interval or two
    pieces meeting at the reference origin can safely give one breakpoint pair.
    """
    parts = []
    for entry in raw.split(";"):
        try:
            a, b = map(int, entry.split("-"))
        except (ValueError, AttributeError):
            raise ValueError(f"Invalid standardized interval: {raw!r}") from None
        if not (0 <= a < chromosome_length and 0 <= b <= chromosome_length):
            raise ValueError(f"Interval outside reference: {raw!r}")
        if a == b or (a > b and not circular):
            raise ValueError(f"Invalid interval: {raw!r}")
        parts.append((a, b))
    if len(parts) == 1:
        return parts[0]
    if circular and len(parts) == 2:
        final = [p for p in parts if p[1] == chromosome_length]
        initial = [p for p in parts if p[0] == 0]
        if len(final) == 1 and len(initial) == 1 and final[0] != initial[0]:
            return final[0][0], initial[0][1]
    return None


def compatible(a, b, length, circular, tolerance):
    if a["endpoints"] is None or b["endpoints"] is None:
        return False
    return all(distance(x, y, length, circular) <= tolerance
               for x, y in zip(a["endpoints"], b["endpoints"]))


def families_for_caller(calls, caller, species, cohort, reference, genomes, callable_ids,
                        length, circular, tolerance):
    """Greedy deterministic complete-link grouping, as in legacy SyRI code."""
    groups = []
    for call in sorted(calls, key=lambda c: (
        c["endpoints"] is None, c["endpoints"] or (length, length),
        c["query_assembly"], c["call_id"],
    )):
        eligible = [g for g in groups if all(
            compatible(call, member, length, circular, tolerance) for member in g
        )]
        if eligible:
            chosen = min(eligible, key=lambda g: (
                sum(sum(distance(x, y, length, circular)
                        for x, y in zip(call["endpoints"], member["endpoints"]))
                    for member in g) / len(g),
                groups.index(g),
            ))
            chosen.append(call)
        else:
            groups.append([call])
    rows, mapping, states = [], [], []
    for index, group in enumerate(groups, 1):
        fid = f"{'L' if caller == 'leen' else 'S'}{index:04d}"
        carriers = {c["query_assembly"] for c in group}
        count = Counter(c["query_assembly"] for c in group)
        geometry = "simple" if group[0]["endpoints"] is not None else "complex_unclustered"
        first, second = group[0]["endpoints"] if geometry == "simple" else ("", "")
        rows.append({
            "species_id": species, "cohort_id": cohort, "caller": caller,
            "family_id": fid, "reference_assembly": reference,
            "ref_start_0": first, "ref_end_0": second,
            "median_reported_length_bp": median(c["reported_length_bp"] for c in group),
            "observed_carriers": len(carriers), "member_calls": len(group),
            "duplicate_calls_within_genomes": sum(v - 1 for v in count.values()),
            "geometry_status": geometry,
            "interpretation": "candidate_family_not_independent_event",
        })
        for call in group:
            call["family_id"] = fid
            mapping.append({k: call[k] for k in (
                "caller", "query_assembly", "call_id", "ref_intervals_0",
                "reported_length_bp",
            )} | {"family_id": fid, "geometry_status": geometry})
        for accession in sorted(genomes):
            states.append({"caller": caller, "family_id": fid,
                           "assembly_id": accession, "state":
                           "?" if accession not in callable_ids else
                           "1" if accession in carriers else "0"})
    return rows, mapping, states


def build_catalogue(manifest, report, outdir, tolerance=5000):
    if tolerance < 0:
        raise ValueError("Breakpoint tolerance must be nonnegative")
    package = load_package(manifest)
    species, cohort = package.species_id, package.cohort_id
    reference = package.reference.assembly_id
    genomes = {g.assembly_id: g for g in (package.reference, *package.queries)}
    ref_length = package.reference.chromosome_length_bp
    circular = package.reference.topology == "circular"
    query_ids = set(genomes) - {reference}
    fields, comparisons = read_table(report)
    if fields != COMPARISON_FIELDS:
        raise ValueError("Module 5 comparison report has unexpected schema")
    indexed = {}
    by_caller = defaultdict(list)
    seen_call_ids = set()
    statuses = {"leen": {}, "syri": {}}
    for row in comparisons:
        caller, query = row["caller"], row["query_assembly"]
        if (row["species_id"], row["cohort_id"], row["reference_assembly"]) != (
            species, cohort, reference
        ) or caller not in statuses or query not in query_ids:
            raise ValueError("Comparison species/cohort/reference/query mismatch")
        key = caller, query
        if key in indexed:
            raise ValueError(f"Duplicate Module 5 comparison: {key}")
        indexed[key] = row
        outcome = row["status"]
        if outcome not in {"success", "unresolved", "unsupported", "failed"}:
            raise ValueError(f"Unknown comparison status: {outcome}")
        statuses[caller][query] = outcome
        columns, calls = read_table(row["inversions_path"])
        if columns != FIELDS:
            raise ValueError("Module 5 inversion table has unexpected schema")
        if outcome == "success":
            if row["call_count"] != str(len(calls)):
                raise ValueError(f"Call-count mismatch for {key}")
        elif calls or row["call_count"]:
            raise ValueError(f"Unsuccessful comparison must not have calls: {key}")
        for call in calls:
            if (call["caller"], call["reference_assembly"], call["query_assembly"]) != (
                caller, reference, query
            ):
                raise ValueError(f"Inversion identity mismatch: {key}")
            value = int(call["reported_length_bp"])
            if not 0 < value <= ref_length:
                raise ValueError("Invalid reported inversion length")
            call_key = (caller, query, call["call_id"])
            if not call["call_id"] or call_key in seen_call_ids:
                raise ValueError(f"Duplicate or empty inversion call ID: {call_key}")
            seen_call_ids.add(call_key)
            parsed = endpoints(call["ref_intervals_0"], ref_length, circular)
            by_caller[caller].append({
                "caller": caller, "query_assembly": query, "call_id": call["call_id"],
                "ref_intervals_0": call["ref_intervals_0"], "reported_length_bp": value,
                "endpoints": parsed,
            })
    expected = {(caller, query) for caller in statuses for query in query_ids}
    if set(indexed) != expected:
        raise ValueError(f"Missing Module 5 comparisons: {sorted(expected - set(indexed))}")

    family_rows, map_rows, state_rows = [], [], []
    for caller in ("leen", "syri"):
        available = {reference} | {g for g, v in statuses[caller].items() if v == "success"}
        families, mapping, states = families_for_caller(
            by_caller[caller], caller, species, cohort, reference, genomes,
            available, ref_length, circular, tolerance,
        )
        family_rows.extend(families)
        map_rows.extend(mapping)
        state_rows.extend(states)

    # Report all plausible same-genome pairings. Uniquely reciprocal candidates
    # are labelled provisional; one-to-many situations remain ambiguous.
    pairing_rows = []
    for query in sorted(query_ids):
        leen = [c for c in by_caller["leen"] if c["query_assembly"] == query]
        syri = [c for c in by_caller["syri"] if c["query_assembly"] == query]
        candidate_edges = [(lc, s) for lc in leen for s in syri
                           if compatible(lc, s, ref_length, circular, tolerance)]
        ld = Counter(id(lc) for lc, _ in candidate_edges)
        sd = Counter(id(s) for _, s in candidate_edges)
        for lc, s in candidate_edges:
            status = "unique_candidate" if ld[id(lc)] == sd[id(s)] == 1 else "ambiguous_candidate"
            pairing_rows.append({
                "query_assembly": query, "leen_call_id": lc["call_id"],
                "syri_call_id": s["call_id"], "leen_family_id": lc["family_id"],
                "syri_family_id": s["family_id"],
                "left_breakpoint_distance_bp": distance(
                    lc["endpoints"][0], s["endpoints"][0], ref_length, circular),
                "right_breakpoint_distance_bp": distance(
                    lc["endpoints"][1], s["endpoints"][1], ref_length, circular),
                "association_status": status,
                "interpretation": "provisional_breakpoint_proximity_not_consensus",
            })

    num_unique = sum(x["association_status"] == "unique_candidate" for x in pairing_rows)
    num_ambiguous = len(pairing_rows) - num_unique
    summary_rows = []
    for caller in ("leen", "syri"):
        success = sum(value == "success" for value in statuses[caller].values()) + 1
        f = [row for row in family_rows if row["caller"] == caller]
        summary_rows.append({
            "species_id": species, "cohort_id": cohort, "caller": caller,
            "total_genomes": len(genomes), "callable_genomes": success,
            "unknown_genomes": len(genomes) - success,
            "candidate_families": len(f),
            "families_with_multiple_carriers": sum(int(x["observed_carriers"]) > 1 for x in f),
            "calls": len(by_caller[caller]),
            "complex_geometry_calls": sum(x["endpoints"] is None for x in by_caller[caller]),
            "candidate_cross_caller_pairs": len(pairing_rows),
            "unique_cross_caller_pairs": num_unique,
            "ambiguous_cross_caller_pairs": num_ambiguous,
            "interpretation": "tree_free_candidate_families_not_independent_events",
        })
    size_rows = [{
        "species_id": species, "cohort_id": cohort, "assembly_id": assembly,
        "chromosome_length_bp": genome.chromosome_length_bp,
        "genome_length_bp": genome.genome_length_bp,
        "is_reference": str(assembly == reference).lower(),
        "leen_callable": str(assembly == reference or statuses["leen"][assembly] == "success").lower(),
        "syri_callable": str(assembly == reference or statuses["syri"][assembly] == "success").lower(),
    } for assembly, genome in sorted(genomes.items())]

    output = Path(outdir)
    write_table(output / "candidate_families.tsv", FAMILY_FIELDS, family_rows)
    write_table(output / "call_family_map.tsv", MAP_FIELDS, map_rows)
    write_table(output / "tip_states.tsv", STATE_FIELDS, state_rows)
    write_table(output / "cross_caller_associations.tsv", ASSOCIATION_FIELDS, pairing_rows)
    write_table(output / "genome_sizes.tsv", SIZE_FIELDS, size_rows)
    write_table(output / "summary.tsv", SUMMARY_FIELDS, summary_rows)
    print(f"Module 6 {species}/{cohort}: {len(family_rows)} candidate families; "
          f"{len(pairing_rows)} provisional cross-caller pairs "
          f"({num_unique} uniquely associated, {num_ambiguous} ambiguous)")
    return summary_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--breakpoint-tolerance", type=int, default=5000)
    args = parser.parse_args()
    build_catalogue(args.manifest, args.report, args.outdir, args.breakpoint_tolerance)


if __name__ == "__main__":
    main()
