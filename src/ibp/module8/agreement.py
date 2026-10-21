"""Auditable, provisional Leen/SyRI agreement from Module 7's exact call IDs.

Agreement is an *intra-genome, candidate breakpoint association*, not biological
truth, independent evolutionary events, or evidence of absent inversions.
"""

from collections import Counter, defaultdict
from pathlib import Path

from ibp.module7.aggregate import OUTPUTS, read_table

AGREEMENT_FIELDS = [
    "species_id", "cohort_id", "joint_callable_query_genomes", "leen_calls",
    "syri_calls", "unique_candidate_pairs", "leen_only_calls", "syri_only_calls",
    "ambiguous_leen_calls", "ambiguous_syri_calls", "unassessable_leen_calls",
    "unassessable_syri_calls", "leen_families_with_unique_support",
    "syri_families_with_unique_support", "interpretation",
]
CALL_FIELDS = [
    "species_id", "cohort_id", "caller", "query_assembly", "call_id", "family_id",
    "reported_length_bp", "category", "other_caller_callable", "possible_partner_call_ids",
]
FAMILY_SUPPORT_FIELDS = [
    "species_id", "cohort_id", "caller", "family_id", "median_reported_length_bp",
    "member_calls", "unique_candidate_calls", "ambiguous_candidate_calls",
    "caller_only_calls", "unassessable_calls", "support_class", "interpretation",
]
SUPPORT_CLASSES = (
    ("has_unique_candidate", "Unique candidate"),
    ("ambiguous_only", "Ambiguous only"),
    ("caller_only_with_evaluable_calls", "No matched call"),
    ("not_evaluable", "Unassessable"),
)


def build_agreement(indir):
    """Return per-species summary, per-call classifications and per-family evidence.

    A unique edge is counted ONCE as a pair, though each caller also has a
    unique-matched call. A many-to-many edge is *not* interpreted as one event.
    Caller-only means unmatched in a genome where *both* comparisons succeeded.
    """
    indir = Path(indir)
    families = read_table(indir / "candidate_families.tsv", OUTPUTS["candidate_families.tsv"])
    mappings = read_table(indir / "call_family_map.tsv", OUTPUTS["call_family_map.tsv"])
    associations = read_table(indir / "cross_caller_associations.tsv", OUTPUTS["cross_caller_associations.tsv"])
    genomes = read_table(indir / "genome_sizes.tsv", OUTPUTS["genome_sizes.tsv"])
    summaries = read_table(indir / "species_summary.tsv", OUTPUTS["species_summary.tsv"])

    species_keys = {(r["species_id"], r["cohort_id"]) for r in summaries}
    species_summary = defaultdict(dict)
    for row in summaries:
        k = row["species_id"], row["cohort_id"]
        if row["caller"] in species_summary[k]:
            raise ValueError(f"Duplicate caller summary for {k}/{row['caller']}")
        species_summary[k][row["caller"]] = row
    if any(set(species_summary[k]) != {"leen", "syri"} for k in species_keys):
        raise ValueError("Each species must have one Leen and one SyRI summary")

    genome_index = {}
    for g in genomes:
        k = g["species_id"], g["cohort_id"], g["assembly_id"]
        if k in genome_index or k[:2] not in species_keys:
            raise ValueError(f"Duplicate or unknown genome: {k}")
        if any(g[x] not in {"true", "false"} for x in ("is_reference", "leen_callable", "syri_callable")):
            raise ValueError(f"Invalid callable flag: {k}")
        genome_index[k] = g

    family_index = {}
    for f in families:
        k = f["species_id"], f["cohort_id"], f["caller"], f["family_id"]
        if k in family_index or k[:2] not in species_keys or k[2] not in {"leen", "syri"}:
            raise ValueError(f"Duplicate/unknown family: {k}")
        family_index[k] = f

    call_index = {}
    by_family = defaultdict(list)
    for m in mappings:
        k = (m["species_id"], m["cohort_id"], m["caller"],
             m["query_assembly"], m["call_id"])
        family_key = (*k[:3], m["family_id"])
        genome_key = (*k[:2], k[3])
        if k in call_index or family_key not in family_index or genome_key not in genome_index:
            raise ValueError(f"Duplicate/orphan mapped call: {k}")
        if genome_index[genome_key]["is_reference"] == "true" or (
            genome_index[genome_key][f"{k[2]}_callable"] != "true"
        ):
            raise ValueError(f"Inversion call in non-callable/reference genome: {k}")
        if float(m["reported_length_bp"]) <= 0:
            raise ValueError(f"Invalid call length: {k}")
        call_index[k] = m
        by_family[family_key].append(k)

    edges = defaultdict(list)
    seen_edges = set()
    for edge in associations:
        s, c, query = (edge[x] for x in ("species_id", "cohort_id", "query_assembly"))
        lkey = (s, c, "leen", query, edge["leen_call_id"])
        skey = (s, c, "syri", query, edge["syri_call_id"])
        pair = lkey, skey
        if pair in seen_edges or lkey not in call_index or skey not in call_index:
            raise ValueError(f"Duplicate/orphan cross-caller edge: {pair}")
        seen_edges.add(pair)
        if (call_index[lkey]["family_id"] != edge["leen_family_id"] or
                call_index[skey]["family_id"] != edge["syri_family_id"]):
            raise ValueError(f"Cross-caller family IDs disagree: {pair}")
        if edge["association_status"] not in {"unique_candidate", "ambiguous_candidate"}:
            raise ValueError(f"Invalid cross-caller status: {pair}")
        edges[lkey].append((skey, edge["association_status"]))
        edges[skey].append((lkey, edge["association_status"]))
    for key, partners in edges.items():
        expected = ("unique_candidate" if len(partners) == 1
                    and len(edges[partners[0][0]]) == 1 else "ambiguous_candidate")
        if any(p[1] != expected for p in partners):
            raise ValueError(f"Non-reciprocal/incorrect candidate edge cardinality: {key}")

    call_rows = []
    counts = defaultdict(Counter)
    family_counts = defaultdict(Counter)
    for k, m in sorted(call_index.items()):
        s, c, caller, query, call_id = k
        other = "syri" if caller == "leen" else "leen"
        genome = genome_index[s, c, query]
        other_callable = genome[f"{other}_callable"] == "true"
        matched = edges.get(k, [])
        if matched and not other_callable:
            raise ValueError(f"Cross-caller pair with unsuccessful comparison: {k}")
        if matched:
            category = matched[0][1]
        else:
            category = "caller_only" if other_callable else "unassessable"
        fkey = (s, c, caller, m["family_id"])
        family_counts[fkey][category] += 1
        counts[s, c][f"{caller}_{category}"] += 1
        call_rows.append({
            "species_id": s, "cohort_id": c, "caller": caller,
            "query_assembly": query, "call_id": call_id, "family_id": m["family_id"],
            "reported_length_bp": m["reported_length_bp"], "category": category,
            "other_caller_callable": str(other_callable).lower(),
            "possible_partner_call_ids": ";".join(sorted(p[0][4] for p in matched)),
        })
    for (s, c), caller_rows in species_summary.items():
        for caller, row in caller_rows.items():
            if int(row["calls"]) != sum(k[:3] == (s, c, caller) for k in call_index):
                raise ValueError(f"Module 7 call total mismatch for {s}/{caller}")

    family_rows = []
    for k, family in sorted(family_index.items()):
        s, c, caller, fid = k
        n = family_counts[k]
        if int(family["member_calls"]) != len(by_family[k]):
            raise ValueError(f"Module 7 family member mismatch: {k}")
        if n["unique_candidate"]:
            support = "has_unique_candidate"
        elif n["ambiguous_candidate"]:
            support = "ambiguous_only"
        elif n["caller_only"]:
            support = "caller_only_with_evaluable_calls"
        else:
            support = "not_evaluable"
        family_rows.append({
            "species_id": s, "cohort_id": c, "caller": caller, "family_id": fid,
            "median_reported_length_bp": family["median_reported_length_bp"],
            "member_calls": family["member_calls"],
            "unique_candidate_calls": n["unique_candidate"],
            "ambiguous_candidate_calls": n["ambiguous_candidate"],
            "caller_only_calls": n["caller_only"],
            "unassessable_calls": n["unassessable"], "support_class": support,
            "interpretation": "provisional_within_genome_caller_support_not_truth",
        })

    unique_edges_by_species = Counter(
        (r["species_id"], r["cohort_id"]) for r in associations
        if r["association_status"] == "unique_candidate"
    )
    summary_rows = []
    for s, c in sorted(species_keys):
        observed = counts[s, c]
        matched = unique_edges_by_species[s, c]
        if observed["leen_unique_candidate"] != matched or observed["syri_unique_candidate"] != matched:
            raise ValueError(f"Unique match count mismatch for {s}")
        both_callable = sum(
            g["is_reference"] == "false" and g["leen_callable"] == "true"
            and g["syri_callable"] == "true"
            for k, g in genome_index.items() if k[:2] == (s, c)
        )
        summary_rows.append({
            "species_id": s, "cohort_id": c,
            "joint_callable_query_genomes": both_callable,
            "leen_calls": int(species_summary[s, c]["leen"]["calls"]),
            "syri_calls": int(species_summary[s, c]["syri"]["calls"]),
            "unique_candidate_pairs": matched,
            "leen_only_calls": observed["leen_caller_only"],
            "syri_only_calls": observed["syri_caller_only"],
            "ambiguous_leen_calls": observed["leen_ambiguous_candidate"],
            "ambiguous_syri_calls": observed["syri_ambiguous_candidate"],
            "unassessable_leen_calls": observed["leen_unassessable"],
            "unassessable_syri_calls": observed["syri_unassessable"],
            "leen_families_with_unique_support": sum(
                f["species_id"] == s and f["cohort_id"] == c and f["caller"] == "leen"
                and f["support_class"] == "has_unique_candidate" for f in family_rows
            ),
            "syri_families_with_unique_support": sum(
                f["species_id"] == s and f["cohort_id"] == c and f["caller"] == "syri"
                and f["support_class"] == "has_unique_candidate" for f in family_rows
            ),
            "interpretation": "unique_breakpoint_candidates_not_confirmed_consensus",
        })
    return summary_rows, call_rows, family_rows


def make_agreement_figures(summary_rows, family_rows, outdir):
    """Draw two new figures; include explicit annotations for sparse datasets."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = Path(outdir) / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    # A readable subset of species is shown in the PNG. The TSV includes all.
    groups = [
        ("leen_only_calls", "Leen-only calls"),
        ("syri_only_calls", "SyRI-only calls"),
        ("unique_candidate_pairs", "Unique pairs"),
        ("ambiguous_leen_calls", "Ambiguous Leen calls"),
        ("ambiguous_syri_calls", "Ambiguous SyRI calls"),
    ]
    ranked = sorted(summary_rows, key=lambda row: (
        -sum(int(row[field]) for field, _ in groups), row["species_id"]
    ))[:25]
    fig, ax = plt.subplots(figsize=(11, max(4, len(ranked) * 0.33 + 1.5)))
    if ranked:
        positions = list(range(len(ranked)))
        left = [0] * len(ranked)
        for field, label in groups:
            counts = [int(r[field]) for r in ranked]
            ax.barh(positions, counts, left=left, label=label)
            left = [x + y for x, y in zip(left, counts)]
        ax.set_yticks(positions, [r["species_id"] for r in ranked], fontsize=8)
        ax.invert_yaxis()
        fig.legend(loc="lower center", bbox_to_anchor=(0.5, 0.095), ncol=3, fontsize=8)
    else:
        ax.text(0.5, 0.5, "No species with valid catalogues", ha="center", va="center", transform=ax.transAxes)
    ax.set(xlabel="Unique candidate pairs / individual unmatched or ambiguous calls",
           title=f"Leen vs SyRI provisional agreement — top {len(ranked)} species")
    fig.text(0.5, 0.01, "Unassessable calls (other caller failed) are excluded; see caller_agreement.tsv."
             " Ambiguous Leen/SyRI counts are calls, not pair counts.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.23, 1, 0.96))
    fig.savefig(figures / "caller_agreement.png", dpi=175, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, caller in zip(axes, ("leen", "syri")):
        categories = [
            [float(f["median_reported_length_bp"]) / 1000 for f in family_rows
             if f["caller"] == caller and f["support_class"] == support]
            for support, _ in SUPPORT_CLASSES
        ]
        nonempty = [(i + 1, values) for i, values in enumerate(categories) if values]
        if nonempty:
            ax.boxplot([v for _, v in nonempty], positions=[i for i, _ in nonempty],
                       widths=0.55, showfliers=False)
            ax.set_yscale("log")
        else:
            ax.text(0.5, 0.5, "No candidate families", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim(0.5, len(SUPPORT_CLASSES) + 0.5)
        ax.set_xticks(list(range(1, len(SUPPORT_CLASSES) + 1)),
                      [f"{label}\n(n={len(values)})" for (_, label), values in zip(SUPPORT_CLASSES, categories)],
                      fontsize=8)
        ax.set_title(caller.upper())
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Family median caller-reported length (kb; log scale)")
    fig.suptitle("Candidate inversion-family lengths by provisional caller support")
    fig.text(0.5, 0.015, "A family is 'unique-supported' if at least one member call has a unique candidate match."
             " See caller_agreement_families.tsv for mixed support.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(figures / "agreement_family_lengths.png", dpi=175, bbox_inches="tight")
    plt.close(fig)
