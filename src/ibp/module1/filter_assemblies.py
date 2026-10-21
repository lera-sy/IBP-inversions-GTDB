"""
Filter the prokaryotic metadata so we only retain long-read and de novo assemblies.

Approach:
First we filter for long reads: sequencing technology has to match known long-read keywords, like "PacBio" or "Oxford Nanopore".
Next, de novo filtering: assembly method must match with known de novo assemblers (we use a whitelist for this)

Assemblies with unknown methods are kept separately for manual inspection (useful depending on how many assemblies are filtered out).
Assemblies with empty method fields are dropped as we cannot determine if they are de novo or not.

Note: Both these filters are empirical, thus not perfect. We may miss some assemblies that are long-read or de novo, and we may include some that are not --> important to manually verify depending on the number of assemblies filtered out.
"""

import argparse

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-c", "--candidates", required=True)
parser.add_argument("-u", "--unknown", required=True)
parser.add_argument("-r", "--reference-guided", required=True)
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t", dtype=str)
print(f"Loaded {len(df)} assemblies from the prokaryotic metadata.")


## Filter conditions

# Define col names
tech = "Assembly Sequencing Tech"
method = "Assembly Assembly Method"

# Long read filter (keywords)
LONGREAD_KEYWORDS = (
    r"pacbio|nanopore|\bont\b|smrt|hifi|hi-fi|\bclr\b"
    r"|minion|promethion|sequel|revio|gridion"
)

#only retain long-read assemblies and keep the hybrid ones identifiable downstream
technology = df[tech].fillna("")
has_long_read = technology.str.contains(LONGREAD_KEYWORDS, case=False, regex=True)
has_short_read = technology.str.contains(r"illumina|bgi|mgiseq|dnbseq|ion torrent|sanger|\b454\b|roche|aviti", case=False, regex=True)
df = df[has_long_read].copy()
df["Sequencing Category"] = "long-read-only"
hybrid_rows = has_short_read.reindex(df.index, fill_value=False)
df.loc[hybrid_rows, "Sequencing Category"] = "hybrid"
print(f"Long-read: {len(df)} assemblies retained after filtering for long-read sequencing technologies.")
print(df["Sequencing Category"].value_counts().to_string())

# Known de novo assemblers --> whitelist
DENOVO_ASSEMBLERS = (
    r"flye|canu|hgap|wtdbg|miniasm|raven"
    r"|hifiasm|peregrine|shasta|nextdenovo|necat"
    r"|falcon|unicycler"
    # Added 2026-10-21 (rerun): long-read de novo tools found among the
    # "unknown" methods. Autocycler/Trycycler build consensus from several
    # de novo assemblies; SMRT Analysis/SMRT Link "Microbial Assembly" is
    # PacBio's HGAP-based de novo assembler.
    r"|autocycler|trycycler|dragonflye|hybracter"
    r"|smrt[\s_-]*(?:analysis|link|portal)|microbial[\s_-]*assembl|pb[\s_-]*assembly"
)

# Reference guided / scaffolding keywords --> blacklist
REFERENCE_BLACKLIST = (
    r"ragtag|ragoo|saute|medusa|mauve"
    r"|reference[\s_-]*(?:guided|based|assisted|assembl|mapper|scaffold)"
    r"|assembled against (?:a )?reference|mapped to (?:a )?reference"
)

method_filled = df[method].fillna("").str.strip()
is_known = method_filled.str.len() > 0
is_ref_guided = method_filled.str.contains(REFERENCE_BLACKLIST, case=False, regex=True)
is_denovo = method_filled.str.contains(DENOVO_ASSEMBLERS, case=False, regex=True) & ~is_ref_guided

# Print some stats
candidates = df[is_denovo]
ref_guided = df[is_ref_guided]
unknown = df[is_known & ~is_denovo & ~is_ref_guided]
print(f"De novo (whitelisted assembler): {len(candidates):,}")
print(f"Reference guided (blacklisted): {len(ref_guided):,}")
print(f"Unknown method (needs manual check): {len(unknown):,}")
print(f"Dropped (empty method field): {len(df) - len(candidates) - len(ref_guided) - len(unknown):,}")

# Write out the filtered assemblies

candidates.to_csv(args.candidates, sep="\t", index=False)
unknown.to_csv(args.unknown, sep="\t", index=False)
ref_guided.to_csv(args.reference_guided, sep="\t", index=False)

print("\nWritten:")
print(f"  {args.candidates}       ({len(candidates):,} rows)")
print(f"  {args.unknown}  ({len(unknown):,} rows)")
print (f"  {args.reference_guided} ({len(ref_guided):,} rows)")
if len(unknown):
    print("\nTop unrecognized assembly-method values (check these):")
    print(unknown[method].value_counts().head(40).to_string())