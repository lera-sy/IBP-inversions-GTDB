"""Run tree-free inversion analysis across independently recoverable species cohorts.

The existing Modules 1-6 run unchanged inside a distinct cohort for each
species. Only fully validated Module 6 catalogues are passed to Module 7;
failures and partial caller coverage are recorded, never treated as zeros.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import yaml

from ibp.module5.results import COMPARISON_FIELDS
from ibp.module6.catalog import SUMMARY_FIELDS
from ibp.module7.aggregate import aggregate, read_table, write_table
from ibp.module8.analyze import analyze

STATUS_FIELDS = [
    "species_id", "taxid", "source_cohort", "status", "stage", "comparisons_success",
    "comparisons_unresolved", "comparisons_unsupported", "comparisons_failed",
    "leen_candidate_families", "syri_candidate_families", "log_file", "error",
]
SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def validate_name(value, label):
    if not isinstance(value, str) or not SAFE_NAME.fullmatch(value) or value in {".", ".."}:
        raise ValueError(f"{label} must use only letters, numbers, '.', '-' and '_' and not include paths")
    return value


def get_species(repo, configuration):
    """Use existing Module 1 selection; never silently widen the species set."""
    selected = pd.read_csv(repo / "data/processed/selected_assemblies.tsv", sep="\t", dtype=str)
    species = pd.read_csv(repo / "data/processed/species.tsv", sep="\t", dtype=str)
    if selected.empty or "Organism Taxonomic ID" not in selected.columns:
        raise ValueError("Missing or empty Module 1 species selection")
    if species[["taxid", "species_id"]].isna().any().any() or species["taxid"].duplicated().any():
        raise ValueError("Malformed/duplicate species IDs in data/processed/species.tsv")
    by_taxid = dict(zip(species["taxid"], species["species_id"]))
    by_id = dict(zip(species["species_id"], species["taxid"]))
    if len(by_id) != len(by_taxid):
        raise ValueError("Duplicate species IDs in data/processed/species.tsv")
    available = set(selected["Organism Taxonomic ID"].dropna())
    requested = configuration.get("species", sorted(available))
    if not isinstance(requested, list) or not requested:
        raise ValueError("species must be a nonempty list, or omitted to use all Module 1 species")
    if "species_file" in configuration:
        raise ValueError("Use a species list or omit it; species_file is not supported by this runner")
    result = []
    for name in requested:
        taxid = by_id.get(str(name), str(name))
        if taxid not in by_taxid:
            raise ValueError(f"TaxID {taxid} is missing from data/processed/species.tsv")
        if taxid not in available:
            raise ValueError(f"Requested species {name} is absent from Module 1 selected assemblies")
        species_id = validate_name(by_taxid[taxid], "species_id")
        result.append((species_id, taxid))
    if len({x[0] for x in result}) != len(result):
        raise ValueError("Repeated species in selection")
    return result


def paths_for(scratch, study, species_id):
    source_cohort = validate_name(f"{study}__{species_id}", "source_cohort")
    directory = scratch / "IBP-main" / "runs" / source_cohort
    summary = directory / "module6" / species_id / source_cohort / "summary.tsv"
    comparisons = directory / "module5" / species_id / source_cohort / "comparisons.tsv"
    return source_cohort, summary, comparisons


def read_comparisons(path):
    rows = read_table(path, COMPARISON_FIELDS)
    counts = {kind: sum(r["status"] == kind for r in rows)
              for kind in ("success", "unresolved", "unsupported", "failed")}
    if not rows or sum(counts.values()) != len(rows):
        raise ValueError(f"Invalid comparison status table: {path}")
    return counts


def row_template(species_id, taxid, source_cohort, log_file):
    return {
        "species_id": species_id, "taxid": taxid, "source_cohort": source_cohort,
        "status": "failed", "stage": "snakemake", "comparisons_success": "",
        "comparisons_unresolved": "", "comparisons_unsupported": "",
        "comparisons_failed": "", "leen_candidate_families": "",
        "syri_candidate_families": "", "log_file": str(log_file), "error": "",
    }


def validate_catalogue(summary, comparisons, source_cohort, species_id, validation_root):
    """Validate each candidate independently, so one damaged output is isolated."""
    rows = read_table(summary, SUMMARY_FIELDS)
    if len(rows) != 2 or {r["caller"] for r in rows} != {"leen", "syri"}:
        raise ValueError(f"Incomplete Module 6 summary: {summary}")
    if {(r["species_id"], r["cohort_id"]) for r in rows} != {(species_id, source_cohort)}:
        raise ValueError(f"Module 6 identity mismatch: {summary}")
    # Module 8 requires a nonzero callable denominator for every caller.
    # If a species has no callable genomes, report it as unanalyzable rather
    # than allowing it to invalidate every otherwise successful species.
    if any(int(r["callable_genomes"]) < 1 for r in rows):
        raise ValueError(f"No callable genomes for at least one caller: {summary}")
    counts = read_comparisons(comparisons)
    with tempfile.TemporaryDirectory(dir=validation_root) as tmp:
        # Reuse Module 7's complete schema, identity and referential checks.
        aggregate([summary], Path(tmp), source_cohort)
    return rows, counts


def run_one(repo, scratch, study_root, study, original_config, species_id, taxid,
            cores, dry_run=False, timeout_seconds=None):
    source_cohort, summary, comparisons = paths_for(scratch, study, species_id)
    config_dir = study_root / "configs"
    log_dir = study_root / "logs"
    config_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{species_id}.log"
    entry = row_template(species_id, taxid, source_cohort, log_path)
    individual = dict(original_config)
    individual.pop("species_file", None)
    individual.update({
        "cohort": source_cohort,
        "species": [species_id],
        "isolate_run": True,
        # Even the full study needs a one-species selector, not the raw
        # multi-species metadata. Zero retains all selected genomes for this species.
        "run_assembly_limit": int(original_config.get("run_assembly_limit", 0)),
    })
    config_path = config_dir / f"{species_id}.yaml"
    config_path.write_text(yaml.safe_dump(individual, sort_keys=False), encoding="utf-8")
    command = [
        "snakemake", "--snakefile", "Snakefile", "--configfile", str(config_path),
        "--cores", str(cores), "--rerun-incomplete", "module6_complete",
    ]
    if dry_run:
        command.insert(-1, "--dry-run")
    with log_path.open("w", encoding="utf-8") as output:
        output.write("Command: " + " ".join(command) + "\n\n")
        output.flush()
        try:
            result = subprocess.run(
                command, cwd=repo, env=os.environ.copy(), stdout=output,
                stderr=subprocess.STDOUT, check=False, timeout=timeout_seconds,
            )
            if result.returncode:
                raise RuntimeError(f"Snakemake exited with code {result.returncode}")
            if dry_run:
                entry.update(status="planned", stage="dry_run")
                return entry, None
            summaries, counts = validate_catalogue(
                summary, comparisons, source_cohort, species_id, study_root
            )
            entry.update(
                status="partial" if sum(counts[k] for k in ("unresolved", "unsupported", "failed")) else "completed",
                stage="module6_complete",
                **{f"comparisons_{k}": str(v) for k, v in counts.items()},
                **{f"{r['caller']}_candidate_families": r["candidate_families"] for r in summaries},
            )
            return entry, summary
        except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
            entry.update(status="failed", error=f"{type(error).__name__}: {error}")
            output.write("\nRunner error: " + entry["error"] + "\n")
            return entry, None


def run_study(repo, scratch, config, cores, species_workers=1, plan=False,
              dry_run=False, timeout_minutes=0):
    repo = Path(repo).resolve()
    scratch = Path(scratch).resolve()
    configuration = yaml.safe_load(Path(config).read_text(encoding="utf-8"))
    if not isinstance(configuration, dict):
        raise TypeError("Study YAML must contain a mapping")
    study = validate_name(configuration.get("cohort"), "cohort")
    if cores < 1 or species_workers < 1 or species_workers > cores:
        raise ValueError("cores >= species_workers >= 1 required")
    if int(configuration.get("run_assembly_limit", 0)) < 0:
        raise ValueError("run_assembly_limit cannot be negative")
    selection = get_species(repo, configuration)
    study_root = scratch / "IBP-studies" / study
    study_root.mkdir(parents=True, exist_ok=True)
    if plan:
        for species_id, taxid in selection:
            source, summary, _ = paths_for(scratch, study, species_id)
            print(f"{species_id} (TaxID {taxid}): {source} -> {summary}")
        print(f"Plan: {len(selection)} species, {cores} total cores, "
              f"{species_workers} concurrent species; no work executed")
        return []
    per_species_cores = cores // species_workers
    timeout_seconds = timeout_minutes * 60 if timeout_minutes else None
    statuses = []
    summaries = {}
    # Independent Snakemake executions and result trees. A failed species
    # cannot prevent the remaining species from running.
    with ThreadPoolExecutor(max_workers=species_workers) as pool:
        futures = {
            pool.submit(run_one, repo, scratch, study_root, study, configuration,
                        species_id, taxid, per_species_cores, dry_run, timeout_seconds): species_id
            for species_id, taxid in selection
        }
        for future in as_completed(futures):
            species_id = futures[future]
            try:
                status, summary = future.result()
            except Exception as error:  # noqa: BLE001 - per-species safety boundary
                # Includes unexpected Python errors; retain them in the manifest.
                source, _, _ = paths_for(scratch, study, species_id)
                taxid = next(t for s, t in selection if s == species_id)
                status = row_template(species_id, taxid, source,
                                      study_root / "logs" / f"{species_id}.log")
                status["error"] = f"Unexpected runner failure: {type(error).__name__}: {error}"
                summary = None
            statuses.append(status)
            if summary is not None:
                summaries[species_id] = summary
            print(f"{species_id}: {status['status']} "
                  f"({status['error'] or 'see per-species log'})", flush=True)
            write_table(study_root / "species_status.tsv", STATUS_FIELDS,
                        sorted(statuses, key=lambda r: r["species_id"]))
    if dry_run:
        print(f"Dry-run complete for {len(statuses)} species. "
              f"See {study_root / 'species_status.tsv'}")
        return statuses
    if not summaries:
        print("No species produced a validated catalogue; aggregation was not run.", file=sys.stderr)
        return statuses
    # Publish only successful/partial species. Module 7 retains its normal
    # fail-closed validation, but checks original per-species cohort identity.
    ordered = sorted(summaries)
    cohorts = {s: paths_for(scratch, study, s)[0] for s in ordered}
    out7 = study_root / "module7" / study
    out8 = study_root / "module8" / study
    aggregate([summaries[s] for s in ordered], out7, study, source_cohorts=cohorts)
    analyze(out7, out8)
    print(f"Study finished: {len(summaries)}/{len(selection)} species with validated "
          f"catalogues. Status report: {study_root / 'species_status.tsv'}")
    print(f"Figures: {out8 / 'figures'}")
    return statuses


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--cores", type=int, default=4, help="Total CPU slots for all species")
    parser.add_argument("--species-workers", type=int, default=1,
                        help="Concurrent isolated species workflows; default 1")
    parser.add_argument("--timeout-minutes", type=int, default=0,
                        help="Optional maximum time per species; 0 disables timeout")
    parser.add_argument("--plan", action="store_true", help="Print species/path plan only")
    parser.add_argument("--dry-run", action="store_true",
                        help="Run each per-species Snakemake DAG in dry-run mode")
    args = parser.parse_args(argv)
    if args.timeout_minutes < 0:
        parser.error("--timeout-minutes must be >= 0")
    scratch = os.environ.get("VSC_SCRATCH")
    if not scratch:
        parser.error("Set VSC_SCRATCH to a writable scratch directory")
    try:
        statuses = run_study(Path.cwd(), Path(scratch), args.configfile,
                             args.cores, args.species_workers, args.plan,
                             args.dry_run, args.timeout_minutes)
    except (TypeError, ValueError, OSError, KeyError) as error:
        parser.error(str(error))
    if statuses and not args.dry_run and all(s["status"] == "failed" for s in statuses):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
