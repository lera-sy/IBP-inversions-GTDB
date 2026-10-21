"""Shared, caller-independent inversion output schema.

Reference coordinates use 0-based, half-open intervals.
Native caller calls.tsv files remain unchanged.
"""

import csv
import json
import sys
from pathlib import Path

FIELDS = [
    "caller", "reference_assembly", "query_assembly", "call_id",
    "ref_start_0", "ref_end_0", "ref_intervals_0",
    "qry_start_0", "qry_end_0", "reported_length_bp",
    "interval_span_bp", "call_type", "spans_reference_origin",
    "coordinate_method",
]


def _check_interval(start, end, length, circular):
    if not (0 <= start < length and 0 <= end <= length):
        raise ValueError(f"Interval outside chromosome: {start}-{end}")
    if start == end or (start > end and not circular):
        raise ValueError(f"Invalid interval: {start}-{end}")
    return end - start if end > start else length - start + end


def standardize_calls(status, native_calls):
    """Convert successful Leen/SyRI calls without changing lengths."""
    if status["status"] != "success":
        raise ValueError("Only successful comparisons can be standardized")

    caller = status["caller"]
    if caller not in {"leen", "syri"}:
        raise ValueError(f"Unknown caller: {caller}")
    if status["call_count"] != len(native_calls):
        raise ValueError("Native call count disagrees with status.json")

    ref_length = int(status["reference_length_bp"])
    qry_length = int(status["query_length_bp"])
    normalized = []

    for index, call in enumerate(native_calls, 1):
        if caller == "leen":
            intervals = []
            for part in call["ref_intervals"].split(";"):
                start, end = map(int, part.split("-"))
                intervals.append((start, end))

            qry_start = qry_end = ""
            call_id = f"leen_{index:06d}"
            call_type = call["call_type"]
            method = "window_alignment_estimate"
            circular = True

        else:
            start = int(call["ref_start"]) - 1
            end = int(call["ref_end"])
            intervals = [(start, end)]

            qry_start = min(
                int(call["qry_start"]), int(call["qry_end"])
            ) - 1
            qry_end = max(
                int(call["qry_start"]), int(call["qry_end"])
            )

            if not (0 <= qry_start < qry_end <= qry_length):
                raise ValueError("SyRI query interval outside chromosome")

            call_id = call["syri_id"]
            call_type = "INV"
            method = "syri_native"
            circular = False

        if not intervals:
            raise ValueError("An inversion must have reference intervals")

        span = sum(
            _check_interval(start, end, ref_length, circular)
            for start, end in intervals
        )

        length = int(call["length"])
        if length <= 0 or length > ref_length:
            raise ValueError(f"Invalid caller-reported length: {length}")

        one_interval = len(intervals) == 1

        normalized.append({
            "caller": caller,
            "reference_assembly": status["reference_assembly"],
            "query_assembly": status["query_assembly"],
            "call_id": call_id,
            "ref_start_0": intervals[0][0] if one_interval else "",
            "ref_end_0": intervals[0][1] if one_interval else "",
            "ref_intervals_0": ";".join(
                f"{a}-{b}" for a, b in intervals
            ),
            "qry_start_0": qry_start,
            "qry_end_0": qry_end,
            "reported_length_bp": length,
            "interval_span_bp": span,
            "call_type": call_type,
            "spans_reference_origin": str(
                any(a > b for a, b in intervals)
            ).lower(),
            "coordinate_method": method,
        })

    return normalized


def write_standardized(path, status, native_calls):
    """Write the canonical inversion table after validation."""
    if status["status"] in {"unresolved", "unsupported", "failed"}:
        rows = []  # Placeholder only; see status.json
    else:
        rows = standardize_calls(status, native_calls)
    path = Path(path)
    temp = path.with_name(path.name + ".tmp")

    with temp.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    temp.replace(path)


COMPARISON_FIELDS = [
    "species_id", "cohort_id", "caller",
    "reference_assembly", "query_assembly",
    "status", "call_count", "reason",
    "inversions_path", "status_path",
]


def write_comparison_report(path, comparisons):
    """Audit completed comparisons and write one cohort index.

    Each comparison is:
    (species, cohort, caller, reference_id, query_id, output_dir).
    """
    rows = []

    for species, cohort, caller, reference, query, folder in comparisons:
        folder = Path(folder)
        status_path = folder / "status.json"
        inversions_path = folder / "inversions.tsv"

        status = json.loads(status_path.read_text())

        for field, expected in (
            ("caller", caller),
            ("reference_assembly", reference),
            ("query_assembly", query),
        ):
            if status.get(field) != expected:
                raise ValueError(
                    f"Comparison identity mismatch: {status_path}: {field}"
                )

        outcome = status.get("status")
        if outcome not in {
            "success", "unresolved", "unsupported", "failed"
        }:
            raise ValueError(f"Unknown comparison status: {outcome}")

        with inversions_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != FIELDS:
                raise ValueError(
                    f"Incorrect standardized columns: {inversions_path}"
                )
            calls = list(reader)

        if outcome == "success":
            count = status.get("call_count")
            if type(count) is not int or count < 0 or len(calls) != count:
                raise ValueError(
                    f"Inversion count mismatch: {status_path}"
                )

            if any(
                call["caller"] != caller
                or call["reference_assembly"] != reference
                or call["query_assembly"] != query
                for call in calls
            ):
                raise ValueError(
                    f"Inversion identity mismatch: {inversions_path}"
                )
        else:
            if calls:
                raise ValueError(
                    f"Unsuccessful comparison has calls: {inversions_path}"
                )
            count = ""

        rows.append({
            "species_id": species,
            "cohort_id": cohort,
            "caller": caller,
            "reference_assembly": reference,
            "query_assembly": query,
            "status": outcome,
            "call_count": count,
            "reason": status.get("reason") or status.get("error") or "",
            "inversions_path": str(inversions_path),
            "status_path": str(status_path),
        })

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")

    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COMPARISON_FIELDS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    temporary.replace(path)

    outcomes = ("success", "unresolved", "unsupported", "failed")
    counts = {
        outcome: sum(row["status"] == outcome for row in rows)
        for outcome in outcomes
    }

    summary = ", ".join(
        f"{outcome}={counts[outcome]}" for outcome in outcomes
    )
    print(f"Module 5 report: {path} | {summary}")

    unsuccessful = len(rows) - counts["success"]
    if unsuccessful:
        print(
            f"WARNING: {unsuccessful} comparisons were not successful. "
            f"Review {path} before downstream analysis.",
            file=sys.stderr,
        )
