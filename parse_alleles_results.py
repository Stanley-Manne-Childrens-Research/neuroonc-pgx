"""Aggregate StellarPGx CYP2D6 reports into calls, frequencies, and Table S1.

Requires Python 3.10 or later and the local sample_star_allele_table.py helper;
uses only the Python standard library. Keep the helper at the repository root
when invoking the implementation under Cyp2d6-Processing/ directly.
Example from the repository root:

    python3 aggregate_stellarpgx_calls.py --input-dir StellarPGxCalls

Upstream files and paths
-----------------------
The input directory must already contain text reports produced by the upstream
StellarPGx calling workflow. This script does not call genotypes, read BAM/VCF
files, or run StellarPGx. The default directory is "StellarPGxCalls" relative to the
working directory; --input-dir selects another location. In this workspace,
StellarPGxCalls is a symlink to Cyp2d6-Processing/alleles and can be used instead.
The chosen directory must exist and contain files matching --pattern (default:
"*.alleles"). Search is nonrecursive unless --recursive is supplied.

Report basenames identify samples. The full suffix
".filtered_ctrl_cyp2d6.alleles" is removed when present; otherwise only a final
".alleles" is removed. Directory names do not contribute to sample IDs. Other
prefixes/suffixes are not normalized, so IDs must already match any allowlist.

--sample-allowlist defaults to sample-cohort-manifest.txt relative to the working
directory. This pre-existing text file must contain one sample ID per line; use
--sample-allowlist PATH to select a different manifest. Blank lines are ignored
and duplicate allowlist entries collapse to one ID. The script does not create
a cohort manifest, inspect metadata, or verify that every
listed sample has a report. The allowlist affects frequencies and Table S1: the
detailed calls CSV includes every parsed report, including disallowed samples.

Report parsing and counting assumptions
--------------------------------------
Result, Activity score, and Metabolizer/Metaboliser status are recognized using
the case-insensitive line patterns in LABEL_PATTERNS. The first matching line
supplies either its inline value or the next nonblank line. That fallback does
not verify that the next line is a value rather than another heading. Absent
values are exported as empty strings. Activity scores and metabolizer statuses
are copied as text, not recalculated or checked for consistency. Undecodable
characters are replaced during file reading rather than raising a decode error.

A result is considered parseable for frequency counting if it contains a slash
and at least one nonempty slash-delimited token contains "*". Every nonempty
token then contributes one haplotype occurrence. This is a permissive text rule:
it does not require exactly two tokens, validate each token as a star allele, or
resolve ambiguous calls. Input Result fields must therefore already represent
the intended calls. Labels are retained verbatim after trimming whitespace;
copy-number labels such as *1x2 remain a single token, without expanding copies
or grouping labels into workbook categories.

Input reports must already be deduplicated. Each matching file is processed
independently, even if multiple files have the same sample ID. Haplotype count
includes repeated occurrences within a report (e.g., *1/*1 contributes two).
sample_count increments once per distinct haplotype per included report; it is
a count of unique individuals only when there is one report per sample.
Frequency = haplotype count / total parsed haplotype occurrences across included
reports. Missing or unparseable results contribute nothing to this denominator;
the denominator is not twice the allowlist size. No explicit zero-frequency
alleles are emitted.

Outputs and failure behavior
----------------------------
--output defaults to alleles_results.csv, with columns sample_id, file_name,
file_path (absolute), result, activity_score, metabolizer_status, and haplotypes
(tokens joined by "|"). There is one row per matched file, in sorted path order.
--haplotype-output defaults to haplotype_frequencies.csv, with columns haplotype,
count, sample_count, and frequency. Rows are ordered by decreasing count and
frequency is rounded to six decimal places. Both CSVs have headers; parent
directories are created and existing outputs are overwritten. Choose distinct
output paths. If no included result is parseable, the frequency CSV has only its
header. Missing input directories, no matching reports, or a missing requested
allowlist raise errors; other file-access failures also propagate.

--sample-star-allele-output defaults to table-s1-sample_star_allele.csv. It uses
the same columns as get_cbtn_allele_frequencies.py: sample_id, gene,
source_diplotype_count, counted_alleles, and labels_preview. Each cohort-eligible
report contributes a CYP2D6 row, including missing/unparseable results. A parsed
result contributes one source diplotype and its number of slash-delimited allele
tokens; an unparseable result contributes zero for both counts. labels_preview
retains the extracted Result text, including unparseable text, or is blank when
no result is found. Duplicate sample reports remain separate, as in the counts.
Existing CYP2D6 rows are replaced; rows for other genes are preserved. The
PharmCAT writer likewise preserves CYP2D6 rows, so either script can run first.
An existing table must have the exact shared header or the writer raises an
error rather than overwriting it. A valid run with no cohort-eligible reports
removes previous CYP2D6 rows. The summary includes the number of CYP2D6 rows.

Printed "Missing result" counts all parsed files; "Unresolved result" counts
only frequency-eligible reports that fail the haplotype parsing rule. Neither
figure is a count of unique individuals unless the inputs are deduplicated.

For downstream workbook updates, explicitly set --haplotype-output to
Cyp2d6-Processing/haplotype_frequencies.csv, the default input expected by
Cyp2d6-Processing/update_cyp2d6_frequency_table.py, or pass that updater --csv.
Likewise, --output Cyp2d6-Processing/alleles_results.csv places the calls where
build_shared_pharmcat_cyp2d6_samples.py expects them. Default output paths are
relative to the working directory, not to this script or its input directory.
This script does not update Excel workbooks or invoke either downstream script.
"""

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

# The root aggregation entry point may be a symlink to the legacy parser under
# Cyp2d6-Processing/. Make the shared root helper available for either invocation.
if Path(__file__).resolve().parent.name == "Cyp2d6-Processing":
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sample_star_allele_table import write_sample_star_allele_table


# Internal output-field names mapped to case-insensitive report-line patterns.
# Patterns are matched against stripped lines from the beginning of the line;
# the colon is optional and capture group 1 contains any inline value. Both
# British (metaboliser) and American (metabolizer) spellings are accepted.
LABEL_PATTERNS = {
    "result": re.compile(r"^result\s*:?\s*(.*)$", re.IGNORECASE),
    "activity_score": re.compile(r"^activity score\s*:?\s*(.*)$", re.IGNORECASE),
    "metabolizer_status": re.compile(
        r"^metaboli[sz]er status\s*:?\s*(.*)$", re.IGNORECASE
    ),
}

# Expected StellarPGx CYP2D6 report suffix, removed in full to recover the sample
# ID. Other report names fall back to removing only a trailing ".alleles".
COMMON_SUFFIX = ".filtered_ctrl_cyp2d6.alleles"

# Default cohort manifest, resolved relative to the working directory. It must
# exist before aggregation; --sample-allowlist can select another manifest.
DEFAULT_SAMPLE_ALLOWLIST = "sample-cohort-manifest.txt"

# Shared supplementary table. Replace CYP2D6 rows on each run and preserve
# PharmCAT rows; --sample-star-allele-output can select a different destination.
DEFAULT_SAMPLE_STAR_ALLELE_OUTPUT = "table-s1-sample_star_allele.csv"


def sample_id_from_name(file_name: str) -> str:
    """Derive an ID from a basename without normalizing other naming aliases."""
    if file_name.endswith(COMMON_SUFFIX):
        return file_name[: -len(COMMON_SUFFIX)]
    if file_name.endswith(".alleles"):
        return file_name[: -len(".alleles")]
    return file_name


def extract_labeled_value(text: str, field_name: str) -> str | None:
    """Read the first matching field inline or from the following nonblank line.

    field_name must be a key in LABEL_PATTERNS. Return None if no matching label
    or following value is found; fallback lines are not validated as field values.
    """
    pattern = LABEL_PATTERNS[field_name]
    lines = text.splitlines()

    for i, line in enumerate(lines):
        match = pattern.match(line.strip())
        if not match:
            continue

        inline_value = match.group(1).strip()
        if inline_value:
            return inline_value

        for candidate in lines[i + 1 :]:
            value = candidate.strip()
            if value:
                return value
        return None

    return None


def extract_haplotypes(result: str) -> list[str]:
    """Return trimmed slash-delimited tokens, or [] for an unparseable result.

    Require a slash and at least one token containing "*". Retain every nonempty
    token once that rule passes, without allele validation or copy expansion.
    """
    if "/" not in result:
        return []

    parts = [part.strip() for part in result.split("/") if part.strip()]
    haplotypes = parts

    # Filter out obvious non-calls (for example free-text warning lines)
    if not haplotypes or not any("*" in hap for hap in haplotypes):
        return []

    return haplotypes


def main() -> None:
    """Parse reports and write calls plus cohort-filtered frequencies and Table S1."""
    parser = argparse.ArgumentParser(
        description='Parse ".alleles" files and extract the value after "Result".'
    )
    parser.add_argument(
        "--input-dir",
        default="StellarPGxCalls",
        help='Directory containing ".alleles" files.',
    )
    parser.add_argument(
        "--pattern",
        default="*.alleles",
        help='Glob pattern used inside --input-dir (default: "*.alleles").',
    )
    parser.add_argument(
        "--output",
        default="alleles_results.csv",
        help="Output CSV path.",
    )
    parser.add_argument(
        "--haplotype-output",
        default="haplotype_frequencies.csv",
        help="Output CSV path for haplotype frequencies.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search recursively under --input-dir.",
    )
    parser.add_argument(
        "--sample-allowlist",
        default=DEFAULT_SAMPLE_ALLOWLIST,
        help=(
            "Text file of sample IDs to include in haplotype frequency calculations "
            f"(default: {DEFAULT_SAMPLE_ALLOWLIST}, relative to the working directory)."
        ),
    )
    parser.add_argument(
        "--sample-star-allele-output",
        default=DEFAULT_SAMPLE_STAR_ALLELE_OUTPUT,
        help=(
            "Shared Table S1 CSV; replace CYP2D6 rows and preserve other genes "
            f"(default: {DEFAULT_SAMPLE_STAR_ALLELE_OUTPUT})."
        ),
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    haplotype_output_path = Path(args.haplotype_output).expanduser().resolve()
    sample_star_allele_output_path = Path(args.sample_star_allele_output).expanduser().resolve()
    if len({output_path, haplotype_output_path, sample_star_allele_output_path}) != 3:
        parser.error("Calls, frequencies, and Table S1 must use distinct output paths.")
    sample_allowlist_path = (
        Path(args.sample_allowlist).expanduser().resolve()
        if args.sample_allowlist
        else None
    )

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    paths = (
        sorted(input_dir.rglob(args.pattern))
        if args.recursive
        else sorted(input_dir.glob(args.pattern))
    )

    if not paths:
        raise FileNotFoundError(
            f'No files matched pattern "{args.pattern}" in {input_dir}'
        )

    allowed_samples = None
    if sample_allowlist_path is not None:
        if not sample_allowlist_path.exists():
            raise FileNotFoundError(f"Sample allowlist not found: {sample_allowlist_path}")
        allowed_samples = {
            line.strip()
            for line in sample_allowlist_path.read_text().splitlines()
            if line.strip()
        }

    # Calls include every report; frequencies and Table S1 use the allowlist.
    rows = []
    sample_star_allele_rows = []
    missing_result_count = 0
    unresolved_result_count = 0
    filtered_out_count = 0
    haplotype_counter: Counter[str] = Counter()
    haplotype_sample_counter: Counter[str] = Counter()

    for path in paths:
        sample_id = sample_id_from_name(path.name)
        content = path.read_text(errors="replace")
        result = extract_labeled_value(content, "result")
        activity_score = extract_labeled_value(content, "activity_score")
        metabolizer_status = extract_labeled_value(content, "metabolizer_status")
        if result is None:
            missing_result_count += 1
            result = ""

        haplotypes = extract_haplotypes(result)
        include_for_frequency = allowed_samples is None or sample_id in allowed_samples
        if include_for_frequency:
            sample_star_allele_rows.append(
                {
                    "sample_id": sample_id,
                    "gene": "CYP2D6",
                    "source_diplotype_count": int(bool(haplotypes)),
                    "counted_alleles": len(haplotypes),
                    "labels_preview": result,
                }
            )
            if not haplotypes:
                unresolved_result_count += 1
            else:
                haplotype_counter.update(haplotypes)
                # Count presence once per report, including homozygous calls.
                for haplotype in sorted(set(haplotypes)):
                    haplotype_sample_counter[haplotype] += 1
        else:
            filtered_out_count += 1

        rows.append(
            {
                "sample_id": sample_id,
                "file_name": path.name,
                "file_path": str(path),
                "result": result,
                "activity_score": activity_score or "",
                "metabolizer_status": metabolizer_status or "",
                "haplotypes": "|".join(haplotypes),
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "file_name",
                "file_path",
                "result",
                "activity_score",
                "metabolizer_status",
                "haplotypes",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    # Normalize by parsed token occurrences, not the number of cohort samples.
    total_haplotype_calls = sum(haplotype_counter.values())
    haplotype_rows = []
    for haplotype, count in haplotype_counter.most_common():
        frequency = count / total_haplotype_calls if total_haplotype_calls else 0.0
        haplotype_rows.append(
            {
                "haplotype": haplotype,
                "count": count,
                "sample_count": haplotype_sample_counter[haplotype],
                "frequency": f"{frequency:.6f}",
            }
        )

    haplotype_output_path.parent.mkdir(parents=True, exist_ok=True)
    with haplotype_output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["haplotype", "count", "sample_count", "frequency"],
        )
        writer.writeheader()
        writer.writerows(haplotype_rows)

    write_sample_star_allele_table(
        sample_star_allele_output_path, sample_star_allele_rows, {"CYP2D6"}
    )

    print(f"Parsed files: {len(rows)}")
    print(f"Filtered out from frequency calculations: {filtered_out_count}")
    print(f"Missing result: {missing_result_count}")
    print(f"Unresolved result for haplotype parse: {unresolved_result_count}")
    print(f"Distinct haplotypes: {len(haplotype_counter)}")
    print(f"Total haplotype calls: {total_haplotype_calls}")
    print(f"Table S1 CYP2D6 rows: {len(sample_star_allele_rows)}")
    if sample_allowlist_path is not None:
        print(f"Using sample allowlist: {sample_allowlist_path}")
    print(f"Wrote: {output_path}")
    print(f"Wrote: {haplotype_output_path}")
    print(f"Updated: {sample_star_allele_output_path}")


if __name__ == "__main__":
    main()
