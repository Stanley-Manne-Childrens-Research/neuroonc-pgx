"""Copy prepared CYP2D6 haplotype frequencies into an existing reference workbook.

Requires Python 3.10 or later and openpyxl. Run from the repository root:

    python3 update_cyp2d6_frequency_table.py

Upstream files and preparation
-----------------------------
The default input is Cyp2d6-Processing/haplotype_frequencies.csv. It must have a
header row with exact "haplotype" and "frequency" column names. Each data row
supplies an allele label and its observed frequency as a proportion, not a
percentage. Extra columns (including count and sample_count) are ignored. Labels
and frequency text are stripped of surrounding whitespace; blank labels are
skipped. Frequencies must be parseable by float(). Finite values in [0, 1] are
assumed upstream: this script does not check bounds, finiteness, or the total.

aggregate_stellarpgx_calls.py produces this CSV from existing StellarPGx .alleles
reports; parse_alleles_results.py is another compatible producer.
When aggregating from the repository root, explicitly select the expected path:

    python3 aggregate_stellarpgx_calls.py --haplotype-output Cyp2d6-Processing/haplotype_frequencies.csv

The aggregator otherwise writes haplotype_frequencies.csv in the working
directory. Alternatively, pass this updater --csv PATH to select that output.
Cohort selection belongs to the upstream aggregation step; the aggregator uses
sample-cohort-manifest.txt by default. This updater reads neither the manifest
nor individual calls and cannot verify cohort membership, counting assumptions,
or whether the CSV is current. It sums supplied frequencies without recomputing
them from counts or renormalizing after matching. Upstream rounding is retained.

The default target, DiplotypeFrequencies/CYP2D6_frequency_table.xlsx, must already
exist. Its "Allele frequency" sheet must have "CYP2D6 allele" in A2 and "CBTN"
in K2. Allele rows start on row 3; labels in column A are stripped and must be
unique. Rows with missing/blank labels are ignored. --workbook, --sheet, and
--column can select another target, but the exact header requirements remain.
The target column is selected by letter, not discovered from the CBTN header.
Reference workbooks and allele nomenclature must be prepared separately; this
script does not download tables, create allele rows, or track reference versions.
All relative paths resolve from the working directory, not the script directory.

Matching and aggregation assumptions
-----------------------------------
Exact allele-label matches take precedence. Otherwise, text after the first tab
is removed and the trimmed label is checked for a match. If still unmatched,
COPY_NUMBER_PATTERN may map a numeric copy-number suffix to an existing x≥3
workbook category (e.g., *1x3 -> *1x≥3). The current pattern requires the numeric
suffix to begin with 3-9: x3 and x30 qualify, but x10 and x20 do not. An exact
workbook match is still accepted for any of these labels. x2 is not grouped.
The grouping renames labels only; frequencies are not multiplied by copy number.

All CSV rows mapping to the same workbook label are summed, including duplicate
input rows. There is no input deduplication, so upstream rows must represent the
intended contributions. Unmatched labels are reported and omitted; they are not
added to the workbook or redistributed among matched alleles.

Outputs and failure behavior
----------------------------
Existing CBTN values for every recognized allele row are cleared before matched
sums are written. Unmatched workbook rows remain blank, not zero; explicit zero
input frequencies are written as zero. Cells on rows without allele labels are
not cleared. A valid header-only CSV therefore clears all recognized CBTN rows.
The workbook is saved in place without creating a backup; other worksheet values
are not intentionally changed. No additional CSV or plot is produced.

The console summary counts parsed nonblank-label CSV rows, exact/normalized
matches (input rows, not distinct workbook alleles), unmatched input rows, and
workbook allele rows left blank. Unmatched labels and blank rows are listed.
Unmatched labels and an empty frequency table do not prevent saving. Successful
completion does not imply that every observed allele was represented or that
the updated frequencies sum to one.

Missing input files, invalid CSV headers/numeric text, missing worksheets,
unexpected workbook headers, and duplicate workbook labels raise ValidationError;
main prints an error to stderr and returns exit status 1. These checks occur
before clearing values. Other exceptions, including unreadable workbooks,
invalid column letters, and save failures, propagate. The workbook must be
writable; this script does not create its parent directory.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.utils import column_index_from_string


# Prepared observed-frequency CSV with haplotype/frequency headers. Relative to
# the working directory; --csv can select another upstream aggregation output.
DEFAULT_CSV = Path("Cyp2d6-Processing/haplotype_frequencies.csv")

# Existing reference workbook to update in place; overridden with --workbook.
DEFAULT_WORKBOOK = Path("DiplotypeFrequencies/CYP2D6_frequency_table.xlsx")

# Exact name of the worksheet containing allele rows; overridden with --sheet.
DEFAULT_SHEET = "Allele frequency"

# Excel column letter for observed CBTN frequencies; overridden with --column.
# The chosen column must carry EXPECTED_CBTN_HEADER on Excel row 2.
DEFAULT_COLUMN = "K"

# Required allele-column header in A2, independent of the target CBTN column.
EXPECTED_ALLELE_HEADER = "CYP2D6 allele"

# Required row-2 header for Children's Brain Tumor Network (CBTN) frequencies.
EXPECTED_CBTN_HEADER = "CBTN"

# Capture an allele prefix and a trailing xN whose first digit is 3-9. Used to
# group unmatched labels into an existing <prefix>x≥3 row. This pattern does
# not recognize x10-x29 or other suffixes beginning with 1/2; exact matches still
# take precedence. Captured copy numbers are not used as frequency multipliers.
COPY_NUMBER_PATTERN = re.compile(r"^(.*)x([3-9][0-9]*)$")


class ValidationError(Exception):
    """Input/schema failure reported by main as an error with exit status 1."""

    pass


def parse_args() -> argparse.Namespace:
    """Read CSV, workbook, sheet, and target-column overrides from the CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Update the CBTN column in the CYP2D6 allele frequency workbook "
            "from haplotype_frequencies.csv."
        )
    )
    parser.add_argument(
        "--csv",
        dest="csv_path",
        default=str(DEFAULT_CSV),
        help="Input haplotype frequency CSV.",
    )
    parser.add_argument(
        "--workbook",
        default=str(DEFAULT_WORKBOOK),
        help="Target workbook to update in place.",
    )
    parser.add_argument(
        "--sheet",
        default=DEFAULT_SHEET,
        help="Workbook sheet to update.",
    )
    parser.add_argument(
        "--column",
        default=DEFAULT_COLUMN,
        help="Workbook column letter for the CBTN values.",
    )
    return parser.parse_args()


def ensure_file(path: Path, label: str) -> None:
    """Require an existing regular file, using label in validation errors."""
    if not path.exists():
        raise ValidationError(f"{label} not found: {path}")
    if not path.is_file():
        raise ValidationError(f"{label} is not a file: {path}")


def load_csv_rows(csv_path: Path) -> list[dict[str, object]]:
    """Read stripped nonblank haplotype labels and float-convertible frequencies.

    Preserve duplicate labels and ignore extra columns. Numerical bounds,
    finiteness, and frequency totals are not validated here.
    """
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValidationError(f"CSV has no header row: {csv_path}")

        missing = {"haplotype", "frequency"} - set(reader.fieldnames)
        if missing:
            missing_text = ", ".join(sorted(missing))
            raise ValidationError(
                f"CSV missing required column(s): {missing_text}"
            )

        rows = []
        for row in reader:
            haplotype = (row.get("haplotype") or "").strip()
            frequency_text = (row.get("frequency") or "").strip()
            if not haplotype:
                continue
            try:
                frequency = float(frequency_text)
            except ValueError as exc:
                raise ValidationError(
                    f"Invalid frequency for haplotype {haplotype}: {frequency_text}"
                ) from exc
            rows.append({"haplotype": haplotype, "frequency": frequency})

    return rows


def load_workbook_state(
    workbook_path: Path, sheet_name: str, column_letter: str
) -> tuple[object, object, dict[str, int], int]:
    """Validate the target and return workbook, worksheet, row lookup, and column.

    Lookup keys are stripped column-A labels; row and column indices are
    one-based. Reject duplicate nonblank labels before changing any CBTN cells.
    """
    workbook = load_workbook(workbook_path)
    if sheet_name not in workbook.sheetnames:
        raise ValidationError(
            f"Workbook sheet not found: {sheet_name!r} in {workbook_path}"
        )

    worksheet = workbook[sheet_name]
    column_index = column_index_from_string(column_letter)

    if worksheet["A2"].value != EXPECTED_ALLELE_HEADER:
        raise ValidationError(
            f"Unexpected A2 header: {worksheet['A2'].value!r}; "
            f"expected {EXPECTED_ALLELE_HEADER!r}"
        )

    header_cell = worksheet.cell(row=2, column=column_index)
    if header_cell.value != EXPECTED_CBTN_HEADER:
        raise ValidationError(
            f"Unexpected {column_letter}2 header: {header_cell.value!r}; "
            f"expected {EXPECTED_CBTN_HEADER!r}"
        )

    row_lookup: dict[str, int] = {}
    for row_index in range(3, worksheet.max_row + 1):
        allele = worksheet.cell(row=row_index, column=1).value
        if allele is None:
            continue
        allele_label = str(allele).strip()
        if not allele_label:
            continue
        if allele_label in row_lookup:
            raise ValidationError(f"Duplicate allele row in workbook: {allele_label}")
        row_lookup[allele_label] = row_index

    return workbook, worksheet, row_lookup, column_index


def normalize_haplotype_label(haplotype: str, row_lookup: dict[str, int]) -> str | None:
    """Try tab-suffix removal, then copy-number grouping into existing rows.

    Exact original-label matching is performed by apply_updates first. Return
    None when neither normalization produces a label present in row_lookup.
    """
    label = haplotype.strip()

    if "\t" in label:
        trimmed = label.split("\t", 1)[0].strip()
        if trimmed in row_lookup:
            return trimmed
        label = trimmed

    match = COPY_NUMBER_PATTERN.match(label)
    if match:
        normalized = f"{match.group(1)}x≥3"
        if normalized in row_lookup:
            return normalized

    return None


def clear_existing_values(
    worksheet: object, row_lookup: dict[str, int], column_index: int
) -> None:
    """Blank target-column cells for labeled allele rows in memory."""
    for row_index in row_lookup.values():
        worksheet.cell(row=row_index, column=column_index).value = None


def apply_updates(
    worksheet: object,
    row_lookup: dict[str, int],
    column_index: int,
    csv_rows: list[dict[str, object]],
) -> tuple[int, int, list[tuple[str, float]], dict[str, float]]:
    """Sum exact/normalized input frequencies by target allele and write in memory.

    Return exact input-row count, normalized input-row count, unmatched
    (label, frequency) pairs, and the per-target sums. Saving is the caller's job.
    """
    exact_matches = 0
    normalized_matches = 0
    unmatched: list[tuple[str, float]] = []
    aggregated_values: dict[str, float] = {}

    for row in csv_rows:
        haplotype = str(row["haplotype"])
        frequency = float(row["frequency"])

        if haplotype in row_lookup:
            aggregated_values[haplotype] = aggregated_values.get(haplotype, 0.0) + frequency
            exact_matches += 1
            continue

        normalized = normalize_haplotype_label(haplotype, row_lookup)
        if normalized is not None:
            aggregated_values[normalized] = aggregated_values.get(normalized, 0.0) + frequency
            normalized_matches += 1
            continue

        unmatched.append((haplotype, frequency))

    for allele, frequency in aggregated_values.items():
        row_index = row_lookup[allele]
        worksheet.cell(row=row_index, column=column_index).value = frequency

    return exact_matches, normalized_matches, unmatched, aggregated_values


def get_blank_rows(
    worksheet: object, row_lookup: dict[str, int], column_index: int
) -> list[str]:
    """List target allele labels whose cells are None or empty; retain zeros."""
    blank_labels = []
    for allele, row_index in row_lookup.items():
        if worksheet.cell(row=row_index, column=column_index).value in (None, ""):
            blank_labels.append(allele)
    return blank_labels


def print_summary(
    total_rows: int,
    exact_matches: int,
    normalized_matches: int,
    unmatched: list[tuple[str, float]],
    blank_rows: list[str],
) -> None:
    """Print input-row match counts, unmatched entries, and blank workbook labels."""
    print("Update summary:")
    print(f"  total_csv_rows={total_rows}")
    print(f"  exact_matches={exact_matches}")
    print(f"  normalized_matches={normalized_matches}")
    print(f"  unmatched_csv_haplotypes={len(unmatched)}")
    print(f"  workbook_rows_left_blank={len(blank_rows)}")

    if unmatched:
        print("\nUnmatched CSV haplotypes:")
        for haplotype, frequency in unmatched:
            print(f"  {haplotype}: {frequency:.6f}")

    if blank_rows:
        print("\nWorkbook allele rows left blank:")
        for allele in blank_rows:
            print(f"  {allele}")


def main() -> int:
    """Validate inputs, replace target frequencies, save in place, and report status."""
    args = parse_args()
    csv_path = Path(args.csv_path).expanduser().resolve()
    workbook_path = Path(args.workbook).expanduser().resolve()

    try:
        ensure_file(csv_path, "CSV")
        ensure_file(workbook_path, "Workbook")
        csv_rows = load_csv_rows(csv_path)
        workbook, worksheet, row_lookup, column_index = load_workbook_state(
            workbook_path, args.sheet, args.column
        )

        clear_existing_values(worksheet, row_lookup, column_index)
        exact_matches, normalized_matches, unmatched, _ = apply_updates(
            worksheet, row_lookup, column_index, csv_rows
        )
        blank_rows = get_blank_rows(worksheet, row_lookup, column_index)
        workbook.save(workbook_path)
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Updated workbook: {workbook_path}")
    print_summary(
        total_rows=len(csv_rows),
        exact_matches=exact_matches,
        normalized_matches=normalized_matches,
        unmatched=unmatched,
        blank_rows=blank_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
