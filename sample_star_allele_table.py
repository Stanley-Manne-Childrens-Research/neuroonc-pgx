"""Write the shared Table S1 without replacing rows produced by other callers.

PharmCAT and StellarPGx use the same five-column schema. Each writer replaces
only its own genes, so sequential runs in either order preserve the other
caller's rows. Repeated runs replace prior rows rather than appending duplicates.
"""

import csv
from pathlib import Path


# Shared supplementary-table columns, in publication output order.
SAMPLE_STAR_ALLELE_FIELDS = [
    "sample_id",
    "gene",
    "source_diplotype_count",
    "counted_alleles",
    "labels_preview",
]


def write_sample_star_allele_table(
    output_path: Path,
    rows: list[dict[str, object]],
    genes_to_replace: set[str],
) -> None:
    """Replace selected genes while retaining all existing rows for other genes.

    Reject incompatible existing headers before opening the output for writing.
    Sort by gene, decreasing source-diplotype count, then sample ID, matching the
    original PharmCAT output. Duplicate input reports remain separate rows.
    """
    retained = []
    if output_path.exists():
        with output_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != SAMPLE_STAR_ALLELE_FIELDS:
                raise ValueError(
                    f"Unexpected Table S1 columns in {output_path}; "
                    f"expected {SAMPLE_STAR_ALLELE_FIELDS}"
                )
            retained = [row for row in reader if row["gene"] not in genes_to_replace]

    combined = sorted(
        retained + rows,
        key=lambda row: (row["gene"], -int(row["source_diplotype_count"]), row["sample_id"]),
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_STAR_ALLELE_FIELDS)
        writer.writeheader()
        writer.writerows(combined)
