"""Count CBTN allele occurrences from PharmCAT reports and update frequency tables.

Run from the repository root: python3 get_cbtn_allele_frequencies.py
Requires Python, openpyxl, and the local sample_star_allele_table.py helper.
All configured paths are relative to the current
working directory; there are no command-line options or automatic path discovery.

Upstream files and preparation
-----------------------------
* PharmCatCalls/*.subset.phenotype.json: existing PharmCAT phenotype
  reports, discovered nonrecursively. Each selected JSON must contain a mapping
  at geneReports -> CPIC, keyed by gene name. A gene report's sourceDiplotypes
  list supplies allele1/allele2 objects with a name field; label is optional and
  used only in debug output. Missing sourceDiplotypes defaults to an empty list.
  The script does not generate PharmCAT reports or read raw VCFs.
* shared_pharmcat_cyp2d6_samples.txt: required, newline-delimited canonical sample
  IDs. Blank lines are ignored and repeated manifest entries collapse to one ID.
  build_shared_pharmcat_cyp2d6_samples.py creates this file from metadata.csv,
  Cyp2d6-Processing/alleles_results.csv, and the PharmCAT JSON filenames. Those
  metadata/CYP2D6 files are needed to build the manifest, not read by this script.
* DiplotypeFrequencies/<GENE>_frequency_table.xlsx: existing reference workbooks
  to update in place. Each needs the "Allele frequency" sheet, an exact "CBTN"
  header on Excel row 2, and allele labels in column A starting on row 3. The
  reference tables and their allele nomenclature must be prepared separately;
  this script does not download tables or create workbooks, sheets, or rows.
  Workbook source/version information is not tracked by this script.

Counting assumptions
--------------------
The report filename, after removing its suffix and an optional _<batch>_ prefix,
identifies the sample. Only IDs in the manifest are included. Input reports must
already be deduplicated: multiple files for the same canonical ID are each
counted, including filename aliases. Debug "samples" counts gene reports rather
than distinct canonical IDs. Manifest membership alone does not validate a call.

Every named allele1/allele2 entry in every source diplotype contributes one
occurrence. Alternative diplotypes are all counted, without selecting one or
weighting ambiguity. There is no additional call-quality or unknown-label filter,
and copy-number labels are not expanded into extra allele copies. Each gene's
denominator is its total counted occurrences, not twice the manifest size;
denominators can differ between genes. Missing genes or allele entries contribute
no occurrences, rather than explicit zero-frequency rows.

For DPYD, bracketed combinations separated by " + " are split into component
labels. Each component receives the combination's count, but the denominator
remains the original count before splitting. Component frequencies can therefore
sum to more than one; they do not describe mutually exclusive haplotypes.

Outputs and update behavior
---------------------------
alleles_freqs.csv is overwritten without a header. Its five columns are gene,
allele, frequency (a proportion), occurrence count, and gene-specific denominator.
The gene summary CSV and table-s1-sample_star_allele.csv include headers and
report counting/ambiguity diagnostics. Table S1 replaces rows for GENES while
preserving CYP2D6 rows produced by aggregate_stellarpgx_calls.py. Repeated runs
replace prior rows for the respective genes; duplicate input reports are still
retained. The existing Table S1 header must match its five-column schema.
Output parent directories must already exist and destinations must be writable.

For each encountered target gene, the workbook's CBTN cells from row 3 onward
are cleared, then filled by exact allele-label matches after the aliases below.
Unmatched workbook labels remain blank, not zero; observed labels without a
workbook row remain in the CSV but are not added to the workbook. Genes absent
from all selected reports are not updated, so earlier workbook values may remain.
CYP2D6 is intentionally excluded: its workbook is populated separately by
Cyp2d6-Processing/update_cyp2d6_frequency_table.py from haplotype_frequencies.csv.

A missing manifest or malformed selected JSON raises an error. No matching JSON
files or no selected samples produces empty frequency output without an error.
Missing/unreadable workbooks, missing sheets, or missing CBTN headers produce
warnings and are skipped; CSV output is written before workbook updates. A run
does not guarantee that every manifest sample or target gene was represented.
"""

import csv
import glob
import json
import os
from pathlib import Path

from openpyxl import load_workbook

from sample_star_allele_table import write_sample_star_allele_table


# PharmCAT gene reports to count and export. CYP2D6 uses a separate calling and
# frequency-update workflow; reports for other genes are ignored.
GENES = [
    "ABCG2",
    "CYP2B6",
    "CYP2C19",
    "CYP2C9",
    "CYP3A5",
    "CYP4F2",
    "DPYD",
    "G6PD",
    "NUDT15",
    "SLCO1B1",
    "TPMT",
    "VKORC1",
    "UGT1A1",
    "CACNA1S",
    "RYR1"
]

# Input report pattern. The basename before .subset.phenotype.json supplies the
# sample ID; _<numeric batch>_ prefixes are removed before manifest lookup.
PHENOTYPE_JSON_GLOB = "PharmCatCalls/*.subset.phenotype.json"

# Main output: headerless, fully quoted CSV with columns
# gene, allele, frequency, occurrence_count, denominator. Overwritten each run.
ALLELE_FREQS_CSV = Path("alleles_freqs.csv")

# Per-gene counts and averages of source diplotypes and named allele entries,
# including the largest diplotype count and its sample ID. Includes a header.
GENE_SUMMARY_CSV = Path("alleles_freqs_gene_summary.csv")

# Per-report/per-gene diagnostics: canonical sample ID, source-diplotype count,
# counted allele entries, and up to ten source labels. Includes a header.
# Shared with StellarPGx aggregation: preserve rows for genes outside GENES.
SAMPLE_STAR_ALLELE_CSV = Path("table-s1-sample_star_allele.csv")

# Existing reference workbooks, named <GENE>_frequency_table.xlsx. Matched CBTN
# values are updated in place; ancestry reference frequencies are preserved. These files are obtained from ClinPGx.org.
DIPLOTYPE_FREQUENCIES_DIR = Path("DiplotypeFrequencies")

# Manifest of speciemns used for this cohort. Required sample allowlist shared
# with the CYP2D6 analysis; one canonical ID per line. This script loads the
# manifest but does not create or refresh it.
SHARED_SAMPLE_MANIFEST = Path("sample-cohort-manifest.txt")

# Exact worksheet name containing allele labels and ancestry/CBTN frequencies.
ALLELE_FREQUENCY_SHEET = "Allele frequency"

# Exact Excel row-2 header identifying the observed Children's Brain Tumor
# Network (CBTN) frequency column; its position is located independently per file. This column is appended to the reference workbooks and is not used for any calculations.
CBTN_HEADER = "CBTN"

# Gene-specific translations from counted PharmCAT labels (keys) to reference
# workbook labels (values). Applied only when filling workbooks; CSV labels stay
# unchanged. Unlisted labels require exact matches. This is a rename operation,
# not aggregation: if several labels map to one target, the last value wins.
WORKBOOK_ALLELE_ALIASES = {
    "ABCG2": {
        "rs2231142 reference (G)": "rs2231142 (G)",
        "rs2231142 variant (T)": "rs2231142 (T)",
    },
    "DPYD": {
        "c.85T>C (*9A)": "*9A",
        "c.1601G>A (*4)": "*4",
        "c.1627A>G (*5)": "*5",
        "c.1905+1G>A (*2A)": "*2A",
        "c.2194G>A (*6)": "*6",
        "c.1129-5923C>G, c.1236G>A (HapB3)": "HapB3",
    },
}


def canonicalize_sample_id(sample_id):
    """Remove a leading _<numeric batch>_ prefix; otherwise preserve the ID."""
    if sample_id.startswith("_"):
        parts = sample_id.split("_", 2)
        if len(parts) == 3 and parts[1].isdigit():
            return parts[2]
    return sample_id


def read_shared_sample_ids():
    """Read unique, nonblank canonical IDs from the required sample manifest."""
    return {
        line.strip()
        for line in SHARED_SAMPLE_MANIFEST.read_text().splitlines()
        if line.strip()
    }


def collect_allele_counts(shared_sample_ids):
    """Return raw per-gene occurrence counts and gene/report debug records.

    All source diplotypes in each selected report contribute counts. Canonical
    IDs control inclusion but are not used to deduplicate input report files.
    """
    by_gene = {}
    gene_debug = {}
    sample_debug_rows = []
    paths = sorted(glob.glob(PHENOTYPE_JSON_GLOB))

    for path in paths:
        raw_patient = os.path.basename(path).split(".subset.phenotype.json")[0]
        patient = canonicalize_sample_id(raw_patient)
        if patient not in shared_sample_ids:
            print(f"SKIP_NON_SHARED\tsample={raw_patient}\tcanonical={patient}")
            continue

        with open(path, "r") as file:
            data = json.load(file)

        for gene in data["geneReports"]["CPIC"]:
            if gene not in GENES:
                continue

            source_diplotypes = data["geneReports"]["CPIC"][gene].get("sourceDiplotypes", [])
            gene_counts = by_gene.setdefault(gene, {})
            counted_alleles = 0
            labels = []
            for diplo in source_diplotypes:
                label = diplo.get("label")
                if label:
                    labels.append(label)
                for allele_key in ("allele1", "allele2"):
                    allele_info = diplo.get(allele_key)
                    if allele_info and "name" in allele_info:
                        allele_name = allele_info["name"]
                        gene_counts[allele_name] = gene_counts.get(allele_name, 0) + 1
                        counted_alleles += 1

            stats = gene_debug.setdefault(
                gene,
                {
                    "samples": 0,
                    "total_source_diplotypes": 0,
                    "total_counted_alleles": 0,
                    "max_source_diplotypes": 0,
                    "max_source_diplotype_sample": "",
                },
            )
            stats["samples"] += 1
            stats["total_source_diplotypes"] += len(source_diplotypes)
            stats["total_counted_alleles"] += counted_alleles
            if len(source_diplotypes) > stats["max_source_diplotypes"]:
                stats["max_source_diplotypes"] = len(source_diplotypes)
                stats["max_source_diplotype_sample"] = patient

            sample_debug_rows.append(
                {
                    "sample_id": patient,
                    "gene": gene,
                    "source_diplotype_count": len(source_diplotypes),
                    "counted_alleles": counted_alleles,
                    "labels_preview": " | ".join(labels[:10]),
                }
            )

    return by_gene, gene_debug, sample_debug_rows


def collapse_dpyd_counts(gene_counts):
    """Split bracketed DPYD combinations while retaining the raw denominator.

    A combination's count is assigned to every component; "Reference" and
    nonbracketed labels are retained as individual labels.
    """
    collapsed_counts = {}
    raw_total = sum(gene_counts.values())

    for allele_name, count in gene_counts.items():
        if allele_name == "Reference":
            labels = ["Reference"]
        elif allele_name.startswith("[") and allele_name.endswith("]"):
            labels = [label.strip() for label in allele_name[1:-1].split(" + ") if label.strip()]
        else:
            labels = [allele_name]

        for label in labels:
            collapsed_counts[label] = collapsed_counts.get(label, 0) + count

    return collapsed_counts, raw_total


def build_frequency_rows(by_gene):
    """Return sorted five-field CSV rows and per-gene allele-frequency mappings.

    Frequencies are occurrence count / gene-specific denominator, with DPYD
    component counts using the denominator before combination splitting.
    """
    rows = []
    frequencies_by_gene = {}

    for gene in sorted(by_gene):
        gene_counts = by_gene[gene]
        total = sum(gene_counts.values())
        if gene == "DPYD":
            gene_counts, total = collapse_dpyd_counts(gene_counts)

        frequencies_by_gene[gene] = {}
        for allele in sorted(gene_counts):
            frequency = float(gene_counts[allele]) / float(total)
            rows.append([gene, allele, frequency, gene_counts[allele], total])
            frequencies_by_gene[gene][allele] = frequency

    return rows, frequencies_by_gene


def write_csv(rows):
    """Overwrite the main frequency CSV, quoting all fields and omitting headers."""
    with open(ALLELE_FREQS_CSV, "w", newline="") as outfile:
        writer = csv.writer(outfile, dialect="excel", quoting=csv.QUOTE_ALL)
        writer.writerows(rows)
    print(f"WROTE\tcsv\t{ALLELE_FREQS_CSV}")


def write_debug_outputs(gene_debug, sample_debug_rows):
    """Write the gene summary, replace PharmCAT Table S1 rows, and print counts."""
    with open(GENE_SUMMARY_CSV, "w", newline="") as outfile:
        writer = csv.DictWriter(
            outfile,
            fieldnames=[
                "gene",
                "samples",
                "total_source_diplotypes",
                "total_counted_alleles",
                "avg_source_diplotypes_per_sample",
                "avg_counted_alleles_per_sample",
                "max_source_diplotypes",
                "max_source_diplotype_sample",
            ],
        )
        writer.writeheader()
        for gene in sorted(gene_debug):
            stats = gene_debug[gene]
            samples = stats["samples"]
            writer.writerow(
                {
                    "gene": gene,
                    "samples": samples,
                    "total_source_diplotypes": stats["total_source_diplotypes"],
                    "total_counted_alleles": stats["total_counted_alleles"],
                    "avg_source_diplotypes_per_sample": (
                        stats["total_source_diplotypes"] / samples if samples else 0.0
                    ),
                    "avg_counted_alleles_per_sample": (
                        stats["total_counted_alleles"] / samples if samples else 0.0
                    ),
                    "max_source_diplotypes": stats["max_source_diplotypes"],
                    "max_source_diplotype_sample": stats["max_source_diplotype_sample"],
                }
            )

    write_sample_star_allele_table(SAMPLE_STAR_ALLELE_CSV, sample_debug_rows, set(GENES))

    print(f"WROTE\tdebug\t{GENE_SUMMARY_CSV}")
    print(f"WROTE\tdebug\t{SAMPLE_STAR_ALLELE_CSV}")

    for gene in sorted(gene_debug):
        stats = gene_debug[gene]
        samples = stats["samples"]
        avg_diplotypes = stats["total_source_diplotypes"] / samples if samples else 0.0
        avg_alleles = stats["total_counted_alleles"] / samples if samples else 0.0
        print(
            "DEBUG_GENE\t"
            f"gene={gene}\t"
            f"samples={samples}\t"
            f"total_source_diplotypes={stats['total_source_diplotypes']}\t"
            f"total_counted_alleles={stats['total_counted_alleles']}\t"
            f"avg_source_diplotypes={avg_diplotypes:.3f}\t"
            f"avg_counted_alleles={avg_alleles:.3f}\t"
            f"max_source_diplotypes={stats['max_source_diplotypes']}\t"
            f"max_sample={stats['max_source_diplotype_sample']}"
        )


def get_cbtn_column_index(sheet):
    """Return the first exact CBTN header's one-based column, or None."""
    for column_index in range(1, sheet.max_column + 1):
        if sheet.cell(row=2, column=column_index).value == CBTN_HEADER:
            return column_index
    return None


def build_workbook_frequency_map(gene, frequencies):
    """Translate allele names for workbook lookup without summing alias collisions."""
    workbook_frequencies = {}
    aliases = WORKBOOK_ALLELE_ALIASES.get(gene, {})

    for allele_name, frequency in frequencies.items():
        workbook_allele_name = aliases.get(allele_name, allele_name)
        workbook_frequencies[workbook_allele_name] = frequency

    return workbook_frequencies


def update_workbook_for_gene(gene, frequencies):
    """Clear and refill an existing workbook's CBTN column, then save in place.

    Missing/unreadable workbooks or required sheet/header omissions are reported
    and skipped. Allele rows are matched exactly; unmatched cells remain blank.
    """
    workbook_path = DIPLOTYPE_FREQUENCIES_DIR / f"{gene}_frequency_table.xlsx"
    if not workbook_path.exists():
        print(f"WARNING\tmissing_workbook\t{gene}\t{workbook_path}")
        return

    try:
        workbook = load_workbook(workbook_path)
    except Exception as exc:
        print(f"WARNING\tinvalid_workbook\t{gene}\t{workbook_path}\t{exc}")
        return

    if ALLELE_FREQUENCY_SHEET not in workbook.sheetnames:
        print(f"WARNING\tmissing_sheet\t{gene}\t{workbook_path}\t{ALLELE_FREQUENCY_SHEET}")
        return

    sheet = workbook[ALLELE_FREQUENCY_SHEET]
    cbtn_column = get_cbtn_column_index(sheet)
    if cbtn_column is None:
        print(f"WARNING\tmissing_cbtn_header\t{gene}\t{workbook_path}")
        return

    workbook_frequencies = build_workbook_frequency_map(gene, frequencies)
    matched = 0
    for row_index in range(3, sheet.max_row + 1):
        allele_name = sheet.cell(row=row_index, column=1).value
        sheet.cell(row=row_index, column=cbtn_column).value = None
        if allele_name in workbook_frequencies:
            sheet.cell(row=row_index, column=cbtn_column).value = workbook_frequencies[allele_name]
            matched += 1

    workbook.save(workbook_path)
    print(
        f"UPDATED\tworkbook\t{gene}\tmatched={matched}\t"
        f"alleles={len(workbook_frequencies)}\tpath={workbook_path}"
    )


def update_workbooks(frequencies_by_gene):
    """Update only genes represented in the collected frequency mappings."""
    for gene, frequencies in sorted(frequencies_by_gene.items()):
        update_workbook_for_gene(gene, frequencies)


def main():
    """Load the cohort, count occurrences, write CSVs, then update workbooks."""
    shared_sample_ids = read_shared_sample_ids()
    print(f"USING_SHARED_SAMPLES\tcount={len(shared_sample_ids)}\tpath={SHARED_SAMPLE_MANIFEST}")
    by_gene, gene_debug, sample_debug_rows = collect_allele_counts(shared_sample_ids)
    rows, frequencies_by_gene = build_frequency_rows(by_gene)
    write_csv(rows)
    write_debug_outputs(gene_debug, sample_debug_rows)
    update_workbooks(frequencies_by_gene)


if __name__ == "__main__":
    main()
