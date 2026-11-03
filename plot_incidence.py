"""Plot cohort phenotype distributions and fractions classified as actionable.

Run from the repository root: python3 PublicationResources/plot_incidence.py
Requires Python 3.10 or later, pandas, numpy, matplotlib, and the adjacent local
phenotype_plot_utils.py helper. There are no command-line options. Despite the
historical filename, these are sample-level proportions, not longitudinal
incidence rates or proportions of distinct patients.

Upstream files and path assumptions
----------------------------------
All input and output paths are relative to this script's resolved directory,
currently PublicationResources/, rather than the working directory. Required
data files must be placed or linked at those locations before running:

* alleles_and_phenotypes.csv: headerless PharmCAT-derived table with at least nine
  columns. Zero-based columns 0, 5, 6, 7, and 8 contain Sample ID, Authority, Gene,
  Diplotype, and Phenotype respectively; intervening demographic fields are
  ignored. The repository's parse_pharmcat_outputs.py produces this format from
  PharmCAT JSON reports and metadata.csv, but writes to its working directory.
  This plotting script does not run that preparation step or parse raw calls.
* shared_pharmcat_cyp2d6_samples.txt: existing cohort manifest with one exact
  sample ID per line. Whitespace and blank lines are removed and duplicates
  collapse to one ID. The legacy manifest filename is still configured here;
  sample-cohort-manifest.txt is not selected automatically. No sample-alias
  canonicalization or patient-level merging occurs in this plotting script.
* Cyp2d6-Processing/alleles_results.csv: per-sample StellarPGx aggregation output
  with sample_id, result, and metabolizer_status headers. The helper reads this
  file when the selected phenotype table contains CYP2D6 rows. With the current
  script location, the expected subdirectory is inside PublicationResources/.

The adjacent phenotype_plot_utils.py defines the gene list, phenotype
normalization, actionability sets, and conflict-resolution priorities. Keep this
helper with the script and retain its version with the analysis. Individual
PharmCAT/StellarPGx report files, allele-frequency workbooks, treatment records,
and Table S1 are not read by this script.

Phenotype selection and counting assumptions
-------------------------------------------
The helper selects DPWG annotations for VKORC1 and CPIC annotations for other
genes; VKORC1 has no CPIC fallback. CYP4F2 uses study-defined categories derived
from CPIC diplotypes (*1/*3, *3/*3, otherwise Normal Metabolizer for valid parsed
calls). Existing CYP2D6 rows are supplemented/replaced with matched StellarPGx
results and mapped metabolizer statuses; the first calls row per sample is used.
Unmatched CYP2D6 rows retain their original annotations. Calls for samples with
no existing CYP2D6 phenotype-table row are not added by this replacement step.

Only manifest-listed samples are retained. The helper resolves multiple rows
to one phenotype per sample/gene: exact duplicates collapse, configured
gene-specific priorities resolve listed conflicts, and unresolved conflicts
retain the first input phenotype with a warning. Missing phenotype text becomes
"No Result"; missing sample/gene combinations are not synthesized.

Actionability follows the helper's explicit study categories, not a live lookup
of prescribing guidelines or treatment exposure. Some genes include Indeterminate
as actionable and others do not. Unrecognized labels produce warnings and are
treated as non-actionable. Review the helper's sets to reproduce this definition.

Per-gene fraction = actionable resolved samples / samples with a resolved row
for that gene, including No Result rows. Denominators may differ by gene and
need not equal the manifest size. Aggregate fraction = samples with at least one
actionable gene / samples represented anywhere in the resolved data. Samples
absent from every gene do not enter the aggregate denominator. In the Boolean
sample table, a missing gene row is False and cannot be distinguished from a
resolved non-actionable phenotype. A nonempty resolved cohort is assumed.

Outputs and display conventions
-------------------------------
Existing outputs in the script directory are overwritten:

* actionable_phenotypes_table.csv: one row per represented Sample ID, Boolean
  columns for every GENES entry, and Aggregate (any actionable gene in GENES).
* phenotype_data_table.csv: per-gene counts, denominators, and fractions, plus
  Aggregate across GENES. Historical "Interesting Phenotypes" column names mean
  actionable phenotypes under the helper's study definition.
* actionable_phenotypes_fraction_with_aggregate_and_stacked.png: per-gene bars
  for INTERACTION_LIST plus their Aggregate; grey segments show the complement.
  INTERACTION_LIST currently includes all GENES, so this matches the CSV
  aggregate. Narrowing it changes only the plotted subset/aggregate, not the
  previously saved sample table or full-gene summary.
* all_phenotype_pie_charts.png and all_phenotype_pie_charts.svg: a fixed 4-by-4
  gene grid. Categories below 0.05% of the original gene total are omitted;
  displayed percentages are renormalized over retained categories and rounded
  to one decimal place. Labels below a displayed 8% move outside with leader
  lines. Actionable categories use tab10 colors; others use greys. Colors are
  assigned independently per gene, not shared by phenotype across all panels.

PNG exports use EXPORT_DPI; SVG is also saved for the pie grid. The fixed grid
supports at most 16 genes. Genes with no resolved rows have hidden pie panels
and zero summary fractions. The output directory must be writable. Missing
files, invalid input schemas, and other loading/plotting errors propagate;
the script does not create inputs or refresh the cohort manifest.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

# Local study definitions and loaders: GENES controls table/pie order, while
# helper functions define actionability and one-phenotype-per-sample resolution.
from phenotype_plot_utils import (
    GENES,
    is_actionable_phenotype,
    load_cpic_phenotypes,
    resolve_entity_gene_phenotypes,
)

# Genes included in the stacked bar chart and its Aggregate calculation. This
# is an independent copy of GENES; CSVs and pies continue to use the full GENES.
INTERACTION_LIST = list(GENES)

# Resolution in dots per inch for both PNG figures; does not set SVG dimensions.
EXPORT_DPI = 300

# Required legacy cohort-manifest filename, resolved alongside this script.
# This is not the newer sample-cohort-manifest.txt used by allele-count scripts.
SHARED_SAMPLE_MANIFEST = "shared_pharmcat_cyp2d6_samples.txt"

# Non-actionable pie-category colors, cycled in label order within each gene.
# Actionable categories instead cycle through Matplotlib's tab10 palette.
NON_ACTIONABLE_GREY_PALETTE = [
    "#d9d9d9",
    "#c7c7c7",
    "#b5b5b5",
    "#a3a3a3",
    "#919191",
]


def phenotype_sort_key(gene: str, phenotype: str) -> tuple[int, str]:
    """Order actionable labels first, then alphabetically within each category."""
    actionable_rank = 0 if is_actionable_phenotype(gene, phenotype) else 1
    return (actionable_rank, phenotype)


def build_color_map(gene: str, labels: list[str]) -> dict[str, object]:
    """Assign per-gene colors to actionable labels and greys to other labels.

    Each palette cycles independently if the number of labels exceeds its size.
    """
    color_map: dict[str, object] = {}
    color_idx = 0
    grey_idx = 0
    for label in labels:
        if is_actionable_phenotype(gene, label):
            color_map[label] = plt.cm.tab10.colors[color_idx % len(plt.cm.tab10.colors)]
            color_idx += 1
        else:
            color_map[label] = NON_ACTIONABLE_GREY_PALETTE[
                grey_idx % len(NON_ACTIONABLE_GREY_PALETTE)
            ]
            grey_idx += 1
    return color_map


def annotate_small_slices(ax, wedges, autotexts) -> None:
    """Move displayed percentages below 8% outside the pie with leader lines.

    Use the rounded label text to select slices, group labels by pie side, and
    apply a minimum vertical gap plus bounded shifts to reduce label overlap.
    """
    small_slices = []
    for idx, autotext in enumerate(autotexts):
        try:
            pct = float(autotext.get_text().replace("%", ""))
        except ValueError:
            pct = 0
        if pct >= 8:
            continue

        angle = (wedges[idx].theta2 + wedges[idx].theta1) / 2
        line_r = 1.03
        angle_rad = np.deg2rad(angle)
        line_x = line_r * np.cos(angle_rad)
        line_y = line_r * np.sin(angle_rad)
        side = 1 if np.cos(angle_rad) >= 0 else -1
        small_slices.append(
            {
                "side": side,
                "line_x": line_x,
                "line_y": line_y,
                "target_y": 1.34 * np.sin(angle_rad),
                "text": autotext.get_text(),
            }
        )
        autotext.set_text("")

    # Layout distances use pie/axes data coordinates; all percentages are labels,
    # not proportions used to recompute the underlying phenotype counts.
    label_x = 1.46
    elbow_x = 1.16
    min_gap = 0.18
    y_min = -1.26
    y_max = 1.26

    for side in (-1, 1):
        side_slices = [slice_info for slice_info in small_slices if slice_info["side"] == side]
        side_slices.sort(key=lambda slice_info: slice_info["target_y"], reverse=True)
        previous_y = None
        for slice_info in side_slices:
            label_y = slice_info["target_y"]
            if previous_y is not None and previous_y - label_y < min_gap:
                label_y = previous_y - min_gap
            slice_info["label_y"] = label_y
            previous_y = label_y

        if side_slices and side_slices[-1]["label_y"] < y_min:
            shift = y_min - side_slices[-1]["label_y"]
            for slice_info in side_slices:
                slice_info["label_y"] += shift
        if side_slices and side_slices[0]["label_y"] > y_max:
            shift = side_slices[0]["label_y"] - y_max
            for slice_info in side_slices:
                slice_info["label_y"] -= shift

        for slice_info in side_slices:
            x = side * label_x
            mid_x = side * elbow_x
            label_y = slice_info["label_y"]
            ax.plot(
                [slice_info["line_x"], mid_x, x],
                [slice_info["line_y"], label_y, label_y],
                color="grey",
                lw=1.5,
            )
            ax.text(
                x,
                label_y,
                slice_info["text"],
                ha="left" if side > 0 else "right",
                va="center",
                fontsize=15,
            )


def format_pct(pct: float) -> str:
    """Format a pie percentage to one decimal; hide labels below 0.05 percent."""
    return f"{pct:.1f}%" if pct >= 0.05 else ""


def read_shared_sample_ids(manifest_path: Path) -> set[str]:
    """Read unique, nonblank sample IDs without further alias normalization."""
    return {
        line.strip()
        for line in manifest_path.read_text().splitlines()
        if line.strip()
    }


def main() -> None:
    """Load cohort phenotypes, resolve/classify calls, and export tables and plots."""
    # Paths follow this script's location, even when launched from another folder.
    base_dir = Path(__file__).resolve().parent
    workbook_path = base_dir / "alleles_and_phenotypes.csv"
    cyp2d6_calls_path = base_dir / "Cyp2d6-Processing" / "alleles_results.csv"
    shared_sample_manifest_path = base_dir / SHARED_SAMPLE_MANIFEST
    shared_sample_ids = read_shared_sample_ids(shared_sample_manifest_path)

    # Historical helper name: VKORC1 actually uses DPWG; other genes use CPIC.
    phenotypes = load_cpic_phenotypes(
        workbook_path,
        genes=GENES,
        cyp2d6_calls_path=cyp2d6_calls_path,
    )
    phenotypes = phenotypes[phenotypes["Sample ID"].isin(shared_sample_ids)].copy()
    print(
        f"Using shared sample cohort from {shared_sample_manifest_path} "
        f"({len(shared_sample_ids)} samples); retained "
        f"{phenotypes['Sample ID'].nunique()} samples in plot input."
    )

    # Collapse duplicate/ambiguous rows to one phenotype per sample and gene.
    resolved = resolve_entity_gene_phenotypes(phenotypes, entity_column="Sample ID")
    resolved["is_actionable"] = resolved.apply(
        lambda row: is_actionable_phenotype(row["Gene"], row["Phenotype"]),
        axis=1,
    )

    # The table includes represented samples only; missing gene calls become
    # False in Boolean columns rather than creating additional phenotype rows.
    sample_ids = sorted(resolved["Sample ID"].unique())
    actionable_table = pd.DataFrame({"Sample ID": sample_ids})

    table_data = []
    for gene in GENES:
        gene_data = resolved[resolved["Gene"] == gene].copy()
        actionable_sample_ids = set(gene_data.loc[gene_data["is_actionable"], "Sample ID"])
        actionable_table[gene] = actionable_table["Sample ID"].isin(actionable_sample_ids)

        # "Interesting" is the historical output name for actionable calls.
        # Denominator: resolved rows for this gene, including No Result.
        interesting_sample_ids = int(gene_data["is_actionable"].sum())
        total_sample_ids = gene_data["Sample ID"].nunique()
        fraction = interesting_sample_ids / total_sample_ids if total_sample_ids else 0
        table_data.append([gene, interesting_sample_ids, total_sample_ids, fraction])

    # Save the full-gene aggregate before the bar-only subset is considered.
    actionable_table["Aggregate"] = actionable_table[GENES].any(axis=1)
    actionable_table.to_csv(base_dir / "actionable_phenotypes_table.csv", index=False)
    print("Actionable phenotypes table saved as: actionable_phenotypes_table.csv")

    total_sample_ids = actionable_table["Sample ID"].nunique()
    aggregate_fraction = (
        actionable_table["Aggregate"].sum() / total_sample_ids if total_sample_ids else 0
    )
    table_data.append(
        ["Aggregate", int(actionable_table["Aggregate"].sum()), total_sample_ids, aggregate_fraction]
    )

    # Bar-chart aggregate uses INTERACTION_LIST. The saved sample CSV and the
    # table_data summary retain their earlier aggregate across all GENES.
    bar_results = [row for row in table_data if row[0] in INTERACTION_LIST]
    actionable_table["Aggregate"] = actionable_table[INTERACTION_LIST].any(axis=1)
    aggregate_interaction_count = int(actionable_table["Aggregate"].sum())
    aggregate_interaction_fraction = (
        aggregate_interaction_count / total_sample_ids if total_sample_ids else 0
    )
    bar_results.append(
        [
            "Aggregate",
            aggregate_interaction_count,
            total_sample_ids,
            aggregate_interaction_fraction,
        ]
    )
    bar_results_df = pd.DataFrame(
        bar_results,
        columns=[
            "Gene",
            "Interesting Phenotypes",
            "Total Sample IDs",
            "Fraction of Interesting Phenotypes",
        ],
    )
    print(
        "Aggregate bar: "
        f"{aggregate_interaction_count} / {total_sample_ids} = {aggregate_interaction_fraction:.3f}"
    )

    # Fractions and their complements form full-height stacked bars (0 to 1).
    fig, ax = plt.subplots(figsize=(10, 6))
    for _, row in bar_results_df.iterrows():
        gene = row["Gene"]
        fraction = row["Fraction of Interesting Phenotypes"]
        color = "darkblue" if gene == "Aggregate" else "skyblue"
        ax.bar(gene, fraction, color=color)
        ax.bar(gene, 1 - fraction, bottom=fraction, color="lightgray")
    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel("Fraction of Actionable Phenotypes", fontsize=12)
    ax.set_xticks(range(len(bar_results_df["Gene"])))
    ax.set_xticklabels(bar_results_df["Gene"], rotation=45, fontsize=10)
    plt.tight_layout()
    plt.savefig(
        base_dir / "actionable_phenotypes_fraction_with_aggregate_and_stacked.png",
        dpi=EXPORT_DPI,
    )
    plt.close()

    # One panel per gene in GENES order; the fixed grid holds at most 16 genes.
    fig, axes = plt.subplots(4, 4, figsize=(20, 19))
    axes = axes.flatten()
    for idx, gene in enumerate(GENES):
        gene_data = resolved[resolved["Gene"] == gene].copy()
        if gene_data.empty:
            axes[idx].axis("off")
            continue

        phenotype_counts = gene_data["Phenotype"].value_counts()
        labels = sorted(phenotype_counts.index.tolist(), key=lambda phenotype: phenotype_sort_key(gene, phenotype))
        phenotype_counts = phenotype_counts.reindex(labels)
        total_count = phenotype_counts.sum()
        # Drop categories below 0.05% before pie plotting. Matplotlib normalizes
        # percentages over the retained counts, not the original gene total.
        phenotype_counts = phenotype_counts[
            phenotype_counts.apply(lambda count: (count / total_count) * 100 >= 0.05)
        ]
        labels = phenotype_counts.index.tolist()
        color_map = build_color_map(gene, labels)
        colors = [color_map[label] for label in labels]

        wedges, _, autotexts = axes[idx].pie(
            phenotype_counts,
            labels=None,
            autopct=format_pct,
            startangle=140,
            colors=colors,
            textprops={"fontsize": 15},
            pctdistance=0.68,
            radius=1.24,
        )
        axes[idx].set_title(gene, fontsize=18)
        annotate_small_slices(axes[idx], wedges, autotexts)

        legend_handles = [Patch(facecolor=color_map[label], label=label) for label in labels]
        axes[idx].legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.14),
            fontsize=10,
            frameon=False,
            ncol=1 if len(labels) <= 3 else 2,
        )

    plt.tight_layout(pad=1.2)
    fig.subplots_adjust(left=0.04, right=0.98, top=0.97, bottom=0.05, wspace=-0.2, hspace=0.72)
    plt.savefig(
        base_dir / "all_phenotype_pie_charts.png",
        bbox_inches="tight",
        dpi=EXPORT_DPI,
    )
    plt.savefig(base_dir / "all_phenotype_pie_charts.svg", bbox_inches="tight")
    plt.close()

    # Full-gene summary: unaffected by narrowing the bar-only INTERACTION_LIST.
    pd.DataFrame(
        table_data,
        columns=[
            "Gene",
            "Interesting Phenotypes",
            "Total Sample IDs",
            "Fraction of Interesting Phenotypes",
        ],
    ).to_csv(base_dir / "phenotype_data_table.csv", index=False)

    print("Analysis complete. Outputs saved.")


if __name__ == "__main__":
    main()
