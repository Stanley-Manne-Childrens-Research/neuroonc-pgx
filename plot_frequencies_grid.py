"""Generate Figure 1 of the manuscript: allele frequencies across PGx genes.

Run from the repository root: python3 plot_frequencies_grid.py
Requires Python, pandas, numpy, matplotlib, and openpyxl for XLSX reading.
This is the manuscript Figure 1 generation script. It writes the 300-dpi PNG
DiplotypeFrequencies/all_genes_allele_frequencies_grid.png, overwriting an
existing file. It does not create or update frequency tables or Table S1.

Upstream files and preparation
-----------------------------
DiplotypeFrequencies/ must already contain the prepared reference workbooks.
The script reads files ending in frequency_table.xlsx, ignoring Excel lock
files beginning with "~$". Discovery is nonrecursive; a missing gene workbook
omits that gene. Gene labels come from the filename before the first underscore.
The directory path is relative to the working directory, not this script.

Each workbook's first worksheet must have headers on Excel row 2, allele labels
in column A, numeric reference-frequency columns, and an exact "CBTN" header.
Frequencies are proportions between 0 and 1, not percentages. The CBTN column
must already contain the intended cohort's observed frequencies: populate it
using get_cbtn_allele_frequencies.py for non-CYP2D6 genes and
update_cyp2d6_frequency_table.py for CYP2D6. The latter consumes the frequencies
from aggregate_stellarpgx_calls.py. Neither individual call files, the cohort
manifest, nor intermediate CSVs are read by this plotting script. It assumes
the workbooks are current and cohort-consistent; it does not verify provenance.

Figure 1 quantities and selection
--------------------------------
Each panel represents one gene. Grey lines with horizontal end marks show the
minimum and maximum positive frequencies among the selected workbook columns.
Open circles show observed CBTN frequencies; open diamonds show ancestry-based
calculated frequencies. The vertical axis is logarithmic. Legends are removed
from the saved figure, so the manuscript caption must define these symbols.

Calculated frequency = sum(ancestry-group proportion * reference frequency *
mapping multiplier). corbett_counts defines the fixed ancestry composition and
mapping connects reference-column names to those groups. Unmapped populations
do not contribute to the calculation; the full AMR weight is assigned to Latino.
Missing/nonnumeric selected values are converted to zero and nonpositive values
are dropped. Available ancestry weights are not renormalized for missing data.

Rows with at least one positive selected value are eligible. Within each gene,
sort by observed CBTN frequency descending (missing values last) and retain at
most ten rows. Draw markers/ranges only when the selected reference maximum is
at least 0.0001; tick labels are retained even for rows below this threshold.
G6PD display labels retain only the text before the first comma. Missing CBTN
values are not imputed, and zero/missing observed or calculated frequencies
cannot be displayed on the logarithmic axis.

Layout assumptions and current limitations
------------------------------------------
Reference columns are selected positionally with data.columns[1:9] (Excel B-I,
eight columns), not by an explicit population-name list. Verify workbook column
order before reproducing Figure 1. Any CBTN column inside that slice also enters
row eligibility and the min/max range, although it is ignored by the ancestry
mapping. Mapped reference columns outside the slice do not contribute.

Panels follow os.listdir order, without sorting; letter-to-gene assignments can
therefore vary with filesystem ordering. The fixed 4-by-4 layout supports at
most 16 included genes, with unused axes removed. Each included gene must have
at least one eligible row; empty ranges are not handled. Panel limits use the
selected min/max ranges, so observed/calculated points outside them may be
clipped. Missing CBTN headers are reported and their workbooks skipped; other
file/schema errors propagate. The output directory must be writable.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Ancestry-group counts defining the mixture used for calculated frequencies;
# these fixed inputs are not derived from the current cohort manifest.
# AFR = African; AMR = admixed American; EAS = East Asian;
# EUR = European; SAS = South Asian. Counts total 1,452 individuals.
corbett_counts =	{
                'AFR' : 153,
                'AMR' : 222,
                'EAS' : 67,
                'EUR' : 968,
                'SAS' : 42,
                }
# Ancestry mixture proportions: each count divided by the common total of 1,452.
# The denominator stays fixed for alleles with incomplete reference values.
corbett_freqs = dict([(k, float(v)/sum((corbett_counts.values()))) for k, v in corbett_counts.items()])

# Exact reference-column name -> (ancestry group, additional multiplier).
# Active multipliers are all 1.0; Latino receives the entire AMR mixture weight.
# The separate American column is intentionally not part of this mapping.
mapping = {
'African American/Afro-Caribbean' : ('AFR', 1.0),
#'American': ('AMR', 0.0),
'Central/South Asian': ('SAS', 1.0),
'East Asian': ('EAS', 1.0),
'European': ('EUR', 1.0),
'Latino': ('AMR', 1.0)
}

# Prepared reference/observed-frequency workbooks and Figure 1 output directory,
# relative to the working directory. No directories or workbooks are created.
directory_path = "DiplotypeFrequencies"  

# Per-gene panel data, accumulated before the manuscript Figure 1 grid is drawn.
plot_data = []

for file_name in os.listdir(directory_path):

    gene_name = file_name.split("_")[0]

    if file_name.endswith("frequency_table.xlsx") and not file_name.startswith("~$"):
        file_path = os.path.join(directory_path, file_name)
        
        # Read the first worksheet, treating Excel row 2 as the column headers.
        data = pd.read_excel(file_path, header=1)

        # Select Excel B-I positionally. Shorter/differently arranged tables may
        # place CBTN inside this slice; see the module's layout assumptions.
        biogeo_columns = data.columns[1:9]
        cbtn_column = "CBTN"  # Exact header for observed cohort frequencies.

        # Ensure the CBTN column is included
        if cbtn_column not in data.columns:
            print(f"The required column '{cbtn_column}' is not present in the file: {file_name}")
            continue

        # Prepare the data
        alleles = data.iloc[:, 0]  # First column (Star Alleles)
        ranges = []
        admixed = []
        cbtn_frequencies = []
        valid_indices = []

        for idx, row in data.iterrows():
            # Coerce mixed/object values to numeric before filling nulls to avoid pandas downcasting warnings.
            frequencies = pd.to_numeric(row[biogeo_columns], errors="coerce").fillna(0.0)
            frequencies = frequencies[frequencies > 0]  # Exclude zero values
            if not frequencies.empty:  # Only include rows with non-zero frequencies
                valid_indices.append(idx)
                cbtn_frequencies.append(row[cbtn_column])
                ranges.append((frequencies.min(), frequencies.max()))
                summer = 0.
                for description, frequency in (frequencies).items():
                    if description in mapping:
                        summer += (corbett_freqs[mapping[description][0]]*frequency*mapping[description][1])

                admixed.append(summer)

        # Filter alleles based on valid indices, then sort by descending CBTN frequency.
        alleles = alleles.iloc[valid_indices]
        sorted_indices = pd.Series(cbtn_frequencies).sort_values(ascending=False).index
        alleles = alleles.iloc[sorted_indices]
        ranges = [ranges[i] for i in sorted_indices]
        cbtn_frequencies = [cbtn_frequencies[i] for i in sorted_indices]
        admixed = [admixed[i] for i in sorted_indices]

        # Figure 1 displays at most ten eligible alleles per gene, ranked by CBTN.
        if len(alleles) > 10:
            alleles = alleles.iloc[:10]
            ranges = ranges[:10]
            cbtn_frequencies = cbtn_frequencies[:10]
            admixed = admixed[:10]
        # Store only the selected rows used in this gene's Figure 1 panel.
        plot_data.append({
            "gene_name": gene_name,
            "alleles": alleles,
            "ranges": ranges,
            "cbtn_frequencies": cbtn_frequencies,
            "admixed": admixed
        })

# Number of included workbooks/panels; the fixed dimensions below permit 16.
n_genes = len(plot_data)
# Manuscript Figure 1 layout: four panel columns and four panel rows.
ncols = 4
nrows = 4

# Figure dimensions are inches; panels are lettered in filesystem input order.
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(28, 28), squeeze=False)

# Panel identifiers used to reference individual genes in the manuscript caption.
panel_labels = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")[:n_genes]

for idx, pdata in enumerate(plot_data):
    row = idx // ncols
    col = idx % ncols
    ax = axes[row][col]

    alleles = pdata["alleles"]
    ranges = pdata["ranges"]
    cbtn_frequencies = pdata["cbtn_frequencies"]
    admixed = pdata["admixed"]
    gene_name = pdata["gene_name"]

    # Abbreviate G6PD labels by retaining only text before the first comma.
    if gene_name.upper() == "G6PD":
        print([str(a) for a in alleles])
        xtick_labels = [" ".join(str(a).split(',')[:1]) for a in alleles]
    else:
        xtick_labels = [str(a) for a in alleles]

    for i, (allele, (min_val, max_val), cbtn, adm) in enumerate(zip(alleles, ranges, cbtn_frequencies, admixed)):
        if max_val >= 0.0001:
            ax.plot([i, i], [min_val, max_val], color="grey", lw=4)  # Range
            # Use marker='_' for horizontal hash marks at min_val and max_val
            ax.scatter([i], [min_val], marker='_', color="grey", s=360, lw=4, zorder=6)
            ax.scatter([i], [max_val], marker='_', color="grey", s=360, lw=4, zorder=6)
            ax.scatter([i], [cbtn], edgecolor="black", facecolor="none", s=220, linewidths=2.5, zorder=5, label="CBTN Cohort" if i == 0 else "")
            ax.scatter([i], [adm], marker="D", edgecolor="black", facecolor="none", s=200, linewidths=2.5, zorder=5, label="Calculated based on ancestry" if i == 0 else "")

    ax.set_xticks(range(len(alleles)))
    ax.set_xticklabels(xtick_labels, rotation=45, ha="right", fontsize=16)
    ax.set_yscale("log")
    ax.set_ylim(bottom=min(r[0] for r in ranges) * 0.9, top=max(r[1] for r in ranges) * 1.1)
    ax.set_ylabel("Frequency", fontsize=20)
    ax.tick_params(axis="y", labelsize=16, width=1.8)
    ax.tick_params(axis="x", width=1.8, labelsize=16)
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)

    # Manuscript Figure 1 panel letter outside the upper left corner.
    ax.annotate(f"{panel_labels[idx]}", xy=(0, 1.05 ), xycoords='axes fraction',
                xytext=(-30, 18), textcoords='offset points',
                fontsize=24, fontweight='bold', va='top', ha='left')

    # Gene name inside, upper right
    ax.text(0.98, 0.92, f"{gene_name}", transform=ax.transAxes,
            fontsize=20, fontweight='bold', va='top', ha='right',
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # Remove panel legends; marker/range definitions belong in the Figure 1 caption.
    if idx == 0:
        leg = ax.legend(fontsize=16, frameon=False)
        if leg:
            leg.remove()
    else:
        ax.legend().set_visible(False)

# Hide unused subplots
for idx in range(n_genes, nrows * ncols):
    row = idx // ncols
    col = idx % ncols
    fig.delaxes(axes[row][col])

plt.tight_layout(rect=(0, 0.03, 1, 1), h_pad=2.0)
fig.subplots_adjust(hspace=0.45, bottom=0.10)
# Manuscript Figure 1 artifact: a 300-dpi PNG, overwritten when the script runs.
output_path = os.path.join(directory_path, "all_genes_allele_frequencies_grid.png")
plt.savefig(output_path, dpi=300)
plt.close(fig)
print(f"Grid plot saved at {output_path}")
