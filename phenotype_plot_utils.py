from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable

import pandas as pd


GENES = [
    "ABCG2",
    "CACNA1S",
    "CYP2B6",
    "CYP2C9",
    "CYP2C19",
    "CYP2D6",
    "CYP3A5",
    "CYP4F2",
    "DPYD",
    "G6PD",
    "NUDT15",
    "RYR1", 
    "SLCO1B1",
    "TPMT",
    "UGT1A1",
    "VKORC1",

]

METABOLIZER_STATUS_MAP = {
    "Normal metaboliser (NM)": "Normal Metabolizer",
    "Intermediate metaboliser (IM)": "Intermediate Metabolizer",
    "Poor metaboliser (PM)": "Poor Metabolizer",
    "Ultrarapid metaboliser (UM)": "Ultrarapid Metabolizer",
    "Indeterminate": "Indeterminate",
}

ACTIONABLE_PHENOTYPES = {
    "CYP2C19": {
        "Indeterminate",
        "Intermediate Metabolizer",
        "Likely Intermediate Metabolizer",
        "Likely Poor Metabolizer",
        "Poor Metabolizer",
        "Rapid Metabolizer",
        "Ultrarapid Metabolizer",
    },
    "CYP2C9": {"Indeterminate", "Intermediate Metabolizer", "Poor Metabolizer"},
    "ABCG2": {"Decreased Function", "Poor Function"},
    "SLCO1B1": {
        "Decreased Function",
        "Poor Function",
        "Possible Decreased Function",
    },
    "NUDT15": {
        "Indeterminate",
        "Intermediate Metabolizer",
        "Poor Metabolizer",
        "Possible Intermediate Metabolizer",
    },
    "CYP3A5": {"Intermediate Metabolizer", "Normal Metabolizer"},
    "CYP2B6": {
        "Indeterminate",
        "Intermediate Metabolizer",
        "Poor Metabolizer",
    },
    "G6PD": {"Deficient", "Indeterminate", "Variable"},
    "CYP4F2": {"*1/*3", "*3/*3"},
    "TPMT": {"Indeterminate", "Intermediate Metabolizer", "Poor Metabolizer"},
    "DPYD": {"Intermediate Metabolizer"},
    "CYP2D6": {
        "Intermediate Metabolizer",
        "Poor Metabolizer",
        "Ultrarapid Metabolizer",
    },
    "UGT1A1": {"Poor Metabolizer" },
    "VKORC1": {"-1639 AA", "-1639 AG"},
    "CACNA1S": {"Malignant Hyperthermia Susceptibility"},
    "RYR1": {"Malignant Hyperthermia Susceptibility"},
}

NON_ACTIONABLE_PHENOTYPES = {
    "CYP2C19": {"No Result", "Normal", "Normal Metabolizer"},
    "CYP2C9": {"No Result", "Normal", "Normal Metabolizer"},
    "ABCG2": {"No Result", "Normal Function"},
    "SLCO1B1": {"No Result", "Increased Function", "Indeterminate", "Normal Function"},
    "NUDT15": {"No Result", "Normal", "Normal Metabolizer"},
    "CYP3A5": {"No Result", "Poor Metabolizer"},
    "CYP2B6": {
        "No Result",
        "Normal Metabolizer",
        "Rapid Metabolizer",
        "Ultrarapid Metabolizer",
    },
    "G6PD": {"No Result", "Normal", "Normal Metabolizer"},
    "CYP4F2": {"No Result", "Normal Metabolizer"},
    "TPMT": {"No Result", "Normal", "Normal Metabolizer"},
    "DPYD": {"No Result", "Normal", "Normal Metabolizer"},
    "CYP2D6": {"Indeterminate", "No Result", "Normal", "Normal Metabolizer"},
    "UGT1A1": {
        "Indeterminate",
        "Intermediate Metabolizer",
        "No Result",
        "Normal Metabolizer",
    },
        "VKORC1": {"-1639 GG"},
    "CACNA1S": {"Uncertain Susceptibility"},
    "RYR1": {"Uncertain Susceptibility"},

}

RESOLUTION_PRIORITIES = {
    "CYP2B6": [
        "Poor Metabolizer",
        "Indeterminate",
        "Intermediate Metabolizer",
        "Ultrarapid Metabolizer",
        "Rapid Metabolizer",
        "Normal Metabolizer",
        "No Result",
    ],
    "TPMT": [
        "Poor Metabolizer",
        "Intermediate Metabolizer",
        "Indeterminate",
        "Normal Metabolizer",
        "No Result",
    ],
    "NUDT15": [
        "Poor Metabolizer",
        "Intermediate Metabolizer",
        "Possible Intermediate Metabolizer",
        "Indeterminate",
        "Normal Metabolizer",
        "No Result",
    ],
    "SLCO1B1": [
        "Poor Function",
        "Decreased Function",
        "Possible Decreased Function",
        "Indeterminate",
        "Normal Function",
        "Increased Function",
        "No Result",
    ],
    "CYP4F2": [
        "*3/*3",
        "*1/*3",
        "Normal Metabolizer",
        "No Result",
    ],
}

_WARNED_UNCATEGORIZED: set[tuple[str, str]] = set()
_WARNED_CONFLICTS: set[tuple[str, tuple[str, ...]]] = set()


def normalize_phenotype(value: object) -> str:
    if pd.isna(value):
        return "No Result"
    phenotype = str(value).strip()
    if not phenotype:
        return "No Result"
    if phenotype.lower() in {"n/a", "na", "nan", "none", "no result"}:
        return "No Result"
    return phenotype


def classify_phenotype(gene: str, phenotype: object) -> str:
    normalized = normalize_phenotype(phenotype)
    if normalized in ACTIONABLE_PHENOTYPES.get(gene, set()):
        return "actionable"
    if normalized in NON_ACTIONABLE_PHENOTYPES.get(gene, set()):
        return "non_actionable"

    warning_key = (gene, normalized)
    if warning_key not in _WARNED_UNCATEGORIZED:
        print(
            f"Warning: uncategorized phenotype for {gene}: {normalized!r}; "
            "treating as non-actionable."
        )
        _WARNED_UNCATEGORIZED.add(warning_key)
    return "non_actionable"


def is_actionable_phenotype(gene: str, phenotype: object) -> bool:
    return classify_phenotype(gene, phenotype) == "actionable"


def _cyp4f2_phenotype_from_diplotype(value: object) -> str:
    """Map CYP4F2 calls to the study's diplotype-based phenotype categories."""
    diplotype = normalize_phenotype(value)
    alleles = [allele.strip() for allele in diplotype.split("/")]
    if len(alleles) != 2 or any(
        re.fullmatch(r"\*\d+(?:\.\d+)?", allele) is None for allele in alleles
    ):
        return "No Result"
    if sorted(alleles) == ["*1", "*3"]:
        return "*1/*3"
    if alleles == ["*3", "*3"]:
        return "*3/*3"
    return "Normal Metabolizer"


def load_cpic_phenotypes(
    workbook_path: str | Path,
    *,
    genes: Iterable[str] = GENES,
    cyp2d6_calls_path: str | Path | None = None,
) -> pd.DataFrame:
    """Load DPWG phenotypes for VKORC1 and CPIC phenotypes for other genes.

    The historical function name is retained for compatibility. VKORC1 rows
    are selected only from DPWG, without a fallback to CPIC.
    CYP4F2 phenotypes are derived from CPIC diplotypes using the study's
    categories rather than the source phenotype annotations.
    """
    print(f"Loading phenotypes from {workbook_path} (DPWG for VKORC1; CPIC for other genes)...")
    data = pd.read_csv(workbook_path, header=None, usecols=[5, 0, 6, 7, 8], dtype=str,     encoding="utf-8-sig")
    data = data.rename(columns={5: "Authority", 0: "Sample ID", 6: "Gene", 8: "Phenotype", 7: "Diplotype"})
    is_vkorc1 = data["Gene"] == "VKORC1"
    selected = data[
        (is_vkorc1 & (data["Authority"] == "DPWG"))
        | (~is_vkorc1 & (data["Authority"] == "CPIC"))
    ].copy()
    selected["Sample ID"] = selected["Sample ID"].astype(str)
    selected["Gene"] = selected["Gene"].astype(str)

    if cyp2d6_calls_path is not None:
        selected = _replace_cyp2d6_calls(selected, Path(cyp2d6_calls_path), data.columns.tolist())

    selected = selected[selected["Gene"].isin(list(genes))].copy()
    is_cyp4f2 = selected["Gene"] == "CYP4F2"
    selected.loc[is_cyp4f2, "Phenotype"] = selected.loc[
        is_cyp4f2, "Diplotype"
    ].apply(_cyp4f2_phenotype_from_diplotype)
    selected["Phenotype"] = selected["Phenotype"].apply(normalize_phenotype)
    return selected


def resolve_entity_gene_phenotypes(
    data: pd.DataFrame,
    *,
    entity_column: str,
) -> pd.DataFrame:
    resolved_input = data[[entity_column, "Gene", "Phenotype"]].copy()
    resolved_input[entity_column] = resolved_input[entity_column].astype(str)
    resolved_input["Gene"] = resolved_input["Gene"].astype(str)
    resolved_input["Phenotype"] = resolved_input["Phenotype"].apply(normalize_phenotype)
    resolved_input["_row_order"] = range(len(resolved_input))
    resolved_input = resolved_input.drop_duplicates(
        subset=[entity_column, "Gene", "Phenotype"],
        keep="first",
    )

    resolved_rows = []
    for (entity_id, gene), group in resolved_input.groupby([entity_column, "Gene"], sort=True):
        ordered_group = group.sort_values("_row_order")
        phenotypes = ordered_group["Phenotype"].tolist()
        unique_phenotypes = list(dict.fromkeys(phenotypes))
        resolved_rows.append(
            {
                entity_column: entity_id,
                "Gene": gene,
                "Phenotype": _resolve_conflicting_phenotypes(gene, unique_phenotypes),
            }
        )

    return pd.DataFrame(resolved_rows)


def _replace_cyp2d6_calls(
    data_cpic: pd.DataFrame,
    cyp2d6_calls_path: Path,
    data_columns: list[str],
) -> pd.DataFrame:
    cyp2d6_existing = data_cpic[data_cpic["Gene"] == "CYP2D6"].copy()
    if cyp2d6_existing.empty:
        return data_cpic

    cyp2d6_calls = pd.read_csv(cyp2d6_calls_path).copy()
    cyp2d6_calls = cyp2d6_calls.drop_duplicates(subset=["sample_id"]).copy()
    cyp2d6_calls["Phenotype"] = (
        cyp2d6_calls["metabolizer_status"]
        .map(METABOLIZER_STATUS_MAP)
        .fillna(cyp2d6_calls["metabolizer_status"])
        .apply(normalize_phenotype)
    )
    cyp2d6_calls["sample_id"] = cyp2d6_calls["sample_id"].astype(str)

    cyp2d6_existing = cyp2d6_existing.merge(
        cyp2d6_calls[["sample_id", "result", "Phenotype"]],
        left_on="Sample ID",
        right_on="sample_id",
        how="left",
    )
    cyp2d6_existing["Diplotype"] = cyp2d6_existing["result"].fillna(cyp2d6_existing["Diplotype"])
    cyp2d6_existing["Phenotype"] = cyp2d6_existing["Phenotype_y"].fillna(
        cyp2d6_existing["Phenotype_x"]
    )
    cyp2d6_existing = cyp2d6_existing[data_columns].copy()

    data_without_cyp2d6 = data_cpic[data_cpic["Gene"] != "CYP2D6"].copy()
    return pd.concat([data_without_cyp2d6, cyp2d6_existing], ignore_index=True)


def _resolve_conflicting_phenotypes(gene: str, phenotypes: list[str]) -> str:
    if len(phenotypes) == 1:
        return phenotypes[0]

    priorities = RESOLUTION_PRIORITIES.get(gene)
    if priorities:
        rank = {phenotype: idx for idx, phenotype in enumerate(priorities)}
        ranked = [phenotype for phenotype in phenotypes if phenotype in rank]
        if ranked:
            return min(ranked, key=lambda phenotype: rank[phenotype])

    conflict_key = (gene, tuple(phenotypes))
    if conflict_key not in _WARNED_CONFLICTS:
        print(
            f"Warning: unresolved phenotype conflict for {gene}: {phenotypes}; "
            f"keeping first value {phenotypes[0]!r}."
        )
        _WARNED_CONFLICTS.add(conflict_key)
    return phenotypes[0]
