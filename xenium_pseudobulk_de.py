"""Xenium pseudobulk differential expression: invasion versus lepidic.

Tested API: PyDESeq2 0.5.4. Use unnormalized, nonnegative integer counts.
Each sample_key must identify an independent biological unit, usually a patient
with one tumor. Repeated ROIs/sections from that unit are summed, not replicated.
The input AnnData is not modified.
"""

from __future__ import annotations

import json
import platform
import warnings
from datetime import datetime, timezone
from importlib.metadata import version
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


TYPES = ["lepidic", "invasion"]
COLORS = {"lepidic": "#01B0F1", "invasion": "#FE0000"}


def _sample_sort_dataframe(sample_ids):
    """Return sortable metadata for sample IDs using the numeric part after tkbadeno.

    Unknown formats are sent to the end and then sorted lexicographically by sample_id.
    """
    sample_ids = pd.Index(sample_ids, dtype="object")
    sample_text = pd.Series(sample_ids.astype(str), index=sample_ids, name="sample_id")
    case_number = pd.to_numeric(
        sample_text.str.extract(r"(?i)tkbadeno(\d+)", expand=False),
        errors="coerce",
    )
    return pd.DataFrame({
        "sample_id": sample_text,
        "_case_number": case_number.fillna(np.inf),
        "_sample_label": sample_text,
    }, index=sample_ids)


def _get_counts(adata, counts_source):
    if counts_source == "X":
        return adata.X, adata.var
    if counts_source == "raw":
        if adata.raw is None:
            raise ValueError("adata.raw is None. Select the actual raw-count source.")
        return adata.raw.X, adata.raw.var
    if counts_source not in adata.layers:
        raise KeyError(
            f"No layer {counts_source!r}. Available layers: {list(adata.layers.keys())}. "
            "Use counts_source='X' or 'raw' only if its provenance is raw counts."
        )
    return adata.layers[counts_source], adata.var


def _validate_integer_counts(matrix):
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    # Chunked checks avoid creating a second large dense array.
    for start in range(0, values.size, 1_000_000):
        block = values[start:start + 1_000_000]
        if not np.isfinite(block).all() or (block < 0).any():
            raise ValueError("Counts contain NaN, infinity, or negative values.")
        if not np.equal(block, np.floor(block)).all():
            raise ValueError(
                "Non-integer values detected. Do not round normalized/log/expected "
                "counts to force DESeq2. Recover an appropriate integer count matrix."
            )


def make_pseudobulk(
    adata,
    *,
    sample_key,
    type_key="type",
    counts_source="counts",
    paired=True,
    min_cells=50,
    min_count=10,
    min_samples=3,
    exclude_genes=None,
):
    """Sum counts for each biological sample x histology.

    counts_source: a layer name, 'X', or 'raw'; no automatic fallback is used.
    paired=True: keep complete pairs after cell-count/total-count QC.
    paired=False: each sample must belong to exactly one histology before QC.
    min_cells is a user-chosen QC setting, not an established Xenium standard.
    A gene needs >= min_count counts in >= min_samples retained pseudobulks.
    This expression filter does not require expression in BOTH histologies.
    """
    for key in [sample_key, type_key]:
        if key not in adata.obs:
            raise KeyError(f"adata.obs has no {key!r}. Columns: {list(adata.obs.columns)}")
    if sample_key == type_key:
        raise ValueError("sample_key and type_key must be different.")
    for key, value in {"min_cells": min_cells, "min_count": min_count,
                       "min_samples": min_samples}.items():
        if not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{key} must be a positive integer.")

    selected = adata.obs[type_key].isin(TYPES).to_numpy()
    if not selected.any():
        raise ValueError("No cells labeled exactly 'lepidic' or 'invasion'.")
    obs = adata.obs.loc[selected, [sample_key, type_key]].copy()
    obs.columns = ["sample_id", "type"]
    if obs["sample_id"].isna().any():
        raise ValueError("Selected cells have missing biological sample IDs.")
    original_n = obs["sample_id"].nunique()
    obs["sample_id"] = obs["sample_id"].astype(str)
    obs["type"] = obs["type"].astype(str)
    if obs["sample_id"].nunique() != original_n:
        raise ValueError("Different sample IDs collide when converted to strings.")
    if obs["sample_id"].str.strip().eq("").any():
        raise ValueError("Empty sample IDs are not allowed.")
    if not paired and (obs.groupby("sample_id")["type"].nunique() > 1).any():
        raise ValueError(
            "paired=False requires independent groups: each biological sample must "
            "have only one type. For overlapping samples, use paired=True or a "
            "separately justified repeated-measures design."
        )

    matrix, var = _get_counts(adata, counts_source)
    if not var.index.is_unique:
        raise ValueError("Gene identifiers must be unique; resolve duplicates upstream.")
    if not (sparse.issparse(matrix) or isinstance(matrix, np.ndarray)):
        raise TypeError(
            "Use an in-memory CPU NumPy/SciPy count matrix. Convert GPU/backed "
            "counts to CPU/in-memory before this function."
        )
    gene_names = pd.Index(var.index.astype(str), name="gene")
    if not gene_names.is_unique:
        raise ValueError("Gene identifiers collide after conversion to strings.")
    # Xenium controls are not biological genes. Prefer explicit feature metadata.
    eligible_gene = np.ones(len(gene_names), dtype=bool)
    if "feature_types" in var.columns:
        feature_types = var["feature_types"].astype(str)
        if feature_types.eq("Gene Expression").any():
            eligible_gene &= feature_types.eq("Gene Expression").to_numpy()
        else:
            warnings.warn("Unrecognized feature_types; inspect controls manually.")
    eligible_gene &= ~gene_names.str.match(
        r"^(NegControlProbe_|NegControlCodeword_|UnassignedCodeword_|DeprecatedCodeword_)"
    )
    if exclude_genes is not None:
        eligible_gene &= ~gene_names.isin(list(exclude_genes))
    if not eligible_gene.any():
        raise ValueError("No biological genes remain after feature exclusion.")
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    matrix = matrix[selected, :][:, eligible_gene]
    _validate_integer_counts(matrix)

    # Sparse one-hot aggregation: never densify the full cell x gene matrix.
    grouping = pd.MultiIndex.from_frame(obs[["sample_id", "type"]])
    codes, groups = pd.factorize(grouping, sort=True)
    membership = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.int64), (codes, np.arange(len(codes)))),
        shape=(len(groups), len(codes)),
    )
    summed = membership @ matrix
    if sparse.issparse(summed):
        summed = summed.toarray()  # Only pseudobulk x gene, usually a small matrix.
    summed = np.asarray(summed, dtype=np.int64)
    meta = groups.to_frame(index=False)
    meta.columns = ["sample_id", "type"]
    meta.index = pd.Index([f"pb{i:04d}" for i in range(len(meta))], name="pseudobulk_id")
    meta["n_cells"] = np.bincount(codes, minlength=len(groups))
    meta["total_counts"] = summed.sum(axis=1)
    counts = pd.DataFrame(summed, index=meta.index, columns=gene_names[eligible_gene])

    qc = meta.copy()
    qc["passes_cell_qc"] = qc["n_cells"].ge(min_cells)
    qc["passes_count_qc"] = qc["total_counts"].gt(0)
    valid = qc["passes_cell_qc"] & qc["passes_count_qc"]
    complete = pd.Series(True, index=qc.index)
    if paired:
        number_of_types = qc.loc[valid].groupby("sample_id")["type"].nunique()
        complete_ids = number_of_types.index[number_of_types.eq(2)]
        complete = qc["sample_id"].isin(complete_ids)
    qc["included"] = valid & complete
    qc["reason"] = np.select(
        [~qc["passes_cell_qc"], ~qc["passes_count_qc"], ~complete],
        ["too_few_cells", "zero_total_counts", "incomplete_pair"],
        default="included",
    )
    counts = counts.loc[qc["included"]].copy()
    meta = meta.loc[counts.index].copy()
    meta["sample_id"] = pd.Categorical(meta["sample_id"])
    meta["type"] = pd.Categorical(meta["type"], categories=TYPES, ordered=True)

    gene_qc = pd.DataFrame(index=gene_names)
    gene_qc["eligible_feature"] = eligible_gene
    gene_qc["n_pseudobulks_ge_min_count"] = 0
    detected = counts.ge(min_count).sum(axis=0)
    gene_qc.loc[detected.index, "n_pseudobulks_ge_min_count"] = detected
    gene_qc["included"] = (
        gene_qc["eligible_feature"] & gene_qc["n_pseudobulks_ge_min_count"].ge(min_samples)
    )
    gene_qc["reason"] = np.select(
        [~gene_qc["eligible_feature"], ~gene_qc["included"]],
        ["non_gene_or_excluded", "low_count"], default="included",
    )
    counts = counts.loc[:, gene_qc.index[gene_qc["included"]]].copy()
    settings = {
        "sample_key": sample_key, "type_key": type_key, "counts_source": counts_source,
        "paired": bool(paired), "min_cells": int(min_cells), "min_count": int(min_count),
        "min_samples": int(min_samples), "n_input_cells": int(adata.n_obs),
        "n_target_cells_before_qc": int(selected.sum()),
        "n_pseudobulks": int(len(meta)), "n_samples": int(meta["sample_id"].nunique()),
        "n_genes_before_filter": int(len(gene_names)), "n_genes_tested": int(counts.shape[1]),
        "exclude_genes": [] if exclude_genes is None else list(exclude_genes),
    }
    print(f"Pseudobulk: {len(meta)} profiles; {meta['sample_id'].nunique()} biological "
          f"samples; {counts.shape[1]} retained genes.")
    print(f"Excluded observed profiles: {int((~qc['included']).sum())} (see sample_qc).")
    return {"counts": counts, "metadata": meta, "sample_qc": qc,
            "gene_qc": gene_qc, "settings": settings}


def _paired_differences(log_counts, metadata):
    frames = {}
    for label in TYPES:
        rows = metadata.index[metadata["type"].eq(label)]
        frame = log_counts.loc[rows].copy()
        frame.index = pd.Index(metadata.loc[rows, "sample_id"].astype(str), name="sample_id")
        frames[label] = frame
    return frames["invasion"].subtract(frames["lepidic"]).sort_index()


def fit_pydeseq2(
    pb,
    *,
    outdir,
    n_cpus=4,
    alpha=0.05,
    lfc_cutoff=0.5,
    size_factors_fit_type="poscounts",
    control_genes=None,
):
    """Fit a Wald test of LFC=0, then apply FDR and a descriptive LFC cutoff.

    Positive log2FoldChange = invasion / lepidic > 1. LFCs are unshrunk MLEs.
    The cutoff is NOT a formal test that |LFC| exceeds the cutoff.
    At least 3 pairs (or 3 independent samples per group) are required here as
    a practical guard; this is not a claim of sufficient statistical power.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    if not 0 < alpha < 1 or lfc_cutoff < 0:
        raise ValueError("Require 0 < alpha < 1 and lfc_cutoff >= 0.")
    if n_cpus < 1:
        raise ValueError("n_cpus must be >= 1.")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    counts, meta = pb["counts"], pb["metadata"].copy()
    pb["sample_qc"].to_csv(outdir / "sample_qc.csv")
    pb["gene_qc"].to_csv(outdir / "gene_qc.csv")
    counts.to_csv(outdir / "pseudobulk_counts.csv")
    meta.to_csv(outdir / "pseudobulk_metadata.csv")
    number_per_type = meta["type"].value_counts().reindex(TYPES, fill_value=0)
    if number_per_type.min() < 3:
        raise ValueError(
            f"Too few biological replicates after QC: {number_per_type.to_dict()}. "
            "This workflow requires >=3 pairs or >=3 independent samples per group. "
            f"Inspect QC tables in {outdir}. Never create replicates by splitting cells."
        )
    if counts.shape[1] < 2:
        raise ValueError("Fewer than two genes pass filtering; inspect gene_qc.csv.")
    if (counts.sum(axis=1) == 0).any():
        raise ValueError("A pseudobulk has no counts after gene filtering. Revise QC.")
    if not counts.index.equals(meta.index):
        raise ValueError("Pseudobulk counts and metadata are misaligned.")
    if control_genes is not None:
        control_genes = list(control_genes)
        missing = sorted(set(control_genes) - set(counts.columns))
        if not control_genes or missing:
            raise ValueError(f"Control genes must be nonempty and pass filtering. Missing: {missing}")
    paired = pb["settings"]["paired"]
    design = "~sample_id + type" if paired else "~type"
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        counts=counts, metadata=meta[["sample_id", "type"]], design=design,
        size_factors_fit_type=size_factors_fit_type, control_genes=control_genes,
        refit_cooks=True, inference=inference, quiet=True,
    )
    design_matrix = np.asarray(dds.obsm["design_matrix"], dtype=float)
    rank = np.linalg.matrix_rank(design_matrix)
    residual_df = design_matrix.shape[0] - rank
    if rank < design_matrix.shape[1] or residual_df < 1:
        raise ValueError("Design is rank-deficient or has no residual degrees of freedom.")
    if residual_df < 3:
        warnings.warn("Fewer than 3 residual degrees of freedom: dispersion estimates are imprecise.")
    dds.deseq2()
    stats = DeseqStats(
        dds, contrast=["type", "invasion", "lepidic"], alpha=alpha,
        cooks_filter=True, independent_filter=True, inference=inference, quiet=True,
    )
    stats.summary()
    res = stats.results_df.copy()
    res.index.name = "gene"
    normalized = pd.DataFrame(dds.layers["normed_counts"], index=counts.index, columns=counts.columns)
    log_counts = np.log2(normalized + 1)
    for label in TYPES:
        res[f"mean_normalized_{label}"] = normalized.loc[meta["type"].eq(label)].mean()
    differences = None
    if paired:
        differences = _paired_differences(log_counts, meta)
        res["n_pairs"] = len(differences)
        res["fraction_pairs_invasion_higher"] = differences.gt(0).mean()
        res["median_paired_log2_difference"] = differences.median()
        res["fraction_pairs_matching_lfc_direction"] = np.where(
            res["log2FoldChange"] > 0,
            differences.gt(0).mean().reindex(res.index),
            differences.lt(0).mean().reindex(res.index),
        )
        differences.to_csv(outdir / "paired_log2_differences.csv")
    else:
        (outdir / "paired_log2_differences.csv").unlink(missing_ok=True)
    res["passes_fdr"] = res["padj"].lt(alpha)
    res["marker_group"] = np.select(
        [res["passes_fdr"] & res["log2FoldChange"].gt(0) & res["log2FoldChange"].ge(lfc_cutoff),
         res["passes_fdr"] & res["log2FoldChange"].lt(0) & res["log2FoldChange"].le(-lfc_cutoff)],
        ["invasion", "lepidic"], default="not_selected",
    )
    res = res.sort_values(["padj", "pvalue"], na_position="last")
    invasion = res.loc[res["marker_group"].eq("invasion")].copy()
    lepidic = res.loc[res["marker_group"].eq("lepidic")].copy()
    for name, table in {
        "deseq2_all_genes": res, "deg_fdr": res.loc[res["passes_fdr"]],
        "markers_invasion": invasion, "markers_lepidic": lepidic,
        "normalized_counts": normalized, "log2_normalized_counts": log_counts,
    }.items():
        table.to_csv(outdir / f"{name}.csv")
    meta["size_factor"] = np.asarray(dds.obs["size_factors"])
    meta.to_csv(outdir / "pseudobulk_metadata.csv")
    settings = {
        **pb["settings"], "design": design, "contrast": ["type", "invasion", "lepidic"],
        "alpha": float(alpha), "lfc_cutoff": float(lfc_cutoff),
        "lfc_estimate": "unshrunk MLE", "test_null": "log2FoldChange = 0",
        "size_factors_fit_type": size_factors_fit_type, "control_genes": control_genes,
        "n_cpus": int(n_cpus), "residual_df": int(residual_df),
        "n_invasion_markers": len(invasion), "n_lepidic_markers": len(lepidic),
        "n_genes_without_padj": int(res["padj"].isna().sum()),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "packages": {p: version(p) for p in ["pydeseq2", "anndata", "numpy", "pandas", "scipy", "matplotlib", "seaborn", "scikit-learn"]},
    }
    (outdir / "analysis_config.json").write_text(json.dumps(settings, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Design: {design}; invasion / lepidic contrast.")
    print(f"FDR < {alpha}: {int(res['passes_fdr'].sum())} genes; "
          f"with |log2FC| >= {lfc_cutoff}: invasion {len(invasion)}, lepidic {len(lepidic)}.")
    return {"results": res, "markers_invasion": invasion, "markers_lepidic": lepidic,
            "normalized_counts": normalized, "log2_normalized_counts": log_counts,
            "paired_log2_differences": differences, "dds": dds, "stats": stats,
            "metadata": meta, "pseudobulk": pb, "settings": settings, "outdir": outdir}


def plot_results(result, *, n_top=20, show=True):
    """Save plots and the matrices used to draw them.

    Each point/column is a pseudobulk biological sample, not a single cell.
    Plots use log2(size-factor-normalized counts + 1), not model residuals.
    No marker heatmaps/profiles are produced when no genes meet both thresholds.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    from sklearn.decomposition import PCA

    if not isinstance(n_top, (int, np.integer)) or n_top < 1:
        raise ValueError("n_top must be a positive integer.")
    outdir = Path(result["outdir"])
    res, meta = result["results"], result["metadata"]
    log_counts = result["log2_normalized_counts"]
    settings = result["settings"]
    paired = settings["paired"]
    files = {}

    def finish(fig, name):
        paths = []
        for ext in ["pdf", "tiff"]:
            path = outdir / f"{name}.{ext}"
            # Render first, then write the complete bytes to the destination.
            with BytesIO() as buffer:
                fig.savefig(buffer, format=ext, dpi=600, bbox_inches="tight")
                payload = buffer.getvalue()
            if ext == "pdf" and not payload.rstrip().endswith(b"%%EOF"):
                raise RuntimeError(f"Incomplete PDF rendering: {name}")
            path.write_bytes(payload)
            paths.append(str(path))
        files[name] = paths
        if show:
            plt.show()
        plt.close(fig)

    with plt.rc_context({"pdf.fonttype": 42, "font.family": "DejaVu Sans", "font.size": 10}):
        # QC includes excluded observed profiles, so missing/low-count pairs remain visible.
        qc = result["pseudobulk"]["sample_qc"]

        sample_order = _sample_sort_dataframe(
            qc["sample_id"].drop_duplicates()
        )
        ordered_sample_ids = (
            sample_order
            .sort_values(["_case_number", "_sample_label"])
            ["sample_id"]
            .tolist()
        )

        qc_plots = [
            (
                "n_cells",
                "Cells per pseudobulk",
                "cells_per_pseudobulk",
            ),
            (
                "total_counts",
                "Raw counts per pseudobulk",
                "raw_counts_per_pseudobulk",
            ),
        ]

        for variable, title, filename in qc_plots:
            fig, ax = plt.subplots(
                figsize=(max(8, 0.35 * len(ordered_sample_ids)), 4.5)
            )

            pivot = (
                qc.pivot(
                    index="sample_id",
                    columns="type",
                    values=variable,
                )
                .reindex(index=ordered_sample_ids, columns=TYPES)
                .fillna(0)
            )

            pivot.plot.bar(
                ax=ax,
                color=[COLORS[t] for t in TYPES],
                width=0.8,
            )

            ax.set_title(title)
            ax.set_ylabel(variable)
            ax.set_xlabel(
                "Biological sample ID (missing profiles shown as zero)"
            )
            ax.legend(title=None, frameon=False)

            if variable == "n_cells":
                ax.axhline(
                    settings["min_cells"],
                    color="0.4",
                    linestyle="--",
                    linewidth=1,
                )

            fig.tight_layout()

            finish(fig, filename)

        # PCA genes selected by variance, not differential-expression significance.
        variances = log_counts.var(axis=0)
        pca_genes = variances[variances.gt(0)].nlargest(min(500, len(variances))).index
        if len(pca_genes) >= 2:
            pca = PCA(n_components=2, svd_solver="full")
            xy = pca.fit_transform(log_counts[pca_genes])
            scores = pd.DataFrame(xy, index=meta.index, columns=["PC1", "PC2"]).join(meta)
            scores.to_csv(outdir / "pca_scores.csv")
            pd.Series(pca_genes, name="gene").to_csv(outdir / "pca_genes.csv", index=False)
            fig, ax = plt.subplots(figsize=(7, 6))
            if paired:
                for _, group in scores.groupby("sample_id", observed=True):
                    ax.plot(group["PC1"], group["PC2"], color="0.8", linewidth=0.8, zorder=1)
            for label in TYPES:
                group = scores.loc[scores["type"].eq(label)]
                ax.scatter(group["PC1"], group["PC2"], color=COLORS[label], label=label, s=45, zorder=2)
            if len(meta) <= 24:
                for _, row in scores.iterrows():
                    ax.annotate(str(row["sample_id"]), (row["PC1"], row["PC2"]),
                                xytext=(3, 3), textcoords="offset points", fontsize=7)
            ax.set(xlabel=f"PC1 ({pca.explained_variance_ratio_[0]:.1%})",
                   ylabel=f"PC2 ({pca.explained_variance_ratio_[1]:.1%})",
                   title=f"Pseudobulk PCA: {len(pca_genes)} variable genes\nlog2(normalized counts + 1); sample effects retained")
            ax.legend(frameon=False)
            sns.despine(ax=ax)
            fig.tight_layout()
            finish(fig, "pca")

        genes = list(dict.fromkeys(
            list(result["markers_invasion"].head(n_top).index)
            + list(result["markers_lepidic"].head(n_top).index)
        ))
        pd.Series(genes, name="gene", dtype=str).to_csv(outdir / "marker_heatmap_genes.csv", index=False)
        if genes:
            # Display order only: lepidic first, then invasion; numeric case order.
            # Extract the digits AFTER tkbadeno, not the run ID at the start.
            # This does not change the fitted model or the differential-expression results.
            sort_meta = meta.copy()
            sort_meta["_type_order"] = (
                sort_meta["type"].astype(str).map({"lepidic": 0, "invasion": 1})
            )
            if sort_meta["_type_order"].isna().any():
                raise ValueError("Heatmap types must be lepidic or invasion.")
            # Extract once per pseudobulk row, keeping its original index.
            # sample_id legitimately repeats across histologies. Looking it up
            # with .loc on a duplicated sample_id index multiplies the rows.
            case_number = sort_meta["sample_id"].astype(str).str.extract(
                r"(?i)tkbadeno(\d+)", expand=False
            )
            if case_number.isna().any():
                bad_ids = (
                    sort_meta.loc[case_number.isna(), "sample_id"]
                    .astype(str).drop_duplicates().tolist()
                )
                raise ValueError(f"Cannot extract tkbadeno case numbers: {bad_ids}")
            sort_meta["_case_number"] = case_number.astype(int)
            order = sort_meta.sort_values(
                ["_type_order", "_case_number", "sample_id"]
            ).index
            matrix = log_counts.loc[order, genes].T
            z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=0).replace(0, 1), axis=0)
            z.to_csv(outdir / "marker_heatmap_zscores.csv")
            fig = plt.figure(figsize=(max(8, len(order) * 0.30), max(4, len(genes) * 0.25 + 2)))
            grid = fig.add_gridspec(2, 2, height_ratios=[0.28, max(2, len(genes) * 0.23)], width_ratios=[1, 0.025],
                                    hspace=0.05, wspace=0.04)
            ax_top = fig.add_subplot(grid[0, 0])
            ax = fig.add_subplot(grid[1, 0])
            cax = fig.add_subplot(grid[1, 1])
            group_codes = meta.loc[order, "type"].astype(str).map({"lepidic": 0, "invasion": 1}).to_numpy()[None, :]
            sns.heatmap(group_codes, ax=ax_top, cmap=ListedColormap([COLORS[t] for t in TYPES]), vmin=0, vmax=1,
                        cbar=False, xticklabels=False, yticklabels=False)
            sns.heatmap(z, ax=ax, cmap="RdBu_r", vmin=-2, vmax=2, center=0, cbar_ax=cax,
                        cbar_kws={"label": "Gene-wise z-score", "extend": "both"}, yticklabels=True, xticklabels=True)
            labels = [f"{meta.at[i, 'sample_id']} | {'L' if meta.at[i, 'type'] == 'lepidic' else 'I'}" for i in order]
            ax.set_xticklabels(labels, rotation=90, fontsize=8)
            ax.tick_params(axis="y", labelsize=9)
            ax.set(xlabel="Biological sample | histology", ylabel="Selected genes")
            # Separate histology groups rather than drawing a line after every pair.
            n_lepidic = int(meta.loc[order, "type"].eq("lepidic").sum())
            if 0 < n_lepidic < len(order):
                ax.axvline(n_lepidic, color="white", linewidth=1.5)
                ax_top.axvline(n_lepidic, color="white", linewidth=1.5)
            ax_top.legend(handles=[Patch(color=COLORS[t], label=t) for t in TYPES],
                          frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.1))
            fig.suptitle("Selected genes: log2(normalized counts + 1), standardized within each gene", y=1.015, fontsize=11)
            finish(fig, "marker_heatmap")

            if paired:
                delta = result["paired_log2_differences"].loc[:, genes].T
                delta.to_csv(outdir / "paired_difference_heatmap_values.csv")
                limit = max(0.1, float(np.quantile(np.abs(delta.to_numpy()), 0.98)))
                fig, ax = plt.subplots(figsize=(max(7, delta.shape[1] * 0.38), max(4, len(genes) * 0.25 + 1.7)))
                sns.heatmap(delta, ax=ax, cmap="RdBu_r", vmin=-limit, vmax=limit, center=0,
                            cbar_kws={"label": "Within-sample log2 difference", "extend": "both"},
                            xticklabels=True, yticklabels=True)
                ax.set(xlabel="Biological sample ID", ylabel="Selected genes",
                       title="Paired differences: log2(N invasion + 1) - log2(N lepidic + 1)\nN = size-factor-normalized counts; descriptive, not model LFC")
                ax.tick_params(axis="x", rotation=90)
                fig.tight_layout()
                finish(fig, "paired_difference_heatmap")

            profile_genes = list(result["markers_invasion"].head(3).index) + list(result["markers_lepidic"].head(3).index)
            ncols = min(3, len(profile_genes))
            nrows = int(np.ceil(len(profile_genes) / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), squeeze=False)
            for ax, gene in zip(axes.ravel(), profile_genes):
                values = meta[["sample_id", "type"]].copy()
                values["expression"] = log_counts[gene]
                if paired:
                    wide = values.pivot(index="sample_id", columns="type", values="expression").reindex(columns=TYPES)
                    for _, row in wide.iterrows():
                        ax.plot([0, 1], row.to_numpy(), color="0.75", linewidth=0.8)
                    for x, label in enumerate(TYPES):
                        ax.scatter(np.repeat(x, len(wide)), wide[label], color=COLORS[label], s=23, zorder=2)
                else:
                    sns.boxplot(data=values, x="type", y="expression", order=TYPES, color="white", showfliers=False, ax=ax)
                    sns.stripplot(data=values, x="type", y="expression", order=TYPES, hue="type", palette=COLORS,
                                  jitter=False, legend=False, ax=ax)
                row = res.loc[gene]
                ax.set(xticks=[0, 1], xticklabels=TYPES, xlabel="", ylabel="log2(normalized counts + 1)",
                       title=f"{gene}\nlog2FC={row['log2FoldChange']:.2f}; padj={row['padj']:.2g}")
                sns.despine(ax=ax)
            for ax in axes.ravel()[len(profile_genes):]:
                ax.set_visible(False)
            fig.tight_layout()
            finish(fig, "top_gene_profiles")
        else:
            print("No genes meet FDR and effect-size criteria; marker plots were skipped.")
    # On reruns, remove only known workflow outputs that are no longer applicable.
    conditional_outputs = {
        "pca": ["pca_scores.csv", "pca_genes.csv"],
        "marker_heatmap": ["marker_heatmap_zscores.csv"],
        "paired_difference_heatmap": ["paired_difference_heatmap_values.csv"],
        "top_gene_profiles": [],
    }
    for name, tables in conditional_outputs.items():
        if name not in files:
            for filename in [f"{name}.pdf", f"{name}.tiff", *tables,]:
                (outdir / filename).unlink(missing_ok=True)
    (outdir / "figure_manifest.json").write_text(json.dumps(files, indent=2), encoding="utf-8")
    return files


def run_pseudobulk_de(
    adata, *, sample_key, type_key="type", counts_source="counts", outdir="xenium_pseudobulk_DE",
    paired=True, min_cells=50, min_count=10, min_samples=3, exclude_genes=None,
    n_cpus=4, alpha=0.05, lfc_cutoff=0.5, size_factors_fit_type="poscounts",
    control_genes=None, n_top=20, show=True,
):
    """Convenience wrapper: aggregate, fit, export tables, and make figures."""
    pb = make_pseudobulk(
        adata, sample_key=sample_key, type_key=type_key, counts_source=counts_source,
        paired=paired, min_cells=min_cells, min_count=min_count, min_samples=min_samples,
        exclude_genes=exclude_genes,
    )
    result = fit_pydeseq2(
        pb, outdir=outdir, n_cpus=n_cpus, alpha=alpha, lfc_cutoff=lfc_cutoff,
        size_factors_fit_type=size_factors_fit_type, control_genes=control_genes,
    )
    result["figures"] = plot_results(result, n_top=n_top, show=show)
    return result
