# Xenium_LUAD

`Xenium_LUAD` is the companion code repository for the manuscript:

"Elastic Fiber and Cytokeratin 7 Staining for Deep Learning–Based Assessment of Lung Adenocarcinoma Invasion: Insights from Spatial Transcriptomics"

This repository contains the Xenium spatial-transcriptomics workflow used in the study, from Proseg-based cell segmentation through global clustering, epithelial-cell subclustering, histology-guided cluster annotation, and paired pseudobulk differential-expression analysis. The deep-learning image-classification and whole-slide inference code is outside the scope of this repository.

> [!IMPORTANT]
> This is not a fully automated annotation pipeline. At multiple stages, clustering results must be imported into Xenium Explorer and reviewed together with the corresponding H&E image. Cluster identifiers and histologic labels in the notebooks are specific to the study dataset and must not be transferred directly to another dataset.

## Workflow overview

Run the files in the following order:

```text
Xenium output bundle
        |
        v
1. proseg.py
   Sopa/Proseg segmentation and creation of an Explorer-compatible output bundle
        |
        v
2. global_clustering.ipynb
   Cross-section integration, Harmony correction, global Leiden clustering, and CSV export
        |
        +--> Review clusters in Xenium Explorer with H&E
        v
3. epithelium_clustering_step1.ipynb
   First epithelial-cell selection and reclustering without Harmony correction
        |
        +--> Review clusters in Xenium Explorer with H&E
        v
4. epithelium_clustering_step2.ipynb
   Second epithelial-cell selection, additional filtering, and final reclustering
        |
        +--> Review clusters in Xenium Explorer with H&E
        v
5. annotation_pseudobulk.ipynb
   Histology-guided cluster annotation and paired lepidic-versus-invasion pseudobulk DE analysis
        |
        v
   Tables and publication-quality figures
```

## Repository contents

| File | Purpose |
|---|---|
| `proseg.py` | Reads one Xenium output bundle with Sopa, runs Proseg segmentation, aggregates transcript counts, stores raw counts in `layers["counts"]`, and writes an Explorer-compatible directory containing `adata.h5ad`. |
| `global_clustering.ipynb` | Loads all Proseg-derived `adata.h5ad` files, filters cells, normalizes expression, performs PCA, Harmony integration, neighbor-graph construction, Leiden clustering, UMAP, marker visualization, and per-section CSV export. |
| `epithelium_clustering_step1.ipynb` | Selects epithelial-enriched global clusters, restores raw counts, renormalizes the selected cells, reclusters them, and exports results for spatial review. |
| `epithelium_clustering_step2.ipynb` | Selects the epithelial clusters retained after the first review, applies the final epithelial filtering and clustering settings, and writes `epithelium.h5ad`. |
| `annotation_pseudobulk.ipynb` | Assigns study-specific histologic labels to final clusters, removes non-tumor sections, constructs paired pseudobulks, fits PyDESeq2, and creates output figures. |
| `xenium_pseudobulk_de.py` | Functions for pseudobulk aggregation, PyDESeq2 fitting, quality-control exports, and visualization. |
| `LICENSE` | GNU General Public License v3.0 only. |

## Software requirements

The segmentation step was run in an environment separate from the downstream clustering and differential-expression analyses. The following core versions were recorded for this workflow:

| Software | Version |
|---|---|
| Python | 3.12 |
| Scanpy | 1.12.1 |
| rapids-singlecell | 0.15.0rc5 |
| Sopa | 2.2.1 |
| Proseg | 3.1.0 |
| PyDESeq2 | 0.5.4 |

### Segmentation environment

- Python 3.12
- Sopa 2.2.1
- Proseg 3.1.0
- Scanpy
- SpatialData and SpatialData-IO

Proseg 3.1.0 was called through Sopa 2.2.1 in `proseg.py`, using an environment separate from the downstream analyses.

### Clustering and differential-expression environment

- Python 3.12
- Scanpy 1.12.1
- rapids-singlecell 0.15.0rc5
- PyDESeq2 0.5.4
- AnnData
- NumPy
- pandas
- SciPy
- Matplotlib
- seaborn
- scikit-learn
- CuPy and RMM compatible with the installed CUDA/RAPIDS stack
- JupyterLab or Jupyter Notebook

The clustering notebooks use GPU-backed CuPy, RMM, and rapids-singlecell operations and therefore require a compatible NVIDIA CUDA environment. Installation of the RAPIDS components is platform- and CUDA-version-specific.

## Input data

`proseg.py` expects a single standard Xenium output bundle as `input_path`. Run the script once for each section. Each run should write an Explorer-compatible output directory such as:

```text
<proseg_outputs>/<sample_id>.explorer/
└── adata.h5ad
```

The downstream notebooks expect the resulting AnnData objects to contain at least:

- raw, nonnegative integer transcript counts in `adata.layers["counts"]`;
- a `cell_id` column in `adata.obs` that can be matched to the corresponding Xenium Explorer bundle; and
- a sample identifier that is added as `adata.obs["batch"]` during concatenation.

Raw Xenium data and tissue images are not distributed with this repository. Refer to the Data Availability statement in the accompanying manuscript.

## Directory configuration

Create the output directories before running the notebooks. One possible layout is:

```text
analysis_root/
├── proseg_outputs/
├── global_clustering_outputs/
├── epithelial_clustering_step1_outputs/
├── epithelial_clustering_step2_outputs/
└── annotation_outputs/
```

For example:

```bash
mkdir -p analysis_root/{proseg_outputs,global_clustering_outputs,epithelial_clustering_step1_outputs,epithelial_clustering_step2_outputs,annotation_outputs}
```

Edit the path variables near the beginning of each file:

```python
input_path = "path/to/the/input_xenium_bundle"
proseg_outputs = "path/to/your/proseg/output/directory"
global_clustering_outputs = "path/to/your/global/clustering/output/directory"
epithelial_clustering_step1_outputs = "path/to/your/epithelial/step1/output/directory"
epithelial_clustering_step2_outputs = "path/to/your/epithelial/step2/output/directory"
annotation_outputs = "path/to/your/annotation_and_pseudobulk/output/directory"
```

In `proseg.py`, `explorer_path` must be the complete output directory for one section, including the `.explorer` suffix, rather than only the parent directory.

The `section_ids` lists in the clustering notebooks must also be replaced with the Xenium bundle directory names used in the analysis. The uploaded study notebooks use a dataset-specific naming convention; using an explicit mapping between each bundle directory and its sample ID is recommended.

## Step-by-step use

### 1. Proseg segmentation

For each Xenium section, edit `input_path` and `explorer_path` in `proseg.py`, then run:

```bash
python proseg.py
```

The script:

1. reads the Xenium bundle with Sopa;
2. creates transcript patches;
3. runs Proseg with 16 threads;
4. aggregates transcripts into the resulting cell boundaries;
5. stores the aggregated raw counts in `layers["counts"]`; and
6. writes an Explorer-compatible output with `save_h5ad=True`.

For additional details on Sopa, see:

> Blampey Q, Mulder K, Gardet M, et al. Sopa: a technology-invariant pipeline for analyses of image-based spatial omics. *Nature Communications*. 2024;15:4981. [https://doi.org/10.1038/s41467-024-48981-z](https://doi.org/10.1038/s41467-024-48981-z)

### 2. Global clustering

Open and execute `global_clustering.ipynb`.

The current notebook applies the following principal settings:

| Operation | Setting |
|---|---|
| Cell filtering | `min_counts=100`; `min_genes=30` |
| Library-size normalization | `target_sum=500` followed by `log1p` |
| PCA | All retained panel genes; `random_state=7` |
| Batch correction | Harmony using `batch` |
| Neighbors | `n_neighbors=15`; representation `X_pca_harmony` |
| Leiden clustering | `resolution=1.0`; `random_state=7` |
| UMAP | `min_dist=0.5`; `spread=1.0`; `random_state=7` |

The principal integrated output is:

```text
global_clustering_outputs/proseg_harmony_integrated.h5ad
```

For every section, the notebook also exports a CSV containing:

```text
cell_id, group, color
```

Import the corresponding CSV into Xenium Explorer and inspect the spatial distribution of each Leiden cluster together with the H&E image. Update the `keep` set in `epithelium_clustering_step1.ipynb` to include only epithelial-enriched clusters.

### 3. First epithelial clustering

Open and execute `epithelium_clustering_step1.ipynb` after updating its `keep` set.

The notebook restores `layers["counts"]` to `X`, normalizes to 500 transcripts per cell, applies `log1p`, and performs clustering without Harmony correction.

| Operation | Setting |
|---|---|
| Neighbors | `n_neighbors=15`; representation `X_pca` |
| Leiden clustering | `resolution=0.8`; `random_state=7` |
| UMAP | `min_dist=0.5`; `spread=1.0`; `random_state=7` |

The principal output is:

```text
epithelial_clustering_step1_outputs/proseg_nocorrection_ep.h5ad
```

Import the per-section CSV files into Xenium Explorer, review the clusters with H&E, and update the `keep` set in `epithelium_clustering_step2.ipynb`.

### 4. Second epithelial clustering

Open and execute `epithelium_clustering_step2.ipynb` after updating its `keep` set.

The uploaded notebook currently uses the following settings:

| Operation | Setting |
|---|---|
| Additional epithelial-cell filtering | `min_counts=200`; `min_genes=80` |
| Library-size normalization | `target_sum=500` followed by `log1p` |
| PCA representation used downstream | First 6 PCs (`X_pca_6`) |
| Neighbors | `n_neighbors=5` |
| Leiden clustering | `resolution=0.45`; `random_state=7` |
| UMAP | `min_dist=0.5`; `spread=1.0`; `random_state=7` |

The principal output is:

```text
epithelial_clustering_step2_outputs/epithelium.h5ad
```

Review the final per-section clustering CSV files in Xenium Explorer. Then edit the cluster-to-histology assignments in `annotation_pseudobulk.ipynb` on the basis of the spatial and histologic review.

### 5. Histologic annotation and paired pseudobulk DE analysis

Open `annotation_pseudobulk.ipynb` and update:

- the final `keep` set;
- the mapping from Leiden cluster to `type`;
- the set of non-tumor sample IDs; and
- the output and analysis parameters.

The study-specific labels include `lepidic`, `invasion`, `both`, `pneumocyte`, and `bronchial_cells`. Only cells labeled exactly `lepidic` or `invasion` are included in the DE comparison.

The uploaded notebook uses:

| Parameter | Value |
|---|---:|
| Biological-sample column | `batch` |
| Histology column | `type` |
| Raw-count source | `layers["counts"]` |
| Paired design | `True` |
| Minimum cells per sample × histology pseudobulk | 100 |
| Gene count threshold | At least 30 counts |
| Minimum pseudobulks meeting the gene threshold | 3 |
| FDR threshold | 0.01 |
| Descriptive absolute log2 fold-change cutoff | 1.0 |
| Size-factor fit type | `poscounts` |
| Top genes displayed per direction | 20 |

For paired analysis, PyDESeq2 is fitted with:

```text
~ sample_id + type
```

and the contrast is:

```text
invasion versus lepidic
```

Therefore, a positive `log2FoldChange` indicates higher expression in the invasive component. The log2 fold-change cutoff is applied after the Wald test and FDR correction; it is not a formal test of the composite null hypothesis that the absolute effect exceeds the cutoff.

`sample_id` must represent an independent biological unit, usually one patient with one tumor. Repeated regions or sections from the same biological unit must not be treated as independent replicates.

## Main pseudobulk outputs

The DE workflow writes quality-control tables, model results, normalized expression matrices, configuration metadata, and figures, including:

```text
sample_qc.csv
gene_qc.csv
pseudobulk_counts.csv
pseudobulk_metadata.csv
deseq2_all_genes.csv
deg_fdr.csv
markers_invasion.csv
markers_lepidic.csv
normalized_counts.csv
log2_normalized_counts.csv
paired_log2_differences.csv
analysis_config.json
figure_manifest.json
cells_per_pseudobulk.pdf / .tiff
raw_counts_per_pseudobulk.pdf / .tiff
pca.pdf / .tiff
marker_heatmap.pdf / .tiff
paired_difference_heatmap.pdf / .tiff
top_gene_profiles.pdf / .tiff
```

Some figures are generated only when the required number of genes or paired samples is available.

## Study-specific elements that must be reviewed

Before applying the workflow to another dataset, review all of the following:

- the `section_ids` lists and sample-ID parsing;
- the global and epithelial `keep` sets;
- the final cluster-to-histology mapping;
- the non-tumor sample list;
- the marker-gene dictionaries;
- the number of PCs, neighbor settings, and Leiden resolutions;
- the cell- and gene-level QC thresholds; and
- whether `batch` truly identifies an independent biological sample for pseudobulk analysis.

The pseudobulk plotting code also contains study-specific natural sorting for sample names containing `tkbadeno<number>`.

## Reproducibility notes

- Random seeds are set to 7 for PCA, neighbor construction, Leiden clustering, and UMAP where supported.
- The DE workflow records the analysis settings and package versions in `analysis_config.json`.
- Raw counts, rather than normalized or log-transformed values, are required for PyDESeq2.
- The notebooks currently contain manually selected cluster identifiers. These identifiers can change if preprocessing, package versions, or clustering parameters change.

## License

The source code in this repository is licensed under the **GNU General Public License v3.0 only** (`GPL-3.0-only`). See `LICENSE` for the full license text.

The license applies to the source code in this repository. Xenium data, histology images, trained models, manuscript text, and other study materials are not included unless explicitly stated and may be subject to separate terms.

## Disclaimer

This repository is provided for research and reproducibility purposes only. It has not been validated or approved for clinical diagnosis, treatment selection, or other patient-care use. The software is provided without warranty, as described in the GPL-3.0 license.

## Citation

Please cite the accompanying manuscript when using this workflow. The bibliographic information will be updated after publication.
