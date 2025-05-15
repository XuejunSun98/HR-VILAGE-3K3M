<p align="right">
  <a href="docs/index.html">Document</a>
</p>

# HR-VILAGE-3K3M: Human Respiratory Viral Immunization Longitudinal Gene Expression

This repository provides the HR-VILAGE-3K3M dataset, a curated collection of human longitudinal gene expression profiles, antibody measurements, and aligned metadata from respiratory viral immunization and infection studies. The dataset includes baseline transcriptomic profiles and covers diverse exposure types (vaccination, inoculation, and mixed exposure). HR-VILAGE-3K3M is designed as a benchmark resource to support the development and evaluation of deep learning models for longitudinal gene expression analysis and to facilitate research into the temporal dynamics of immune responses.

![Flowchart](./Flowchart.png)
**Fig**: Overview of HR-VILAGE-3K3M. (a) HR-VILAGE-3K3M construction workflow. (b) Distribution of sample timepoints for vaccine and inoculation studies, , shown separately for bulk RNA-seq and single-cell RNA-seq datasets. (c) Composition of the dataset, stratified by platform, tissue type, study type, and pathogen, including both bulk and single-cell transcriptomic studies.

# Dataset Description:
- This repo contains 66 studies (59 bulk and 7 single cell studies), comprising 3178 subjects and 14,136 observations along with 2,557,942 single cells.
- We provide preprocessed and normalized gene expression data, raw gene expression data, metadata and antibody data.

# Data Structure:
```
HR-VILAGE-3K3M/
├── README.md
├── study_meta.csv
├── bulk_gene_expr/
│    └── <study_ID>_gene_expr.csv
├── singel_cell_gene_expr/
│    └── <study_ID>.h5ad
├── meta/
│    └── <study_ID>_meta.csv
├── bulk_gene_expr_raw/
│    └── <study_ID>_raw.csv
└── antibody/
     └── <study_ID>_antibody.csv
```
- **study_meta.csv**: Contains study-level metadata (e.g., platform, tissue type, study type) and serves as an overview of the repository. Users can use this file to filter and select a subset of studies based on custom criteria.
- **bulk_gene_expr/**: Processed gene expression matrices for each study (sample-by-gene). All files share the same 41,667 gene columns, with 9,004 genes non-missing across all studies.
- **singel_cell_gene_expr/**: Processed single cell gene expression matrices for each study in h5ad format. Raw count matrices are stored in the .X attribute, and cell-level metadata is stored in the .obs dataframe. 
- **meta/**: Study-specific metadata files (.csv), aligned by row names with the corresponding expression data. All metadata files share the same column structure; missing values are marked as NA.
- **bulk_gene_expr_raw/**: Raw probe-level expression matrices (probe-by-sample), provided when available to support custom preprocessing workflows.
- **antibody/**: Raw antibody measurements with sample IDs matching those in the metadata, enabling integration with gene expression data at the subject level.

# How to start:
To use the HR-VILAGE-3K3M dataset, please access the dataset via [Hugging Face](https://huggingface.co/datasets/xuejun72/HR-VILAGE-3K3M).

Users could directly download all files and read files locally.

Alternatively, the following provides (partially) loading the dataset into Python using `dataset` package.
```python
repo_id = "xuejun72/HR-VILAGE-3K3M"
import pandas as pd
from datasets import load_dataset
```
## Bulk gene expression data

Bulk gene expression data can be loaded and combined using two alternative approaches.

1. Use our predefined configuration name, and pass to `name` argument in `load_dataset()`. `trust_remote_code=True` and `revision="script"` are required.

   Example 1, to download study_meta:
   ```python
   study_meta = load_dataset(repo_id, name="study_meta", trust_remote_code=True, revision="script")["train"].to_pandas()
   ```
   Example 2, to download and combine **all** meta datasets:
   ```python
   meta_dict = load_dataset(repo_id, name = "meta", trust_remote_code=True, revision="script")
   meta = meta_dict["train"].to_pandas()
   ```
   Example 3, to download and combine **all** bulk gene expression datasets. However, this is highly **NOT** recommended since their size are too large and the execution time will be long.
   ```python
   # Not recommended!
   gene_expr_dict = load_dataset(repo_id, name = "bulk_gene_expr", trust_remote_code=True, revision="script")
   ```

   In addition, we provide a study filter function before downloading and loading, which works for **meta** and **bulk gene expression** datasets.
   `split_filter` argument is designed for this filter, which is optional. By default, `split_filter=None` will download all datasets as shown before.
   `split_filter` is a `dict` Python object where `key` is filter factors taking values from `['study_type','platform','tissue','pathogen','vaccine']` and `value` is a list of categories for each `key`. `value` should be exact same as that in study_meta.csv.
   Some examples of a valid `split_filter`:
   ```python
   split_filter = {"study_type": ["inoculation","inoculation"]}
   split_filter = {"study_type": ["vaccine"], "vaccine": ["Influenza TIV"]}
   split_filter = {"study_type": ["vaccine"], "platform": ["RNA-seq"], "tissue": ["PBMC","nasal swab"], "pathogen": []}
   ```
   Example 4, to download and combine a customized filtered meta dataset:
   ```python
   split_filter = {"study_type": ["vaccine"], "platform": ["RNA-seq"], "tissue": ["PBMC","nasal swab"], "pathogen": []}
   meta_filtered_dict = load_dataset(repo_id, name = "meta", trust_remote_code=True, split_filter=split_filter, revision="script")
   for _, value in meta_filtered_dict.items(): 
     meta_filtered = value.to_pandas().set_index("row_name")
   meta_filtered
   ```
   Example 5, to download and combine a customized filtered bulk gene expression dataset:
   ```python
   split_filter = {"study_type": ["vaccine"], "platform": ["RNA-seq"], "tissue": ["PBMC","nasal swab"], "pathogen": []}
   gene_expr_filtered_dict = load_dataset(repo_id, name = "bulk_gene_expr", trust_remote_code=True, split_filter=split_filter, revision="script")
   gene_names = gene_expr_filtered_dict["gene_expr_colnames"][0]["gene_names"]
   all_row_names = []
   all_matrix = []
   for batch in next(iter(gene_expr_filtered_dict.items()))[1]:
     all_row_names.extend(batch["row_names"])
     all_matrix.extend(batch["matrix"])
   gene_expr_filtered = pd.DataFrame(all_matrix, index=all_row_names, columns=gene_names)
   gene_expr_filtered
   ```
2. Use exact path of one csv file, and pass to `data_files` argument in `load_dataset()`.

   Example 1, to download study_meta:
   ```python
   study_meta = load_dataset(repo_id, data_files = "study_meta.csv")["train"].to_pandas()
   ```
   Example 2, to download antibody dataset for GSE194378:
   ```python
   antibody = load_dataset(repo_id, data_files = "antibody/GSE194378_antibody.csv")["train"].to_pandas()
   ```
   Example 3, to download raw gene expression dataset for GSE194378:
   ```python
   raw = load_dataset(repo_id, data_files = "bulk_gene_expr_raw/GSE194378_raw.csv")["train"].to_pandas()
   ```
   Note: for **antibody** and **raw gene expression** datasets, since different study has different columns which cannot be simply combined, loading them must using `data_files` argument and be one-by-one.

## Single cell gene expression data
Single-cell gene expression data can be downloaded and accessed using the anndata package.

  ```python
    import anndata as ad
    # Load the GSE195673 dataset
    GSE195673 = ad.read_h5ad("./GSE195673_processed.h5ad")
    # View cell-level metadata
    GSE195673.obs
  ```
A merged dataset containing all 7 studies is also provided, comprising 2,557,942 cells and 13,589 common genes:
  ```python
    import anndata as ad
    # Load the combined dataset
    combined = ad.read_h5ad("./combined.h5ad")
    # View cell-level metadata
    combined.obs
  ```
The combined object includes detailed cell-level metadata such as sample_id, cell_type, sex, donor_id, time_point_day, dataset, covid_status, age, tissue, and study_type. It also contains dimensionality reductions (X_pca, X_umap) and graph-based neighbor information in obsp for downstream analysis.


# Example code:

Single Cell Visualization
```python
import scanpy as sc
# Visualize UMAP colored by dataset
sc.pl.umap(combined, color='dataset', save='_combined_by_dataset.png', show=True)
```

The codes for data processing and reproducing evaluation results in the paper are in [Documents](https://xuejunsun98.github.io/HR-VILAGE-3K3M/docs/index.html). 

# Data Relations:
The following duplicate studies (right nodes) are not included in HR-VILAGE-3K3M, while their source studies (left leaves) are included in HR-VILAGE-3K3M.
```
GSE73072_H3N2_DEE2 ──┐
                      ├── GSE52428
GSE73072_H1N1_DEE4 ──┘

GSE73072_H3N2_DEE2    ├── GSE30550

GSE73072_HRV_UVA   ──┐
GSE73072_RSV_DEE1  ── ├── GSE17156
GSE73072_H3N2_DEE2 ──┘
```

The studies grouped below originate from the same research group and were typically processed using consistent experimental and technical pipelines. Consequently, these studies are expected to exhibit minimal batch effects compared to those conducted independently. We provide this grouping information to assist users interested in pooling data across studies—combining studies from the same group can help mitigate confounding technical variability and enhance statistical power in downstream analyses such as differential expression. Notably, GSE201533 and GSE201534 represent paired bulk RNA-seq and single-cell RNA-seq data from the same subjects, while GSE246937 is paired with GSE246525, offering a valuable opportunity for multi-resolution analyses.
![Study relation](./relation.png)

# How to cite:
