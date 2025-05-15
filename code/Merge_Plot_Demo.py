"""
Preprocessing and Merging Single-Cell Datasets
Author: Yiran Song
Date: May 10, 2025

This script preprocesses and integrates all .h5ad single-cell RNA-seq datasets. 

The pipeline includes:
1. Identifying common genes across all datasets
2. Harmonizing metadata and ensuring consistent formats
3. Ensuring .X matrix exists and is in sparse format
4. Concatenating all datasets into a single AnnData object
5. Saving UMAP plots colored by dataset and cell type

Outputs:
- combined_filtered_concat.h5ad: merged AnnData object
- merged_umap_plots.pdf: UMAP visualizations
- common_genes.txt: list of intersected gene names
"""

import os
import gc
import scanpy as sc
import anndata as ad
import pandas as pd
import scipy.sparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# --- Configuration ---
input_dir = "./meta_processed"
output_dir = input_dir
raw_data_paths = {
    fname.replace("_processed.h5ad", ""): os.path.join(input_dir, fname)
    for fname in os.listdir(input_dir) if fname.endswith("_processed.h5ad")
}

# --- Helper functions ---
def fix_var_names(adata):
    adata.var_names = pd.Index(adata.var_names.astype(str))
    adata.var_names_make_unique()
    return adata

def ensure_sparse_X(adata, sample_name=None):
    if adata.X is None:
        print(f" .X is None in {sample_name or 'unknown sample'}")
        if 'counts' in adata.layers:
            print("→ Using adata.layers['counts']")
            adata.X = adata.layers['counts']
        elif 'raw_counts' in adata.layers:
            print("→ Using adata.layers['raw_counts']")
            adata.X = adata.layers['raw_counts']
        elif 'log_norm' in adata.layers:
            print("→ Using adata.layers['log_norm']")
            adata.X = adata.layers['log_norm']
        elif adata.raw is not None:
            print("→ Using adata.raw.X")
            adata.X = adata.raw.X
        else:
            raise ValueError(f"Cannot find source for .X in {sample_name or 'unknown sample'}")
    
    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    return adata

# --- Step 1: Determine common genes across all datasets ---
gene_sets = []
for sample_name, file_path in raw_data_paths.items():
    adata = ad.read_h5ad(file_path)
    gene_sets.append(set(adata.var_names))
    del adata
    gc.collect()

common_genes = set.intersection(*gene_sets)
print(f"Number of common genes: {len(common_genes)}")

with open(os.path.join(output_dir, "common_genes.txt"), "w") as f:
    for gene in sorted(common_genes):
        f.write(f"{gene}\n")
print(f"Saved {len(common_genes)} common genes to file.")

# --- Step 2: Harmonize and prepare each AnnData object ---
common_fields = [
    "donor_id", "sample_id", "time_point_day", "covid_status",
    "sex", "age", "tissue", "cell_type", "dataset"
]
adata_list = []

for sample_name, file_path in raw_data_paths.items():
    adata = ad.read_h5ad(file_path)

    # Subset to common genes
    filtered = adata[:, list(common_genes)].copy()

    # Ensure .X exists and is sparse
    filtered = ensure_sparse_X(filtered, sample_name)

    # Make obs names unique
    filtered.obs_names = [f"{idx}_{sample_name}" for idx in filtered.obs_names]

    # Ensure all common_fields exist and are categorical or strings
    for col in common_fields:
        if col not in filtered.obs.columns:
            filtered.obs[col] = pd.Series([pd.NA] * filtered.shape[0], dtype="category")
        elif pd.api.types.is_numeric_dtype(filtered.obs[col]):
            filtered.obs[col] = filtered.obs[col].astype(str)
        else:
            filtered.obs[col] = filtered.obs[col].astype("category")

    # Remove duplicated columns in .obs
    filtered.obs = filtered.obs.loc[:, ~filtered.obs.columns.duplicated()]

    # Convert all obs columns to string to avoid dtype mismatch
    filtered.obs = filtered.obs.astype(str)

    adata_list.append(filtered)
    print(f"Loaded and filtered {sample_name}")

# --- Final check: all have .X defined ---
for i, adata in enumerate(adata_list):
    if adata.X is None:
        print(f"Missing .X in dataset #{i} with obs_names like {adata.obs_names[:3].tolist()}")
        print(f"Available layers: {adata.layers.keys()}")
        print(f"Has .raw: {adata.raw is not None}")
        raise ValueError("One or more AnnData objects still lack .X.")

# --- Step 3: Concatenate all ---
adata_all = ad.concat(
    adata_list,
    join='inner',
    label='dataset',
    index_unique=None
)

# --- Save result ---
final_path = os.path.join(output_dir, "combined_filtered_concat.h5ad")
adata_all.write(final_path)
print(f"Saved merged object: {final_path}")
print(f"Final shape: {adata_all.shape}")

# --- Load data ---
# adata = sc.read_h5ad("combined_filtered_concat_processed.h5ad")
# print(f"Loaded AnnData: {adata.shape}")

# --- Plotting to PDF ---
with PdfPages("merged_umap_plots.pdf") as pdf:
    print("📊 Plotting UMAP colored by 'dataset'...")
    sc.pl.umap(adata, color="dataset", show=False)
    pdf.savefig()  # Save current figure
    plt.close()

    if "cell_type" in adata.obs.columns:
        print("📊 Plotting UMAP colored by 'cell_type'...")
        sc.pl.umap(adata, color="cell_type", show=False)
        pdf.savefig()
        plt.close()

    print("UMAP plots saved to merged_umap_plots.pdf")
