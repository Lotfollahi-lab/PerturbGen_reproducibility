# %%
import os
import sys

import seaborn as sns
import scanpy as sc
import celloracle as co
import datetime
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from matplotlib import style
import argparse

parser = argparse.ArgumentParser(description="Build celloracle GRN")
parser.add_argument('--data_path', type=str, help='Path to preprocessed anndata object')
parser.add_argument('--atac_grn_path', type=str, default=None, help='Path to atac GRN')
parser.add_argument('--n_cells_downsample', type=int, default=None, help='Number of cells to downsample to (if more than this number of cells are present)')
parser.add_argument('--celltype', nargs='+', type=str, default=None,
                    help='List of cell types to include. If not provided, all cell types are included.')
parser.add_argument('--n_jobs', type=int, default=1, help='Number of parallel jobs to use')
parser.add_argument('--file_name', type=str, help='Base name for output files')
parser.add_argument('--outdir', type=str, help='Output directory')

args = parser.parse_args()


# %%
date = datetime.datetime.now().strftime("%Y%m%d")

# %%
style.use('default')
style.use(
    '/lustre/scratch126/cellgen/lotfollahi/kl11/'
    'T_perturb/perturbgen/pp/mpl_style.mplstyle'
)


# %%
# mk dir

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

# %%
# load preprocessed data
full_data = sc.read_h5ad(args.data_path)

# filter per celltype
if args.celltype is not None:
    full_data = full_data[full_data.obs['celltype_v2'].isin(args.celltype)].copy()

# %%
# # Random downsampling into 30K cells if the anndata object include more than 30 K cells.
if args.n_cells_downsample is not None:
    if full_data.shape[0] > args.n_cells_downsample:
        # Let's dowmsample into 30K cells
        sc.pp.subsample(full_data, n_obs=args.n_cells_downsample, random_state=123)

# %%
print(f"Cell number is :{full_data.shape[0]}")

# %%
base_GRN = co.data.load_human_promoter_base_GRN()


# %%
# Instantiate Oracle object
oracle = co.Oracle()

# %%
# read with .csv if csv file
if args.atac_grn_path is not None:
    Paul_15_data = None
    if args.atac_grn_path.endswith(".csv"):
        tfinfo_df = pd.read_csv(args.atac_grn_path)
    elif args.atac_grn_path.endswith(".parquet"):
        tfinfo_df = pd.read_parquet(args.atac_grn_path)
    else:
        raise ValueError("Unsupported file format for atac_grn_path. Use .csv or .parquet")
else:
    print("No atac_grn_path provided. Using default Paul 2015 TF target data.")
    Paul_15_data = pd.read_csv("/nfs/team361/am74/Cytomeister/notebooks/HSPC_validation_perturbseq/celloracle/TF_data_in_Paul15.csv")
    tfinfo_df = None

# %%
# Check data in anndata
print("Metadata columns :", list(full_data.obs.columns))
print("Dimensional reduction: ", list(full_data.obsm.keys()))

# %%
full_data.X = full_data.layers["raw_count"].copy()

# %%
# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=full_data,
                                   cluster_column_name="celltype_v2",
                                   embedding_name="X_draw_graph_fa")
oracle.import_TF_data(TF_info_matrix=base_GRN)

# %%
if Paul_15_data is not None:
    # Make dictionary: dictionary key is TF and dictionary value is list of target genes.
    TF_to_TG_dictionary = {}

    for TF, TGs in zip(Paul_15_data.TF, Paul_15_data.Target_genes):
        # convert target gene to list
        TG_list = TGs.replace(" ", "").split(",")
        # store target gene list in a dictionary
        TF_to_TG_dictionary[TF] = TG_list

    # We invert the dictionary above using a utility function in celloracle.
    TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)
    TG_to_TF_dictionary = {
        tg.upper(): [tf.upper() for tf in tfs]
        for tg, tfs in TG_to_TF_dictionary.items()
    }

# %%
    oracle.addTFinfo_dictionary(TG_to_TF_dictionary)
else:
    # You can load TF info dataframe with the following code.
    oracle.import_TF_data(TF_info_matrix=tfinfo_df)
# Perform PCA
oracle.perform_PCA()

# %%
# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.savefig(f"{args.outdir}/{args.file_name}_explained_variance_ratio.pdf", dpi=300, bbox_inches='tight')
print(n_comps)
n_comps = min(n_comps, 50)

# %%
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

# %%
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

# %%
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=args.n_jobs)

# %%
# Calculate GRN for each population in "louvain_annot" clustering unit.
# This step may take some time.(~30 minutes)
links = oracle.get_links(cluster_name_for_GRN_unit="celltype_v2", alpha=10,
                         verbose_level=10)
links.links_dict.keys()


# %%
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

# %%
links.plot_degree_distributions(plot_model=True, save=f"{args.outdir}/degree_distribution/")
links.get_network_score()

links.merged_score.head()
links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10,
                              use_cluster_specific_TFdict=True)

# %%
oracle.to_hdf5(f'{args.outdir}/{args.file_name}.celloracle.oracle')
links.to_hdf5(file_path=f'{args.outdir}/{args.file_name}.celloracle.links')
oracle.adata.write_h5ad(f'{args.outdir}/{args.file_name}.h5ad')
