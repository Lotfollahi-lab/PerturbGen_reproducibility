# %%
import scanpy as sc
import anndata as ad
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
import os
# from perturbgen.configs import ROOT
from pathlib import Path
import datetime
import random
from matplotlib import style
# import f1 score and accuracy
from sklearn.metrics import (
    f1_score, 
    accuracy_score, 
    precision_score, 
    recall_score,
    matthews_corrcoef,
    balanced_accuracy_score
)
import re

parser = argparse.ArgumentParser()

parser.add_argument('--model', type=str, default='perturbgen', help='Model name to include in output filenames')
parser.add_argument('--oracle_path', type=str, default=None, help='Path to CellOracle model directory')
parser.add_argument('--h5ad_files', type=str, default=None, help='Path to directory containing h5ad files')
parser.add_argument('--tf_excl_celltype', type=str, nargs='+', default=['Early GMP', 'MEP', 'Megakaryocyte Precursor', 'EoBasoMast Precursor'], help='Cell types to exclude from analysis')
parser.add_argument('--hspc_incl_celltype', type=str, nargs='+', default=['Early_Ery', 'Late_Ery'], help='Cell types to include from analysis')
parser.add_argument('--logfc_threshold', type=float, default=0.25, help='Log fold change threshold for DEGs')
parser.add_argument('--confusion_fig_name', type=str, default='confusion_matrix_stem_topN_degs_grouped', help='Name of the output confusion matrix figure')
parser.add_argument('--res_dir', type=str, default='T_perturb/res/hspc/tf_perturbation/', help='Directory to save results')
parser.add_argument('--corr_fig_name', type=str, default='mse_erythro_stem_logfc_correlation', help='Name of the output figure')
parser.add_argument('--table_name', type=str, default='tf_stem_perturbation_evaluation', help='Name of the output table')
parser.add_argument('--direction_fig_name', type=str, default='precision_f1_stem_topN_degs_grouped', help='Name of the output figure for each gene')
parser.add_argument('--n_prop', type=str, default='precision_f1_stem_topN_degs_grouped', help='Name of the output figure for each gene')
args = parser.parse_args()
# %%
# os.chdir(ROOT)
print(f'Current working directory: {os.getcwd()}')

if args.model == 'celloracle':
    import celloracle as co
    print(f'Using CellOracle version: {co.__version__}')

# %%
date = datetime.datetime.now().strftime("%Y%m%d")

# %%
np.random.seed(42)
random.seed(42)

# %%
# os.chdir(ROOT)
# print(f'Current working directory: {os.getcwd()}')

# %%
style.use('default')
style.use(
    '/lustre/scratch126/cellgen/lotfollahi/kl11/'
    'T_perturb/perturbgen/pp/mpl_style.mplstyle'
)

# %%
res_dir = args.res_dir
# create directory if it does not exist
if not os.path.exists(res_dir):
    os.makedirs(res_dir)

# %%
tf_data = sc.read_h5ad('/nfs/team361/am74/Cytomeister/Evaluation_datasets/HSPC_invitro_perturbseq/processed_data.h5ad')

# %%
if args.model == 'perturbgen':
    if args.h5ad_files is None:
        raise ValueError("Please provide --h5ad_files path for perturbgen model")
    else:
        h5ad_path = Path(args.h5ad_files)
        files_h5ad = sorted([str(p) for p in h5ad_path.glob("*.h5ad")])
        print(f"Found {len(files_h5ad)} files ending with .h5ad")
        # read first annData object to subset tf_data to common genes
        adata = sc.read_h5ad(files_h5ad[0])
        tf_data = tf_data[:, tf_data.var_names.isin(adata.var_names)]
        del adata
elif args.model == 'celloracle':
    oracle_data = co.load_hdf5(args.oracle_path)
    tf_data = tf_data[:, tf_data.var_names.isin(oracle_data.adata.var_names)]
# %%
tf_filtered = tf_data[~tf_data.obs['new_CellType'].isin(args.tf_excl_celltype)].copy()
# # log normalize on subsetted data 
tf_filtered.X = tf_filtered.layers['counts'].copy()
sc.pp.normalize_total(tf_filtered)
sc.pp.log1p(tf_filtered)
# %%
sc.tl.rank_genes_groups(
    tf_filtered,
    groupby="target",
    reference='NT',
    method="wilcoxon"
)
tf_deg_df = sc.get.rank_genes_groups_df(tf_filtered, group=None)
degs_filtered_df = tf_deg_df[(tf_deg_df['pvals_adj'] < 0.05) & tf_deg_df['logfoldchanges'].abs() > 0.25]

# %%
# plot number of DEGs per target
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
degs_count = degs_filtered_df['group'].value_counts().sort_index()
# sort values by number of DEGs
degs_count = degs_count.sort_values(ascending=True)
degs_count.plot(kind='barh')
plt.xlabel('Number of DEGs')
plt.ylabel('Target')
# Add dotted vertical line at 50
plt.axvline(x=50, linestyle='--', color='grey')
# plt.title('Number of DEGs per Target')
plt.tight_layout()
plt.savefig(f'{res_dir}/num_degs_per_target.pdf', dpi=300, bbox_inches='tight')

# %%
# filter for groups with at least 50 DEGs
degs_count_filtered = degs_count[degs_count >= 50]

# random baseline to beat
def generate_random_baseline(n_draws=1000, draw_length=40, seed=42):
    """Generate random direction draws for a baseline.

    Returns:
        rng: numpy Generator
        random_directions: single draw (length draw_length)
        random_draws: array of shape (n_draws, draw_length)
    """
    rng = np.random.default_rng(seed)
    # random_directions = rng.choice([-1, 1], size=draw_length)
    random_draws = rng.choice([-1, 1], size=(n_draws, draw_length))
    return rng, random_draws

def read_perturbed_h5ad_files(
    gene_name,
    tf_filtered,
): 
    '''
    Read perturbed h5ad files for a given gene.

    Returns:
        adata: AnnData object
    '''
    pattern = rf"_g{re.escape(gene_name)}_(?:.*)\.h5ad$"
    regex = re.compile(pattern, flags=re.IGNORECASE)
    match_files = [f for f in files_h5ad if regex.search(os.path.basename(f))]
    
    if len(match_files) == 0:
        print(f"No matching file found for gene: {gene_name}")
        return None
    elif len(match_files) > 1:
        print(f"Warning: Multiple matching files found for gene: {gene_name}. Using the first one.")
        adata = sc.read_h5ad(match_files[0])
        adata = adata[:, adata.var_names.isin(tf_filtered.var_names)]   
        return adata
    else:
        print(f"Found matching file: {match_files[0]}")
        adata = sc.read_h5ad(match_files[0])
        adata = adata[:, adata.var_names.isin(tf_filtered.var_names)]
        return adata
    
def compute_predicted_DEG_response(
        adata
):
    '''
    Compute DEGs between perturbed and predicted cells in adata.

    Returns:
        deg_df: DataFrame with DEGs
    '''
    pred = adata.copy()
    pred.X = pred.layers['pred_counts'].copy()
    pred.obs['status'] = 'predicted'
    adata.obs['status'] = 'perturbed'
    adata = adata.concatenate(pred)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    adata = adata[adata.obs['celltype_v2'].isin(args.hspc_incl_celltype)] 
    
    group_ad = "perturbed"
    ref_ad = "predicted"

    sc.tl.rank_genes_groups(
        adata,
        groupby="status",
        groups=[group_ad],
        reference=ref_ad,
        method="wilcoxon"
    )
    deg_df = sc.get.rank_genes_groups_df(adata, group=group_ad)
    return deg_df

def perturb_celloracle(
        oracle_data,
        gene_to_perturb,
        tf_filtered=tf_filtered,
):
    '''
    induce in silico gene perturbation in CellOracle model

    '''
    # Option A: from links (common)
    available_regs = set(oracle_data.TFdict.keys())

    if gene_to_perturb not in available_regs:
        print(f"Skip {gene_to_perturb}: not in base GRN regulators")
        return None
    else:
        try:
            oracle_data.simulate_shift(
                perturb_condition={gene_to_perturb: 0.0},
                n_propagation=int(args.n_prop)
            )
        except ValueError as e:
            msg = str(e)
            if "not included in the base GRN" in msg:
                print(f"Skip {gene_to_perturb}: {msg}")
                return None
            raise  # re-raise any other unexpected ValueError
        memp_data = oracle_data.adata[oracle_data.adata.obs['celltype_v2'].isin(args.hspc_incl_celltype)].copy()
        memp_data = memp_data[:, memp_data.var_names.isin(tf_filtered.var_names)].copy()
        
        # stack imputed_count and simulated_count, copy obs and add single status column
        orig_obs = memp_data.obs.copy().reset_index(drop=True)
        pert_obs = memp_data.obs.copy().reset_index(drop=True)

        orig_obs['status'] = 'imputed'
        pert_obs['status'] = 'perturbed'

        orig_obs.index = [f'orig_{i}' for i in range(memp_data.n_obs)]
        pert_obs.index = [f'pert_{i}' for i in range(memp_data.n_obs)]

        stacked_obs = pd.concat([orig_obs, pert_obs])
        stacked_obs['status'] = stacked_obs['status'].astype('category')

        perturbed_adata = ad.AnnData(
            X = np.vstack([memp_data.layers['imputed_count'], memp_data.layers['simulated_count']]),
            obs = stacked_obs,
            var = memp_data.var.copy()
        )
        # perform DEG analysis
        sc.tl.rank_genes_groups(
            perturbed_adata,
            group='perturbed',
            groupby='status',
            reference='imputed',
            method='wilcoxon',
        )
        deg_results = sc.get.rank_genes_groups_df(perturbed_adata, group='perturbed')
    return deg_results

def compute_correlation_metrics(
    tf_deg_df,
    pred_deg_df,
    direction_dict,
    alpha,
    logfc_threshold
):
    '''
    Compute correlation metrics between tf_filtered_deg and pred_deg_df.

    Returns:
        corr_metrics: DataFrame with correlation metrics
    '''

    tf_deg_df_filtered = tf_deg_df[(tf_deg_df["pvals_adj"] < alpha)]
    merged = pd.merge(
        pred_deg_df[['names', 'logfoldchanges']].rename(columns={'logfoldchanges': 'l2fc_insilico'}),
        tf_deg_df_filtered[['names', 'logfoldchanges']].rename(columns={'logfoldchanges': 'l2fc_invitro'}),
        on='names'
    )
    merged_filtered = merged[(merged['l2fc_invitro'].abs() > logfc_threshold)]
    # mask = np.isfinite(merged_filtered['l2fc_insilico']) & np.isfinite(merged_filtered['l2fc_invitro'])
    # merged_filtered = merged_filtered[mask]
    x = merged_filtered['l2fc_invitro']
    y = merged_filtered['l2fc_insilico']
    if len(x) >= 2:
        spearman_val, _ = spearmanr(x, y)
        direction_dict['spearman'].append(spearman_val)
        pearson_val, _  = pearsonr(x, y)
        direction_dict['pearson'].append(pearson_val)
        direction_dict['n'].append(len(x))
    else:
        spearman_val, pearson_val = np.nan, np.nan
        direction_dict['spearman'].append(np.nan)
        direction_dict['pearson'].append(np.nan)
        direction_dict['n'].append(np.nan)
    
    return direction_dict, merged_filtered, merged

def plot_correlation(
    merged_filtered,
    direction_dict,
    gene,
    fig_name
):
    '''
    Plot correlation between in vitro and in silico log fold changes.

    '''
    x = merged_filtered['l2fc_invitro']
    y = merged_filtered['l2fc_insilico']
    spearman_val = direction_dict['spearman'][-1]
    pearson_val = direction_dict['pearson'][-1]
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, alpha=0.7)
    # cbar = plt.colorbar(ax.collections[0], ax=ax)
    # cbar.set_label('mse')
    
    ax.axhline(0, color="grey", linestyle="--", linewidth=1)
    ax.axvline(0, color="grey", linestyle="--", linewidth=1)
    ax.set_xlabel(f"log2fc_in_vitro")
    ax.set_ylabel(f"log2fc_in_silico")
    ax.set_title(f"n={len(x)} | Spearman={spearman_val:.2f}, Pearson={pearson_val:.2f}")
    plt.tight_layout()
    plt.savefig(os.path.join(res_dir, f"{fig_name}.pdf"), bbox_inches='tight', dpi=300)
    plt.close()

def deg_direction_metrics(
    df,
    true_col="l2fc_invitro",
    pred_col="l2fc_insilico",
    handle_zeros="as_negative",  # {"as_negative", "drop", "separate"}
):
    """
    Directionality metrics for DEG sign prediction.

    Assumes:
      y_true = sign(true_col) is in {-1, +1} (e.g., after filtering abs(true_col)>thresh)
      y_pred = sign(pred_col) is in {-1, 0, +1}

    handle_zeros:
      - "as_negative": treat pred==0 as incorrect sign (counts as not-up AND not-down)
                       (binary metrics treat it as "not-up" for Up-positive computations)
      - "drop": evaluate only where pred!=0 (conditional-on-calling)
      - "separate": same as "drop" for computed metrics, but still returns coverage stats
    """

    y_true = np.sign(df[true_col].to_numpy())
    y_pred = np.sign(df[pred_col].to_numpy())

    # Optional safety: if true zeros exist, drop them (usually you filtered invitro abs>thresh)
    mask_true = (y_true != 0)
    y_true = y_true[mask_true]
    y_pred = y_pred[mask_true]

    if handle_zeros in {"drop", "separate"}:
        mask = (y_pred != 0)
        y_true_f = y_true[mask]
        y_pred_f = y_pred[mask]
    elif handle_zeros == "as_negative":
        y_true_f = y_true
        y_pred_f = y_pred
    else:
        raise ValueError("handle_zeros must be one of {'as_negative','drop','separate'}")

    # Confusion components for Up = +1 as "positive"
    tp = int(np.sum((y_true_f ==  1) & (y_pred_f ==  1)))
    fp = int(np.sum((y_true_f == -1) & (y_pred_f ==  1)))
    fn = int(np.sum((y_true_f ==  1) & (y_pred_f !=  1)))   # -1 or 0 (if present in eval set)
    tn = int(np.sum((y_true_f == -1) & (y_pred_f !=  1)))   # -1 or 0

    # Basic rates for Up (+1)
    precision_up = tp / (tp + fp) if (tp + fp) else 0.0
    recall_up    = tp / (tp + fn) if (tp + fn) else 0.0

    # Directional accuracy = exact sign match in evaluated set
    accuracy = float(np.mean(y_true_f == y_pred_f)) if len(y_true_f) else np.nan

    # Binary labels for sklearn (Up as positive)
    yt_up = (y_true_f == 1).astype(int)
    yp_up = (y_pred_f == 1).astype(int)

    # Up-positive metrics (these equal precision_up/recall_up, but computed robustly)
    precision_up_skl = float(precision_score(yt_up, yp_up, zero_division=0))
    recall_up_skl    = float(recall_score(yt_up, yp_up, zero_division=0))
    f1_up            = float(f1_score(yt_up, yp_up, zero_division=0))

    # Down-positive metrics (treat -1 as positive)
    yt_down = (y_true_f == -1).astype(int)
    yp_down = (y_pred_f == -1).astype(int)
    precision_down = float(precision_score(yt_down, yp_down, zero_division=0))
    recall_down    = float(recall_score(yt_down, yp_down, zero_division=0))
    f1_down        = float(f1_score(yt_down, yp_down, zero_division=0))

    # Balanced accuracy / MCC for the Up-vs-notUp framing
    bal_acc = float(balanced_accuracy_score(yt_up, yp_up)) if len(y_true_f) else np.nan
    mcc = float(matthews_corrcoef(yt_up, yp_up)) if len(np.unique(yt_up)) > 1 or len(np.unique(yp_up)) > 1 else 0.0

    # Macro/Weighted F1 across the two direction classes (Up vs Down)
    # (works only if y_pred_f has no zeros; if zeros exist, they're neither class)
    # For "as_negative", zeros can exist; for "drop"/"separate" they don't.
    if np.all(np.isin(y_pred_f, [-1, 1])) and len(y_true_f):
        y_true_2 = (y_true_f == 1).astype(int)  # 1=Up, 0=Down
        y_pred_2 = (y_pred_f == 1).astype(int)
        f1_macro = float(f1_score(y_true_2, y_pred_2, average="macro", zero_division=0))
        f1_weighted = float(f1_score(y_true_2, y_pred_2, average="weighted", zero_division=0))
    else:
        # still meaningful to return macro from per-class F1s (Up/Down) even if zeros existed,
        # but note: Down F1 computed with yp_down counts zeros as "not-down" (i.e. FN for Down)
        f1_macro = float((f1_up + f1_down) / 2.0) if len(y_true_f) else np.nan
        # weighted by true class support
        sup_up = int(np.sum(y_true_f == 1))
        sup_down = int(np.sum(y_true_f == -1))
        denom = sup_up + sup_down
        f1_weighted = float((f1_up * sup_up + f1_down * sup_down) / denom) if denom else np.nan

    # Coverage stats always on original (true-nonzero-masked) y_pred
    coverage = float(np.mean(y_pred != 0)) if len(y_pred) else np.nan
    abstain_rate = 1.0 - coverage if not np.isnan(coverage) else np.nan

    return {
        "n_total": int(len(y_true)),
        "n_eval": int(len(y_true_f)),
        "coverage_pred_nonzero": coverage,
        "abstain_rate_pred_zero": abstain_rate,
        "handle_zeros": handle_zeros,

        "tp": tp, "fp": fp, "fn": fn, "tn": tn,

        "accuracy": accuracy,
        "balanced_accuracy": bal_acc,
        "mcc": mcc,

        "precision_up": precision_up_skl,
        "recall_up": recall_up_skl,
        "f1_up": f1_up,

        "precision_down": precision_down,
        "recall_down": recall_down,
        "f1_down": f1_down,

        "f1_macro": f1_macro,
        "f1_weighted": f1_weighted,
    }



# %%

full_perturb_res = {}
correlation_dict = {
    'pearson': [],
    'spearman': [],
    'gene': [],
    'n': []
}
direction_metrics_per_gene = []
logfc_thresholds = [0.25, 0.5, 0.75, 1.0]

for gene in degs_count_filtered.index:
    # flexible regex-based matching for files_h5ad
    # mask can be provided via --mask (e.g. "g{gene}_ssrc_tmask.h5ad", ".*{gene}.*\\.h5ad$", "{}_ssrc.h5ad", "<GENE>_ssrc.h5ad")
    # if not provided, default to matching any filename containing the gene and ending with .h5ad
    if args.model == 'perturbgen':
        adata = read_perturbed_h5ad_files(gene, tf_filtered)
        if adata is None:
            continue
        else:
            pred_deg_df = compute_predicted_DEG_response(adata)

            alpha = 0.05
            logfc_threshold = args.logfc_threshold
            correlation_dict['gene'].append(gene)

            correlation_dict, merged_filtered, merged = compute_correlation_metrics(
                degs_filtered_df[degs_filtered_df['group'] == gene],
                pred_deg_df,
                correlation_dict,
                alpha,
                logfc_threshold
            )
            full_perturb_res[gene] = merged
            for thresh in logfc_thresholds:
                mf = merged_filtered[merged_filtered["l2fc_invitro"].abs() > thresh].copy()

                plot_correlation(
                    mf,
                    correlation_dict,
                    gene,
                    f"{gene}_{args.corr_fig_name}_lfc{thresh}"
                )
                if len(mf) > 0:

                    direction = deg_direction_metrics(
                        mf,
                        true_col="l2fc_invitro",
                        pred_col="l2fc_insilico",
                        handle_zeros="as_negative"  # or "drop"
                    )

                    # store per gene
                    # store one row per (gene, thresh)
                    direction_metrics_per_gene.append({
                        "gene": gene,
                        "model": args.model,
                        "logfc_threshold": thresh,
                        **direction
                    })
                    
                    print(f"\nDirectionality summary for {gene} @ |logFC|>{thresh}")
                    print(f"n_total={direction['n_total']}, n_eval={direction['n_eval']}, coverage={direction['coverage_pred_nonzero']:.3f}")
                    print(f"Accuracy: {direction['accuracy']:.3f}")
                    print(f"Balanced accuracy: {direction['balanced_accuracy']:.3f}")
                    print(f"MCC: {direction['mcc']:.3f}")
                    print(f"Macro-F1: {direction['f1_macro']:.3f}")

            
    elif args.model == 'celloracle':
        deg_results = perturb_celloracle(oracle_data, gene)
        if deg_results is None:
            continue
        else:
            alpha = 0.05
            logfc_threshold = args.logfc_threshold
            correlation_dict['gene'].append(gene)

            correlation_dict, merged_filtered, merged = compute_correlation_metrics(
                degs_filtered_df[degs_filtered_df['group'] == gene],
                deg_results,
                correlation_dict,
                alpha,
                logfc_threshold
            )
            full_perturb_res[gene] = merged
            
            for thresh in logfc_thresholds:
                mf = merged_filtered[merged_filtered["l2fc_invitro"].abs() > thresh].copy()

                plot_correlation(
                    mf,
                    correlation_dict,
                    gene,
                    f"{gene}_{args.corr_fig_name}_lfc{thresh}"
                )
                if len(mf) > 0:

                    direction = deg_direction_metrics(
                        mf,
                        true_col="l2fc_invitro",
                        pred_col="l2fc_insilico",
                        handle_zeros="as_negative"  # or "drop"
                    )

                    # store per gene
                    # store one row per (gene, thresh)
                    direction_metrics_per_gene.append({
                        "gene": gene,
                        "model": args.model,
                        "logfc_threshold": thresh,
                        **direction
                    })
                    
                    print(f"\nDirectionality summary for {gene} @ |logFC|>{thresh}")
                    print(f"n_total={direction['n_total']}, n_eval={direction['n_eval']}, coverage={direction['coverage_pred_nonzero']:.3f}")
                    print(f"Accuracy: {direction['accuracy']:.3f}")
                    print(f"Balanced accuracy: {direction['balanced_accuracy']:.3f}")
                    print(f"MCC: {direction['mcc']:.3f}")
                    print(f"Macro-F1: {direction['f1_macro']:.3f}")

# %%
# concatenate list of dataframes full_perturb_res
for gene, df in full_perturb_res.items():
    df['gene'] = gene
full_perturb_res = pd.concat(full_perturb_res.values(), ignore_index=True)
correlation_df = pd.DataFrame(correlation_dict)
direction_df = pd.DataFrame(direction_metrics_per_gene)
direction_df.to_csv(os.path.join(res_dir, f"{date}_{args.table_name}_direction_per_gene.csv"), index=False)

correlation_df.to_csv(os.path.join(res_dir, f"{date}_{args.table_name}_correlation_per_gene.csv"), index=False)

# %%
# compute precision, recall, f1 score and accuracy for different logfc thresholds
# compared against random baseline

logfc_thresholds = [0.25, 0.5, 0.75, 1.0]
# direction_dict = {}
random_baseline_dict = {}
# for thresh in logfc_thresholds:
#     direction_dict[f'accuracy_top_{thresh}'] = []
#     direction_dict[f'F1_score_top_{thresh}'] = []
#     direction_dict[f'precision_top_{thresh}'] = []
#     direction_dict[f'recall_top_{thresh}'] = []
#     direction_dict[f'MCC_top_{thresh}'] = []
#     random_baseline_dict[f'accuracy_top_{thresh}'] = []
#     random_baseline_dict[f'F1_score_top_{thresh}'] = []
#     random_baseline_dict[f'precision_top_{thresh}'] = []
#     random_baseline_dict[f'recall_top_{thresh}'] = []
#     random_baseline_dict[f'MCC_top_{thresh}'] = []
direction_rows = []   # list of dicts, one per threshold

for thresh in logfc_thresholds:

    full_perturb_res_thresh = full_perturb_res[
        full_perturb_res["l2fc_invitro"].abs() > thresh
    ].copy()

    if len(full_perturb_res_thresh) < 50:
        print(f"Only {len(full_perturb_res_thresh)} genes found for threshold {thresh} in gene {gene}")
        direction_rows.append({
            "gene": gene,
            "logfc_threshold": thresh,
            "n": len(full_perturb_res_thresh),
            "status": "too_few_genes"
        })
        continue

    m = deg_direction_metrics(
        full_perturb_res_thresh,
        true_col="l2fc_invitro",
        pred_col="l2fc_insilico",
        handle_zeros="as_negative"   # or "drop"/"separate"
    )

    # keep the loop, but store *all* metrics
    direction_rows.append({
        "gene": gene,
        "logfc_threshold": thresh,
        **m
    })

# later:
direction_df = pd.DataFrame(direction_rows)
# drop completely empty columns
direction_df = direction_df.dropna(axis=1, how='all')

# columns to keep as identifiers
id_vars = ["gene", "logfc_threshold"]

# everything else is a metric
value_vars = [c for c in direction_df.columns if c not in id_vars]

# long format
long_df = direction_df.melt(
    id_vars=id_vars,
    value_vars=value_vars,
    var_name="metric",
    value_name="value"
)

# save
long_df.to_csv(
    os.path.join(res_dir, f"{date}_{args.table_name}_direction_all.csv"),
    index=False
)

        

        # # call random generator
        # _, random_draws = generate_random_baseline(n_draws=1000, draw_length=len(full_perturb_res_thresh), seed=42)
        # random_accuracies = []
        # random_precisions = []
        # random_recalls = []
        # random_f1s = []
        # for draw in random_draws:

        #     random_yp = (draw == 1).astype(int)
        #     random_accuracy = np.mean(
        #         np.sign(full_perturb_res_thresh['l2fc_invitro']) == np.sign(draw)
        #     )
        #     random_precision = precision_score(
        #         yt,
        #         random_yp,
        #         average='binary'
        #     )
        #     random_recall = recall_score(
        #         yt,
        #         random_yp,
        #         average='binary'
        #     )
        #     random_f1 = f1_score(
        #         yt,
        #         random_yp,
        #         average='binary'
        #     )
        #     random_mcc = matthews_corrcoef(
        #         yt,
        #         random_yp
        #     )
        #     random_accuracies.append(random_accuracy)
        #     random_precisions.append(random_precision)
        #     random_recalls.append(random_recall)
        #     random_f1s.append(random_f1)
        # random_baseline_dict[f'accuracy_top_{thresh}'].append(np.mean(random_accuracies))
        # random_baseline_dict[f'F1_score_top_{thresh}'].append(np.mean(random_f1s))
        # random_baseline_dict[f'precision_top_{thresh}'].append(np.mean(random_precisions))
        # random_baseline_dict[f'recall_top_{thresh}'].append(np.mean(random_recalls))
        # random_baseline_dict[f'MCC_top_{thresh}'].append(np.mean(random_mcc))
        # print('Random summary:')
        # print(f"Accuracy: {np.mean(random_accuracies):.2f}")
        # print(f"Precision: {np.mean(random_precisions):.2f}")
        # print(f"Recall: {np.mean(random_recalls):.2f}")
        # print(f"F1 Score: {np.mean(random_f1s):.2f}")
        # print(f"MCC: {np.mean(random_mcc):.2f}")
        # # raise

# %%
# # create dataframe from direction_dict
# direction_df = pd.DataFrame(direction_dict)
# # drop nan cols
# direction_df = direction_df.dropna(axis=1, how='all')
# # create longdf by splitting col names at '_top_'
# long_df = pd.melt(direction_df,var_name='metric_thresh', value_name='value')
# long_df[['metric', 'logfc_threshold']] = long_df['metric_thresh'].str.rsplit('_top_', n=1, expand=True)
# long_df['logfc_threshold'] = long_df['logfc_threshold'].astype(float)
# long_df = long_df.drop(columns=['metric_thresh'])
# # save long_df
# long_df.to_csv(os.path.join(res_dir, f"{date}_{args.table_name}_direction.csv"), index=False)
# # plot the same for random baseline
# random_baseline_df = pd.DataFrame(random_baseline_dict)
# random_baseline_df = random_baseline_df.dropna(axis=1, how='all')
# random_long_df = pd.melt(random_baseline_df,var_name='metric_thresh', value_name='value')
# random_long_df[['metric', 'logfc_threshold']] = random_long_df['metric_thresh'].str.rsplit('_top_', n=1, expand=True)
# random_long_df['logfc_threshold'] = random_long_df['logfc_threshold'].astype(float)
# random_long_df = random_long_df.drop(columns=['metric_thresh'])

# %%
# # concantenate direction_df and random_baseline_df
# combined_df = pd.concat([long_df.assign(type='model'), random_long_df.assign(type='random')])

# for metric in combined_df['metric'].unique():
#     plt.figure(figsize=(10, 6))
#     sns.barplot(
#         data=combined_df[combined_df['metric'] == metric],
#         x='logfc_threshold', y='value', hue='type'
#     )
#     plt.xlabel('log2fc threshold')
#     plt.ylabel(f'{metric}')
#     # plot legend outside of box
#     plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
#     plt.savefig(os.path.join(res_dir, f"{date}_{args.direction_fig_name}_{metric}.pdf"), bbox_inches='tight', dpi=300)