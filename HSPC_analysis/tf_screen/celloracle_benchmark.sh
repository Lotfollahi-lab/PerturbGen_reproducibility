#make a date directory if it does not exist
#!/bin/bash
#BSUB -q normal # run CPU job
#BSUB -n 2 # number of cores
#BSUB -G team361 # groupname for billing
#BSUB -cwd /lustre/scratch126/cellgen/lotfollahi/kl11/ # working directory
#BSUB -o TRACE-reproducibility/logs/celloracle_benchmark_%J.out # output file
#BSUB -e TRACE-reproducibility/logs/celloracle_benchmark_%J.err # error file
#BSUB -M 75000  # RAM memory part 2. Default: 100MB
#BSUB -R 'select[mem>75000] rusage[mem=75000]' # RAM memory part 1. Default: 100MB
#BSUB -J celloracle_benchmark # job name

# activate python environment
module load HGI/common/conda/module
conda init bash && source ~/.bashrc
conda deactivate
conda activate /software/cellgen/team361/am74/envs/celloracle
cwd=$(pwd)

echo "--- Start running CellOracle"

python3 $cwd/TRACE-reproducibility/HSPC_TF_screen/4_plot_benchmark.py \
--model 'celloracle' \
--oracle_path '/lustre/scratch126/cellgen/lotfollahi/kl11/T_perturb/res/hspc/celloracle/memp_prog_baeoma.celloracle.oracle' \
--tf_excl_celltype 'Early GMP' 'MEP' 'EoBasoMast Precursor' \
--hspc_incl_celltype 'Early_Ery' 'Late_Ery' 'BaEoMa' 'MK' \
--corr_fig_name 'celloracle_5k_n5_mouse_ATAC_mk_erythro_baeoma_logfc_correlation' \
--confusion_fig_name 'celloracle_5k_n5_mouse_ATAC_mk_erythro_baeoma_confusion_mtx' \
--table_name 'celloracle_5k_n5_mouse_ATAC_mk_erythro_baeoma_metrics' \
--direction_fig_name 'celloracle_5k_n5_mouse_ATAC_mk_erythro_baeoma_direction' \
--res_dir 'T_perturb/res/hspc/celloracle/' \
--n_prop 5
echo '--- Finished plotting benchmark'
