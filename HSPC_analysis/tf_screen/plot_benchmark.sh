#make a date directory if it does not exist
#!/bin/bash
#BSUB -q normal # run CPU job
#BSUB -n 1 # number of cores
#BSUB -G team361 # groupname for billing
#BSUB -cwd /lustre/scratch126/cellgen/lotfollahi/kl11/ # working directory
#BSUB -o TRACE-reproducibility/logs/plot_benchmark_%J.out # output file
#BSUB -e TRACE-reproducibility/logs/plot_benchmark_%J.err # error file
#BSUB -M 20000  # RAM memory part 2. Default: 100MB
#BSUB -R 'select[mem>20000] rusage[mem=20000]' # RAM memory part 1. Default: 100MB
#BSUB -J plot_benchmark # job name

# activate python environment
source /nfs/team361/cytomeister/.cytomeister/bin/activate
cwd=$(pwd)

echo "--- Start plotting benchmark"

python3 $cwd/TRACE-reproducibility/HSPC_TF_screen/4_plot_benchmark.py \
--model 'perturbgen' \
--h5ad_files 'T_perturb/res/hspc/perturbation_5k_delete' \
--tf_excl_celltype 'Early GMP' 'MEP' 'EoBasoMast Precursor' \
--hspc_incl_celltype 'Early_Ery' 'Late_Ery' 'MK' 'BaEoMa' \
--corr_fig_name 'intermediate_delete_ep14_memp_prog_logfc_correlation' \
--confusion_fig_name 'intermediate_delete_ep14_memp_prog_confusion_mtx' \
--table_name 'intermediate_delete_ep14_memp_prog_metrics' \
--direction_fig_name 'intermediate_delete_ep14_memp_prog' \
--res_dir 'T_perturb/res/hspc/tf_perturbation/'
echo '--- Finished plotting benchmark'
