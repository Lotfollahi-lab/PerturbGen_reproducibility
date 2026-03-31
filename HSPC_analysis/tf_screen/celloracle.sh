#make a date directory if it does not exist
#!/bin/bash
#BSUB -q normal # run CPU job
#BSUB -n 4 # number of cores
#BSUB -G team361 # groupname for billing
#BSUB -cwd /lustre/scratch126/cellgen/lotfollahi/kl11/ # working directory
#BSUB -o TRACE-reproducibility/logs/celloracle_test_%J.out # output file
#BSUB -e TRACE-reproducibility/logs/celloracle_test_%J.err # error file
#BSUB -M 50000  # RAM memory part 2. Default: 100MB
#BSUB -R 'select[mem>50000] rusage[mem=50000]' # RAM memory part 1. Default: 100MB
#BSUB -J celloracle_test # job name

# activate python environment
module load HGI/common/conda/module
conda init bash && source ~/.bashrc
conda deactivate
conda activate /software/cellgen/team361/am74/envs/celloracle
cwd=$(pwd)

echo "--- Start running CellOracle"

python3 $cwd/TRACE-reproducibility/HSPC_TF_screen/2.2_celloracle_grn.py \
--data_path '/lustre/scratch126/cellgen/lotfollahi/kl11/T_perturb/res/hspc/celloracle/full_data.h5ad' \
--atac_grn_path '/nfs/team361/am74/Cytomeister/notebooks/gse96769_HSPC_ATAC/base_GRN_dataframe.parquet' \
--celltype 'Early_Ery' 'Late_Ery' 'MK' 'BaEoMa' \
--n_jobs 4 \
--file_name '5k_n5_human_ATAC_mk_erythro_baeoma' \
--outdir '/lustre/scratch126/cellgen/lotfollahi/kl11/T_perturb/res/hspc/celloracle' \

echo '--- Finished running CellOracle'
