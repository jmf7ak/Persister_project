#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard
#SBATCH --output=/scratch/jmf7ak/persister/outputs/ATP_open/untreated_24B_81017d/riptide_output.log
#SBATCH --error=/scratch/jmf7ak/persister/outputs/ATP_open/untreated_24B_81017d/riptide_error.log

module load miniforge
module load bioconda
source activate riptide_cobra_311
module load gurobi/10.0.1
export PYTHONPATH=/lib/python3.11_utf32
python run_riptide_open_ATP.py "/scratch/jmf7ak/persister/counts_batchcorrected/untreated_24B_81017d.csv" "/scratch/jmf7ak/persister/iPau21.xml" "/scratch/jmf7ak/persister/outputs/ATP_open/untreated_24B_81017d"
