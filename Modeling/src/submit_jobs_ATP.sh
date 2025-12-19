#!/bin/bash

COUNTS_DIR="/scratch/jmf7ak/persister/counts_batchcorrected"
MODEL_FILE="/scratch/jmf7ak/persister/iPau21.xml"
OUTPUT_BASE="/scratch/jmf7ak/persister/outputs/ATP_open"


for COUNTS_FILE in "$COUNTS_DIR"/*.csv; do
    SAMPLE_NAME=$(basename "$COUNTS_FILE" .csv)
    OUTPUT_DIR="$OUTPUT_BASE/$SAMPLE_NAME"
    mkdir -p "$OUTPUT_DIR"

    # Create a separate SBATCH script for each job
    JOB_SCRIPT="${OUTPUT_DIR}/run_riptide_${SAMPLE_NAME}.sh"

    cat <<EOT > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard
#SBATCH --output=${OUTPUT_DIR}/riptide_output.log
#SBATCH --error=${OUTPUT_DIR}/riptide_error.log

module load miniforge
module load bioconda
source activate riptide_cobra_311
module load gurobi/10.0.1
export PYTHONPATH=$EBROOTGUROBI/lib/python3.11_utf32
python run_riptide_open_ATP.py "$COUNTS_FILE" "$MODEL_FILE" "$OUTPUT_DIR"
EOT

    # Submit the job
    sbatch "$JOB_SCRIPT"
done


echo "Submitted Riptide jobs for all counts files in $COUNTS_DIR"

