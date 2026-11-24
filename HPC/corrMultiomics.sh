#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=corrMO
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ks575@exeter.ac.uk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=corrMO.all_%j.log
#SBATCH --account=m.d.sharma

set -euo pipefail
# One BLAS thread per Python block worker avoids nested oversubscription.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

module purge
module load slurm/18.08.4 DefaultModules shared Workspace/v1
module load Mamba/Python3
mamba activate /nobackup/beegfs/home/ISAD/ks575/meantools

work_dir=/nobackup/beegfs/workspace/ks575/Data/WUR_Mar/nwo-groot/coexpression/mate
script_dir=/nobackup/beegfs/workspace/ks575/Data/software/meantools
# Install the full checkout so the runner can import the shared mutual_ranks.py.
python_script="$script_dir/HPC/corrMultiomics_HPC.py"
if [[ ! -f "$python_script" ]]; then
    python_script="$script_dir/corrMultiomics_HPC.py"
fi
scratch_base="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
job_scratch=$(mktemp -d "$scratch_base/corrMO.${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf -- "$job_scratch"' EXIT
export SQLITE_TMPDIR="$job_scratch"

# Select at submission time, e.g. sbatch --export=ALL,MR_STORAGE=ram,MR_MEMORY_GB=20 corrMultiomics.sh
mr_args=(--mr-storage "${MR_STORAGE:-disk}")
if [[ -n "${MR_MEMORY_GB:-}" ]]; then
    mr_args+=(--mr-memory-gb "$MR_MEMORY_GB")
fi

printf 'Job=%s Node=%s CPUs=%s Start=%s\n' "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST" "$SLURM_CPUS_PER_TASK" "$(date -Is)"
# Use the same file for both inputs to construct a unimodal coexpression network.
# For metabolite vs transcript analysis, change only -ft to the metabolite matrix.
# MR construction requires a new database, as in the standard runner.
srun --cpu-bind=cores /usr/bin/time -v python -u "$python_script" \
    -ft "$work_dir/MATE_test1_transcripts_log2z_filtered.csv" \
    -qm "$work_dir/MATE_test1_transcripts_log2z_filtered.csv" \
    -dn "$work_dir/drosophila_TR_${SLURM_JOB_ID}.sqlite" \
    -tn "droso_coex_${SLURM_JOB_ID}" \
    -mad -r -mr -cl -c 0.2 -mdr 10 25 \
    -t "$SLURM_CPUS_PER_TASK" "${mr_args[@]}" --block-size 128 \
    --scratch-dir "$job_scratch" \
    --clusterone-jar "$script_dir/cluster_one-1.0.jar" --java-heap 8g -v
