#!/bin/bash -l
#$ -N mbreak_cv
#$ -cwd
#$ -j y
#$ -V
#$ -l h_rt=24:00:00
#$ -pe omp 14
#$ -m ea
#$ -M cheolju@bu.edu

module load R/4.3.1

mkdir -p cv_output/scc

export BP_CV_MODE="${BP_CV_MODE:-bp}"
export BP_CV_REP="${BP_CV_REP:-10000}"
export BP_CV_NGRID="${BP_CV_NGRID:-1000}"
export BP_CV_CORES="${BP_CV_CORES:-${NSLOTS:-1}}"
export BP_CV_SEED="${BP_CV_SEED:-202}"
export BP_CV_FIT="${BP_CV_FIT:-true}"
export BP_CV_OUT="${BP_CV_OUT:-cv_output/scc}"

extra_args=()
if [ -n "${BP_CV_Q_MIN:-}" ]; then
  extra_args+=(--q-min="${BP_CV_Q_MIN}")
fi
if [ -n "${BP_CV_Q_MAX:-}" ]; then
  extra_args+=(--q-max="${BP_CV_Q_MAX}")
fi

Rscript run_bp_cv_simulation.R \
  --mode="${BP_CV_MODE}" \
  --rep="${BP_CV_REP}" \
  --n-grid="${BP_CV_NGRID}" \
  --cores="${BP_CV_CORES}" \
  --seed="${BP_CV_SEED}" \
  --fit="${BP_CV_FIT}" \
  --out="${BP_CV_OUT}" \
  "${extra_args[@]}"
