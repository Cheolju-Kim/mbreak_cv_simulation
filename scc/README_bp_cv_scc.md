# BU SCC Workflow for mbreak Critical-Value Simulation

This follows the older SCC workflow used in `C:/Users/econp/Desktop/MRT/ssc`.

## 1. Local

Commit and push the current project to GitHub.

```bash
git add ver_1.0.3/cv_simulation.R ver_1.0.3/run_bp_cv_simulation.R ver_1.0.3/scc
git commit -m "Add BP critical value simulation runner"
git push
```

## 2. SCC Setup

On SCC, clone or pull the repository.

```bash
cd ~
mkdir -p mbreak_cv_work
cd mbreak_cv_work
git clone <your-repo-url>
cd <repo>/ver_1.0.3
```

If already cloned:

```bash
cd ~/mbreak_cv_work/<repo>/ver_1.0.3
git pull
```

Load R and test the runner.

```bash
module load R/4.3.1
Rscript run_bp_cv_simulation.R --mode=smoke --rep=20 --n-grid=100 --cores=2 --fit=false
```

## 3. Submit Jobs

Run the BP 1998 Table I/II replication.

```bash
qsub scc/run_bp_cv_simulation.qsub.sh
```

The qsub script defaults to:

```bash
BP_CV_MODE=bp
BP_CV_REP=10000
BP_CV_NGRID=1000
BP_CV_CORES=$NSLOTS
```

For the final expanded mbreak response-surface simulation:

```bash
qsub -v BP_CV_MODE=final,BP_CV_SEED=303 scc/run_bp_cv_simulation.qsub.sh
```

The `final` grid uses `q = 1:20, 25, 30, ..., 60, 70, 80, 90, 100`,
`trm = 0.05, 0.075, ..., 0.30`, `k_cap = 10`, `l_cap = 10`, and
`double_m_cap = 5`.

Because BU SCC jobs should stay within 24 hours, split the full `final` run by
q ranges rather than submitting the entire grid as one job. For example:

```bash
qsub -v BP_CV_MODE=final,BP_CV_REP=10000,BP_CV_NGRID=1000,BP_CV_SEED=303,BP_CV_Q_MIN=1,BP_CV_Q_MAX=20 scc/run_bp_cv_simulation.qsub.sh
qsub -v BP_CV_MODE=final,BP_CV_REP=10000,BP_CV_NGRID=1000,BP_CV_SEED=304,BP_CV_Q_MIN=25,BP_CV_Q_MAX=40 scc/run_bp_cv_simulation.qsub.sh
qsub -v BP_CV_MODE=final,BP_CV_REP=10000,BP_CV_NGRID=1000,BP_CV_SEED=305,BP_CV_Q_MIN=45,BP_CV_Q_MAX=60 scc/run_bp_cv_simulation.qsub.sh
qsub -v BP_CV_MODE=final,BP_CV_REP=10000,BP_CV_NGRID=1000,BP_CV_SEED=306,BP_CV_Q_MIN=70,BP_CV_Q_MAX=70 scc/run_bp_cv_simulation.qsub.sh
qsub -v BP_CV_MODE=final,BP_CV_REP=10000,BP_CV_NGRID=1000,BP_CV_SEED=307,BP_CV_Q_MIN=80,BP_CV_Q_MAX=80 scc/run_bp_cv_simulation.qsub.sh
qsub -v BP_CV_MODE=final,BP_CV_REP=10000,BP_CV_NGRID=1000,BP_CV_SEED=308,BP_CV_Q_MIN=90,BP_CV_Q_MAX=90 scc/run_bp_cv_simulation.qsub.sh
qsub -v BP_CV_MODE=final,BP_CV_REP=10000,BP_CV_NGRID=1000,BP_CV_SEED=309,BP_CV_Q_MIN=100,BP_CV_Q_MAX=100 scc/run_bp_cv_simulation.qsub.sh
```

For a shorter pilot final run:

```bash
qsub -v BP_CV_MODE=final,BP_CV_REP=1000,BP_CV_NGRID=500,BP_CV_SEED=303 scc/run_bp_cv_simulation.qsub.sh
```

## 4. Monitor and Results

Check the queue.

```bash
qstat -u $USER
```

When the job finishes, outputs are under:

```text
cv_output/scc/
```

Useful files:

```text
cv_<mode>.csv
cv_<mode>.rds
bp1998_table_comparison.csv
bp1998_table_comparison_summary.csv
response_surface_coefficients_<mode>.csv
response_surface_fits_<mode>.rds
```

The qsub log is created in the submission directory and is named from the job
name, for example:

```text
mbreak_cv.o<jobid>
```
