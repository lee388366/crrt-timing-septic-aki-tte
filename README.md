# CRRT Timing in Septic Shock-Associated Severe AKI — Target Trial Emulation

Analysis and figure-generation code for the study:

**Early Continuous Renal Replacement Therapy Initiation and 28-Day Mortality in Critically Ill
Patients with Septic Shock and Severe Acute Kidney Injury: A Target Trial Emulation Using
Clone-Censor-Weighting and Inverse Probability of Censoring Weights.**

Menglin Liang, Chen Li (co-first authors), Zhi Mao, Feihu Zhou (corresponding).

## What this repository contains

- `analysis/` — R and SQL for cohort extraction, covariate build, clone-censor-weighting with
  stabilized IPCW (CCW-IPCW), Cox/RMST primary analysis, competing-risk renal-recovery models,
  sensitivity/bias analyses (E-value, tipping point, anchor/grace-window/era/rescue), and the
  prospectively registered Bayesian reanalysis.
- `figures/` — R scripts that render main-text Figures 1–3 and supplementary figures.

**No patient-level data are included.** The study uses MIMIC-IV v3.1, which is available only to
credentialed users through PhysioNet (https://physionet.org/content/mimiciv/3.1/) under its data
use agreement. All `.csv`/`.rds`/`.xlsx` and other derived data are intentionally excluded.

## Reproducing

1. Obtain MIMIC-IV v3.1 via PhysioNet and load it into a local PostgreSQL database (default here:
   database `mimiciv31` on `localhost:5432`, user `postgres`).
2. Set the database password via environment variable before running R:
   ```bash
   export PGPASSWORD=your_local_password
   ```
3. Run scripts from the repository root. Paths are relative to the repo root; the pipeline writes
   intermediate objects to `analysis/data/` (not tracked).

## Registration

The Bayesian reanalysis was prospectively registered on the Open Science Framework
(DOI: 10.17605/OSF.IO/BCD52).

## License

MIT License — see `LICENSE`.
