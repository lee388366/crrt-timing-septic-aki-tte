# Early versus deferred CRRT in vasopressor-supported sepsis with KDIGO stage 3 AKI: target trial emulation

Analysis code for the manuscript *Early Continuous Renal Replacement Therapy Initiation and 28-Day All-Cause Mortality in Vasopressor-Supported Sepsis with Physiologic KDIGO Stage 3 Acute Kidney Injury: A Target Trial Emulation* (BMJ Open, under consideration).

## Repository layout

| Folder | Contents |
|---|---|
| [`revision_2026_11/`](revision_2026_11/) | **Current analysis.** Corrected clone-censor-weighting implementation, all estimates, tables and figures reported in the revised manuscript and Additional file 1. Start here. |
| [`upstream/`](upstream/) | Scripts that extracted the cohort and built the private input datasets (baseline covariates, hourly covariates, alternative-eligibility and urgent-indication branches, secondary-outcome dataset, Table 1 queries) read by `revision_2026_11/`. Estimation steps inside these scripts belong to the superseded implementation and are not used for reported results. |

The code as cited in the original October 2026 submission is preserved under the tag [`v2026.10.17-submission`](../../tree/v2026.10.17-submission). That version contained errors in the weighting implementation (time-fixed weights, an early-strategy weight model that could use post-initiation covariates, and baseline measurements recorded after time zero); the corrections and all changes are listed in Supplemental Table S15 of the revised manuscript and in [`revision_2026_11/README.md`](revision_2026_11/README.md).

## Data access

This is a code-only repository. No patient-level MIMIC-IV data, identifiers, database credentials or patient-level model objects are distributed. MIMIC-IV v3.1 is available to credentialed users through PhysioNet (https://physionet.org/content/mimiciv/3.1/) under its data use agreement. Scripts read local paths from environment variables (see `revision_2026_11/README.md`); patient-level intermediate files are written only to a local `private_data/` directory.

## Licence

MIT (see `LICENSE`).
