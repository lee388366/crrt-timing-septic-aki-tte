# CRRT timing target trial emulation

## Current BMJ Open analysis source snapshot — 17 October 2026

The current source package is in [bmj_open_2026_10_17](bmj_open_2026_10_17/). Earlier root-level scripts and the July 2026 commit remain available for provenance; they are not the current numerical source.

This version publishes the 146-file source snapshot reviewed for the revised manuscript. The corresponding local ZIP SHA-256 is `16b33c77062c04ef312f46e692ef9f71caa5c42b2e75b21851df713bedf5a202`. Every file in that ZIP is represented in the versioned directory. The package README's “LOCAL REVIEW CANDIDATE” and “not uploaded” statements record its build-time status; this root notice records publication of that unchanged snapshot, not a claim of a fully portable validated release.

The principal correction converts the vasopressin dose term by dividing by 60 before applying the existing coefficient, retaining the original cohort, eligibility times, strategies and estimation specifications. The revised primary weighted HR is 1.05144237278742 (95% CI 0.92985279471207–1.18893126909961); the full-cohort HR is 1.11186346301827 (1.00869995484273–1.22557788811216). Source changes and object mappings are documented within the versioned package. Use its current sources and matching private inputs rather than historical estimates embedded in legacy scripts.

### Verification and limits

Original local branches underwent old-result reproduction before the minimal correction and revised runs. The packaged portable adapters passed static parsing and source-equivalence checks; a fresh-environment end-to-end execution has **not** been performed. This repository does not claim that all figures and tables can be regenerated with one command. Complete offline reproduction of the rich/FULL Table 1 and overall Table S1 remains limited by missing historical augmentation inputs or generators; their documented NEE rows have separate verification. See the package README, object/source map, patch log and environment records for details.

### Data access

This is a code-only repository. No patient-level MIMIC-IV/eICU data, clinical identifiers, database credentials or private model objects are distributed. Access to the source databases requires the applicable PhysioNet credentials and data use agreements. Configure private input paths separately and do not commit generated clinical data.
