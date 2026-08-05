# Analysis workspace layout

The public code repository does not contain the analysis data. Set `RADIOGENOMICS_PROJECT_DIR` to the main workspace and `LUNG1_DATA_DIR` to the Lung1 workspace.

```text
RADIOGENOMICS_PROJECT_DIR/
├── 01_raw_data/
├── 02_metadata/
├── 03_processed_images/
├── 04_radiomics_features/
├── 05_molecular_scores/
├── 06_models/
├── 07_results/
├── 08_figures_final/
└── 10_logs/

LUNG1_DATA_DIR/
├── TCIA_DICOM/
├── clinical/
├── metadata/
├── NRRD_full/
├── NRRD_tumor_only_fast/
├── NRRD_tumor_only_final/
├── radiomics_features/
├── external_RRS/
├── survival_validation/
├── domain_shift_diagnostics/
└── signed_log_rrs/
```

The scripts create result directories when needed. Source-data directories and filenames should be retained as documented by the corresponding public repositories.

Do not place credentials, access tokens, raw imaging files, participant-level tables, or local environment files in the public code repository. The repository `.gitignore` excludes common data, image, model, workspace, cache, and log formats.
