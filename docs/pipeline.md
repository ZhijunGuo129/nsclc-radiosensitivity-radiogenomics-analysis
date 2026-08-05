# Analysis pipeline

Set the environment variables in `config/environment.example` before running any analysis script. The sequence below assumes that source data have been placed in the workspace described in `docs/workspace_layout.md`.

## Stage 1: cell-line CRS development

1. `R/01_cell_line/01_download_radiogx_data.R`
2. `R/01_cell_line/02_prepare_radiogx_rnaseq.R`
3. `R/01_cell_line/03_train_crs_model.R`
4. `R/06_sensitivity_analyses/03_permute_crs_model.R`

The permutation analysis repeats the complete nested model-building process, including training-fold variable-gene selection and elastic-net tuning.

## Stage 2: patient molecular phenotype

1. `R/02_nsclc_radiogenomics/01_prepare_transcriptome.R`
2. `R/02_nsclc_radiogenomics/02_map_crs_features.R`
3. `R/02_nsclc_radiogenomics/03_compute_patient_scores.R`

## Stage 3: NSCLC-Radiogenomics imaging and RRS development

1. `R/03_radiomics/01_prepare_patient_identifiers.R`
2. `R/03_radiomics/02_query_tcia_metadata.R`
3. `R/03_radiomics/03_verify_dicom_download.R`
4. `R/03_radiomics/04_build_ct_seg_pairs.R`
5. `python/01_nsclc_radiomics/extract_radiomic_features.py`
6. `R/03_radiomics/05_merge_radiomics_and_molecular_data.R`
7. `R/03_radiomics/06_train_rrs_model.R`
8. `R/03_radiomics/07_summarize_training_cohort.R`
9. `R/06_sensitivity_analyses/04_permute_rrs_model.R`
10. `R/06_sensitivity_analyses/02_analyze_tumor_size_adjustment.R`

The RRS permutation analysis uses 20 repeated outer five-fold splits and inner five-fold tuning. The outer assignments are recovered from the archived training outputs. New reproducible inner-fold assignments are frozen and used identically for the observed and all permuted analyses.

## Stage 4: Lung1 preprocessing and exploratory analyses

1. `R/04_lung1/01_prepare_clinical_data.R`
2. `python/02_lung1_preprocessing/01_inventory_dicom_series.py`
3. `python/02_lung1_preprocessing/02_convert_dicom_seg_to_nrrd.py`
4. `python/02_lung1_preprocessing/03_validate_nrrd_masks.py`
5. `python/02_lung1_preprocessing/04_extract_primary_tumor_masks.py`
6. `python/02_lung1_preprocessing/05_prepare_mask_selection.py`
7. `python/02_lung1_preprocessing/06_export_selected_masks_in_slicer.py`
8. `python/02_lung1_preprocessing/07_validate_tumor_masks.py`
9. `python/03_lung1_radiomics/extract_radiomic_features.py`
10. `R/04_lung1/02_apply_rrs_model.R`
11. `R/04_lung1/03_analyze_overall_survival.R`
12. `R/04_lung1/04_assess_domain_shift.R`
13. `R/04_lung1/05_train_signed_log_sensitivity_model.R`
14. `R/04_lung1/06_analyze_stage_adjusted_survival.R`
15. `R/06_sensitivity_analyses/01_prepare_tumor_size_data.R`

Scripts 02 and 06 in the Lung1 Python preprocessing stage require 3D Slicer Python. The signed-log model is a scale-sensitivity analysis and is not a second primary model.

## Stage 5: pathway analyses and Lung3 replication

1. `R/05_pathways/01_score_predefined_modules.R`
2. `R/05_pathways/02_score_hallmark_pathways.R`
3. `R/05_pathways/03_prepare_lung3_data.R`
4. `R/05_pathways/04_compute_lung3_crs.R`
5. `R/05_pathways/05_evaluate_lung3_replication.R`

## Stage 6: figures and supplementary files

Main figures:

1. `R/07_figures/01_create_crs_figure.R`
2. `R/07_figures/02_create_molecular_replication_figure.R`
3. `R/07_figures/03_create_rrs_figure.R`
4. `R/07_figures/04_create_lung1_figure.R`

Supplementary files:

1. `R/08_supplementary/01_export_supplementary_tables.R`
2. `R/08_supplementary/02_create_supplementary_figure_panels.R`

Figure scripts read analysis outputs and write separate source-data files. They do not substitute for the inferential analysis scripts.
