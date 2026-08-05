# NSCLC radiosensitivity radiogenomics analysis

This repository contains the analysis code accompanying a study of a cell line-anchored radiosensitivity-associated molecular score (CRS), its pathway-level molecular architecture in non-small cell lung cancer, and a computed-tomography radiomic surrogate (RRS).

The repository is organized as a reproducible research-code release rather than a copy of the original interactive analysis directory. Superseded scripts, exploratory tests, temporary workspaces, raw data, patient-level derived data, and machine-specific paths have been removed. The curation criteria are recorded in `docs/code_curation.md`.

## Scope of the code

The code supports the following analyses:

1. construction and repeated nested-cross-validation of the cell-line CRS;
2. application of the CRS to the NSCLC-Radiogenomics transcriptomic cohort;
3. predefined-module and Hallmark pathway analyses;
4. preprocessing of CT images and tumor segmentations;
5. PyRadiomics feature extraction and RRS model development;
6. full-workflow permutation analyses for CRS and RRS performance;
7. tumor-size adjustment of the RRS-CRS association;
8. Lung3 molecular replication;
9. exploratory Lung1 survival analyses and radiomic domain-shift diagnostics;
10. generation of main and supplementary figure source data.

The CRS is a cell line-anchored radiosensitivity-associated molecular phenotype. The RRS is evaluated as a CT radiomic surrogate of that phenotype. The Lung1 analyses are exploratory and do not constitute external clinical validation of the RRS.

## Repository structure

```text
R/
├── 01_cell_line/                 RadioGx acquisition and CRS development
├── 02_nsclc_radiogenomics/       patient transcriptome preparation and CRS calculation
├── 03_radiomics/                 CT metadata, data integration, and RRS development
├── 04_lung1/                     Lung1 application, survival, and domain shift
├── 05_pathways/                  pathway scoring and Lung3 replication
├── 06_sensitivity_analyses/      permutation and tumor-size analyses
├── 07_figures/                   main-figure generation
└── 08_supplementary/             supplementary tables and figure panels

python/
├── 01_nsclc_radiomics/           NSCLC-Radiogenomics feature extraction
├── 02_lung1_preprocessing/       DICOM, SEG, NRRD, and tumor-mask processing
└── 03_lung1_radiomics/           Lung1 feature extraction

environment/                      dependency installation and environment capture
config/                           environment-variable example
data/                             data-access documentation only
docs/                             workflow and release documentation
tools/                            repository validation
```

## Configuration

The code repository and the analysis workspace are intentionally separate. Set the following environment variables before running the analysis:

```text
RADIOGENOMICS_PROJECT_DIR=/absolute/path/to/analysis-workspace
LUNG1_DATA_DIR=/absolute/path/to/lung1-workspace
```

R example:

```r
Sys.setenv(
  RADIOGENOMICS_PROJECT_DIR = "/absolute/path/to/analysis-workspace",
  LUNG1_DATA_DIR = "/absolute/path/to/lung1-workspace"
)
```

Shell example:

```bash
export RADIOGENOMICS_PROJECT_DIR=/absolute/path/to/analysis-workspace
export LUNG1_DATA_DIR=/absolute/path/to/lung1-workspace
```

The expected workspace structure is documented in `docs/workspace_layout.md`.

## Software dependencies

Install the R dependencies from the repository root:

```r
source("environment/install_r_dependencies.R")
```

Install the standard Python dependencies with:

```bash
python -m pip install -r environment/requirements.txt
```

Several Lung1 segmentation scripts are designed to run in the Python environment bundled with 3D Slicer. These scripts import `slicer`, `vtk`, or `DICOMLib` and should not be run in a standard Python interpreter. The analysis used 3D Slicer 5.12.2.

After completing the analysis, capture the software environment with:

```r
source("environment/capture_environment.R")
```

```bash
python environment/capture_environment.py
```

## Execution order

The full execution order and the principal inputs and outputs are listed in `docs/pipeline.md`. Scripts are numbered within each analysis stage. Figure scripts should be run only after their corresponding analysis outputs have been generated.

The computationally intensive permutation scripts are resume-safe. The RRS permutation analysis recovers the archived outer-fold assignments, defines reproducible inner folds, and uses the same frozen split system for the observed and permuted analyses. The prespecified primary permutation test uses the positive upper tail; a doubled-tail two-sided sensitivity result is also reported.

## Data availability

No raw or derived participant-level data are included. Public data sources and access requirements are summarized in `data/README.md`. The source repositories retain control over their data-use terms.

## Repository validation

Run the static repository checks from the repository root:

```bash
python tools/validate_repository.py
```

The validation script checks repository structure, English-only text and filenames, excluded file types, local absolute paths, Python syntax, basic R delimiter balance, dependency declarations, and common release artifacts. A successful static check does not replace end-to-end execution in the original data environment. The completed static review is summarized in `VALIDATION_REPORT.md`; `MANIFEST.tsv` records file sizes and SHA-256 checksums for the release contents.

## Releasing the code

The recommended manuscript-release workflow is:

1. create a public GitHub repository;
2. complete and rename `CITATION.cff.template` to `CITATION.cff`;
3. complete `.zenodo.json.template` and `CODE_AVAILABILITY.md`;
4. run the repository validator and the full analysis in the author environment;
5. create a tagged GitHub release, such as `v1.0.0`;
6. archive that release through the Zenodo GitHub integration;
7. cite the version-specific DOI in the manuscript.

Detailed release instructions are in `docs/github_zenodo_release.md`.

## License

The analysis code is released under the MIT License. The datasets, software packages, and third-party resources used by the analysis remain subject to their own licenses and terms of use.
