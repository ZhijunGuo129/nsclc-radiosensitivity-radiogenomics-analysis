# Validation Report

Release: `v1.0.0`

## Repository curation

The complete source archives were reviewed before this release was assembled. The release retains 34 R analysis scripts and 9 Python scripts. R workspace files, history files, test-only scripts, superseded figure revisions, obsolete low-permutation analyses, temporary utilities, raw data, derived participant-level files, and machine-specific settings were excluded.

## Static checks completed

- Repository filenames are ASCII and English.
- Tracked text contains no CJK characters.
- No raw imaging files, NRRD or NIfTI files, R workspaces, spreadsheets, result CSV files, archives, caches, or model objects are included.
- No local workstation paths are present in R or Python source files.
- Python source files pass abstract-syntax-tree parsing and byte-code compilation.
- R source files pass delimiter and quoted-string balance checks.
- R package references are covered by `environment/install_r_dependencies.R`.
- Standard Python imports are covered by `environment/requirements.txt`; 3D Slicer-only imports are documented separately.
- Script paths listed in `docs/pipeline.md` exist.
- The obsolete absolute-correlation calculation for the RRS two-sided permutation result is absent. The release uses a doubled-tail calculation based on the smaller empirical tail.
- GitHub Actions runs the repository validator and Python compilation on pushes and pull requests.

## Scientific code review points

- The CRS permutation script repeats the complete nested model-building workflow, including training-fold variable-gene selection.
- The RRS permutation script uses the recovered outer-fold assignments and a reproducible frozen inner-fold system for both observed and permuted analyses.
- The primary RRS permutation test is the positive upper-tail test.
- Main figure scripts read analysis outputs rather than substituting manuscript values for the underlying analyses.
- Lung1 stage adjustment uses the auditable `Overall.Stage` complete-case model. No obsolete TNM-adjusted result is included.
- The signed-log RRS remains a scale-sensitivity analysis rather than a co-primary model.

## Checks not performed in this build environment

The build environment did not contain R, the public source datasets, DICOM images, or 3D Slicer. Therefore, the complete R workflow, imaging conversion, radiomic extraction, and numerical reproduction were not executed during package assembly. The authors must run the entire pipeline in the original analysis environment and compare all outputs with the submitted manuscript before creating the public release.

## Release metadata still required

Before publication, the authors must complete `CITATION.cff.template`, `.zenodo.json.template`, and `CODE_AVAILABILITY.md` with the final author list, ORCID identifiers, affiliations, GitHub URL, and version-specific Zenodo DOI.
