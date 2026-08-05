# Release checklist

## Repository content

- [ ] All filenames and tracked text are English.
- [ ] No raw imaging, NRRD, NIfTI, RDS, spreadsheet result files, patient-level tables, or model objects are tracked.
- [ ] No credentials, access tokens, local absolute paths, temporary files, caches, or editor files are tracked.
- [ ] Only the final analysis and figure scripts are included.
- [ ] Script names and the pipeline documentation agree.

## Technical validation

- [ ] `python tools/validate_repository.py` completes without errors.
- [ ] `python -m compileall -q python environment tools` completes without errors.
- [ ] All R scripts parse and execute in the original R environment.
- [ ] 3D Slicer scripts execute in the documented Slicer version.
- [ ] Captured R and Python environment files are retained with the release records.

## Scientific validation

- [ ] Full-workflow permutation results agree with the manuscript.
- [ ] RRS AUC and DeLong confidence interval are reproduced from patient-level out-of-fold predictions.
- [ ] Tumor-size-adjusted results agree with the manuscript.
- [ ] Lung3 all-Hallmark and selected-Hallmark summaries agree with the manuscript.
- [ ] Lung1 primary, signed-log, and stage-adjusted Cox results agree with the manuscript.
- [ ] Obsolete TNM-adjusted results are absent.
- [ ] Figure and supplementary source data are regenerated from analysis outputs.

## Public release

- [ ] Author names, ORCID identifiers, repository URL, and manuscript metadata are entered in `CITATION.cff`.
- [ ] Creator metadata are entered in `.zenodo.json`.
- [ ] The GitHub repository is enabled in the Zenodo GitHub integration.
- [ ] A GitHub release tagged `v1.0.0` is created from the reviewed commit.
- [ ] Zenodo archives the release and assigns a version-specific DOI.
- [ ] The DOI is entered in `CODE_AVAILABILITY.md` and the manuscript.
