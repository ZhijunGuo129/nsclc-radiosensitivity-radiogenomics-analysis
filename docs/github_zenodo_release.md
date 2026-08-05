# GitHub and Zenodo release workflow

## Prepare the repository

1. Create an empty public GitHub repository.
2. Copy the complete contents of this code package into the repository root.
3. Complete `CITATION.cff.template`, rename it to `CITATION.cff`, and remove the template file.
4. Complete `.zenodo.json.template`, rename it to `.zenodo.json`, and remove the template file.
5. Replace the placeholders in `CODE_AVAILABILITY.md` after the repository URL and DOI are known.
6. Confirm that no data, credentials, or machine-specific files are staged for commit.

## Validate before the first commit

```bash
python tools/validate_repository.py
python -m compileall -q python environment tools
```

Run the complete R and Python workflow in the original analysis environment. Compare the generated numerical results with the submitted manuscript and retain the captured session information.

## Create the manuscript release

1. Commit the reviewed repository state.
2. Push the default branch to GitHub.
3. Enable the repository in the Zenodo GitHub integration.
4. Create an annotated Git tag for the manuscript version:

```bash
git tag -a v1.0.0 -m "Manuscript submission release"
git push origin v1.0.0
```

5. Create a GitHub Release from tag `v1.0.0`.
6. Confirm that Zenodo has archived the release and assigned a version-specific DOI.
7. Add the DOI badge and version DOI to the GitHub README and release notes.
8. Insert the version-specific DOI in the manuscript Code Availability statement.

## Versioning policy

Do not alter an archived release. Corrections should be committed to the repository and issued under a new semantic version, for example `v1.0.1` or `v1.1.0`. The manuscript should cite the exact release used for the reported analysis.
