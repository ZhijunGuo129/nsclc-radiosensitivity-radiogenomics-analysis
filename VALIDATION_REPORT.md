# Validation Report

## Repository-level validation

The public repository passed GitHub Actions static validation. The checks included required repository files and directories, repository metadata, absence of machine-specific absolute paths, absence of placeholder metadata, absence of CJK characters in public code/text files, and Python syntax compilation.

## Local targeted validation

A targeted local verification was performed for the tumor-size adjustment analysis using the private local analysis workspace.

The reproduced values were:

- Complete cases: 117
- Size-only R²: 0.067406
- Size + strict RRS R²: 0.147464
- Strict RRS ΔR²: 0.080058
- Strict RRS standardized beta: 0.328
- Strict RRS P value: 0.001484
- Size + signed-log RRS R²: 0.154506
- Signed-log RRS ΔR²: 0.087099
- Signed-log RRS standardized beta: 0.345
- Signed-log RRS P value: 0.000896

These values match the rounded values reported in the manuscript.

## Scope and limitations

Full end-to-end execution was not performed within GitHub Actions because the repository does not redistribute raw imaging data, transcriptomic matrices, clinical records, or patient-level processed files within GitHub Actions because the repository does not redistribute raw imaging data, transcriptomic matrices, clinical. Reproduction of all analyses requires users to obtain the source datasets described in the manuscript and configure local data paths using environment variables.

This validation report should therefore be interpreted as repository-level static validation plus targeted local verification, not as a complete end-to-end rerun of every analysis step.
