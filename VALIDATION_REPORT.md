# Validation Report

## Repository-level validation

The public repository passed GitHub Actions static validation.

The validation checks included:

- Required repository files and directories
- Repository metadata files
- Absence of machine-specific absolute paths
- Absence of unresolved placeholder metadata
- Absence of CJK characters in public code and text files
- Python syntax compilation

## Scope and limitations

This repository provides the analysis code, configuration files, and documentation associated with the manuscript.

Full end-to-end execution was not performed within GitHub Actions because the repository does not redistribute raw imaging data, transcriptomic matrices, clinical records, or patient-level processed files.

Reproduction of all analyses requires users to obtain the source datasets described in the manuscript and configure local data paths using environment variables.

This validation report should therefore be interpreted as repository-level static validation, not as a complete end-to-end rerun of every analysis step.
