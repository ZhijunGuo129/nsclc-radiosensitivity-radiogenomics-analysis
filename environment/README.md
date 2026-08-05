# Software environments

`install_r_dependencies.R` installs the CRAN and Bioconductor packages required by the R analysis scripts.

`requirements.txt` lists packages used by scripts that run in a standard Python environment. Scripts importing `slicer`, `vtk`, or `DICOMLib` must instead be run in the Python environment bundled with 3D Slicer.

Use `capture_environment.R` and `capture_environment.py` after the final analysis run to record package versions. Environment captures are analysis records and are not generated during repository validation.
