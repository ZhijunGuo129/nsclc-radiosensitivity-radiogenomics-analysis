# Code curation record

The release package was prepared from the complete R and Python analysis archives used during manuscript development.

The following materials were retained:

- the final data-acquisition and preprocessing scripts;
- the final CRS and RRS model-development scripts;
- the complete model-building permutation analyses;
- the Lung1 survival, scale-sensitivity, and domain-shift analyses;
- the pathway and Lung3 molecular-replication analyses;
- the final main-figure scripts;
- supplementary table and figure-source scripts;
- environment and repository-validation utilities.

The following materials were excluded:

- R workspace and history files;
- test-only and partial-run scripts;
- superseded figure-layout revisions;
- obsolete low-permutation analyses replaced by the 2,000-permutation workflows;
- manuscript-integration utilities that only copied values into draft tables;
- temporary outputs, downloaded data, derived patient-level files, images, model objects, and spreadsheets;
- local directory paths and workstation-specific settings.

Scripts were renamed by analysis stage and function. Scientific calculations were retained or rewritten against the final analysis specification where an obsolete implementation was identified. In particular, the RRS two-sided permutation sensitivity result is calculated by doubling the smaller empirical tail rather than by applying an absolute-correlation threshold.
