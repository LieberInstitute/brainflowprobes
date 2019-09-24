# brainflowprobles 0.99.5

NEW FEATURES

* Added a `NEWS.md` file for describing changes across versions.
* Added the `OUTDIR` argument to `check_pdf()` and related functions. Updated
the documentation to reflect this change and also renamed the `PATH` argument
in `region_info()` to `OUTDIR` and changed the default from `'.'` to
`tempdir()` for consistency across functions. This resolves the task suggested
by @Liubuntu at https://github.com/Bioconductor/Contributions/issues/1191#issuecomment-530115005.

# brainflowprobes 0.99.3

NEW FEATURES

* Created three utility functions: `check_pdf()`, `get_nearest_annotation()`
and `get_region_cov()` and resolved other tasks suggested by
@Liubuntu at https://github.com/Bioconductor/Contributions/issues/1191#issuecomment-522762540
during the Bioconductor submission process.

# brainflowprobes 0.99.0

NEW FEATURES

* Created the `four_panels()`, `region_info()` and `plot_coverage()` functions.
* My first R/Bioconductor package! ^^
