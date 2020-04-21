## Output examination:
# A list with one element per element in brainflowprobes::pd
# For each dataset, brainflowprobes_cov() returns a list of region
# coverage data.frames. In this example, there was a single input region.
# Then each data.frame itself has 1 row per genome base-pair in the region
# and one column per sample in the dataset

test_that("brainflowprobes_cov", {
    expect_true(is.list(four_panels_example_cov))
    expect_equal(names(four_panels_example_cov), names(brainflowprobes::pd))
    expect_true(all(
        sapply(four_panels_example_cov, length) ==
            length(
                GenomicRanges::GRanges("chr20:10286777-10288069:+")
            )
    ))
    expect_true(all(
        sapply(four_panels_example_cov, function(x) {
            nrow(x[[1]])
        }) ==
            GenomicRanges::width(
                GenomicRanges::GRanges("chr20:10286777-10288069:+")
            )
    ))
    expect_equal(
        sapply(four_panels_example_cov, function(x) {
            ncol(x[[1]])
        }),
        sapply(pd, nrow)
    )
})
