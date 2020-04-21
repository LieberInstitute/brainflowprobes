test_that("get_region_cov", {
    expect_equal(
        get_region_cov(COVERAGE = four_panels_example_cov),
        four_panels_example_cov
    )
    ## Order of the names doesn't matter as long as all the data is there
    expect_equal(
        get_region_cov(
            COVERAGE = four_panels_example_cov[
                rev(names(brainflowprobes::pd))
            ]
        ),
        four_panels_example_cov
    )
    expect_error(
        get_region_cov(COVERAGE = list("hola" = 1))
    )
    expect_error(
        get_region_cov(
            COVERAGE = list(
                "Sep" = matrix(),
                "Deg" = matrix(),
                "Cell" = matrix(),
                "Sort" = matrix()
            )
        ),
        "region coverage data.frames"
    )
    expect_error(
        get_region_cov(
            COVERAGE = list(
                "Sep" = list("1" = matrix()),
                "Deg" = list("1" = matrix()),
                "Cell" = list("1" = matrix()),
                "Sort" = list("1" = matrix())
            )
        ),
        "region coverage data.frames"
    )
    expect_error(
        get_region_cov(
            COVERAGE = list(
                "Sep" = list("1" = data.frame()),
                "Deg" = list("1" = data.frame()),
                "Cell" = list("1" = data.frame()),
                "Sort" = list("1" = data.frame())
            )
        ),
        "should have a column per each of the samples in"
    )
    expect_error(
        get_region_cov(
            COVERAGE = list(
                "Sep" = list("1" = data.frame(matrix(ncol = 23, nrow = 2))),
                "Deg" = list("1" = data.frame(matrix(ncol = 40))),
                "Cell" = list("1" = data.frame(matrix(ncol = 466))),
                "Sort" = list("1" = data.frame(matrix(ncol = 12)))
            )
        ),
        "should have the same number of rows"
    )
})
