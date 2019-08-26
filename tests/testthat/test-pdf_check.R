## Create some fake pdf files
tmp_dir <- tempdir()
pdf(file.path(tmp_dir, 'four_panels.pdf'))
plot(1, 1)
dev.off()
pdf(file.path(tmp_dir, 'regionCoverage_fractionedData.pdf'))
plot(1, 1)
dev.off()

test_that('pdfs are not overwritten', {
    expect_error(four_panels("chr20:10286777-10288069:+",
        PDF = file.path(tmp_dir, 'four_panels.pdf')))
    expect_error(four_panels("chr20:10286777-10288069:+",
        PDF = file.path(tmp_dir, 'four_panels')))
    expect_error(plot_coverage("chr20:10286777-10288069:+",
        PDF = file.path(tmp_dir, 'regionCoverage_fractionedData.pdf')))
    expect_error(plot_coverage("chr20:10286777-10288069:+",
        PDF = file.path(tmp_dir, 'regionCoverage_fractionedData')))
})


## Choose a random PDF file
PDF <- file.path(tempdir(),
    paste0('test_', stats::runif(1, max = 1e10))
)
PDF

## Initially this works because the output PDF does not exist.
## It also adds the PDF extension if the user didn't supply it.
test_that('check_pdf part own', {
    expect_identical(check_pdf(PDF), paste0(PDF, '.pdf'))
})
PDF <- check_pdf(PDF)

## Create a dummy PDF file
pdf(file = PDF)
plot(1, 1)
dev.off()

## Now it doesn't work since the PDF file already exists.
test_that('check_pdf part two', {
    expect_error(check_pdf(PDF), 'already exists! Rename or erase it')
})
