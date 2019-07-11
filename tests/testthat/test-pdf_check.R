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
