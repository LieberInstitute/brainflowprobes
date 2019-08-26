gr <- GenomicRanges::GRanges('chr10:135379301-135379311:+')
ann1 <- get_nearest_annotation(gr)
ann2 <- get_nearest_annotation(gr, CODING_ONLY = TRUE)


test_that("get_nearest_annotation", {
    expect_true(ann2$distance > ann1$distance)
    expect_equal(ann1$distance, 0)
    expect_equal(as.character(ann1$region), 'inside')
    expect_equal(as.character(ann2$region), 'downstream')
})
