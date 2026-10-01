test_that('soyface biomass data tables all have the same column names', {
    cnames <- colnames(soyface_biomass[[1]])

    for (i in seq(2, length(soyface_biomass))) {
        expect_equal(
            cnames,
            colnames(soyface_biomass[[i]])
        )
    }
})
