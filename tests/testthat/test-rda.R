library(testthat)
library(McMasterPandemic)
library(ggplot2)

test_that("identical_predict", {
    skip_if_not_installed("anytime")
    require(anytime)

    ont_cal_2brks <- fix_stored(ont_cal_2brks)
    ont_cal1 <- fix_stored(ont_cal1)

    new_pred_1 <- predict(ont_cal1)
    new_pred_2brks <- predict(ont_cal_2brks)

    load(system.file("testdata", "ONcalib_2020Jun01.rda", package = "McMasterPandemic"))
    load(system.file("testdata", "ONcalib_2brks_2020Jun01.rda", package = "McMasterPandemic"))

    ## re-hack after loading old stuff ...
    ont_cal_2brks <- fix_stored(ont_cal_2brks)
    ont_cal1 <- fix_stored(ont_cal1)

    old_pred_1 <- predict(ont_cal1)
    old_pred_2brks <- predict(ont_cal_2brks)

    ## hack new_preds a little bit to match old ...
    fix_pred <- function(x) {
        ff <- attr(x, "forecast_args")
        ff$debug_hist <- NULL
        attr(x, "forecast_args") <- ff
        return(x)
    }

    new_pred_1 <- fix_pred(new_pred_1)
    new_pred_2brks <- fix_pred(new_pred_2brks)

    ## waldo::compare(old_pred_1, new_pred_1)
    expect_equal(old_pred_1, new_pred_1)

    ## waldo::compare(old_pred_2brks, new_pred_2brks)
    ## tweak
    attr(attr(old_pred_2brks, "forecast_args")$base_params, "description")[14] <-
        attr(attr(new_pred_2brks, "forecast_args")$base_params, "description")[14]
    expect_equal(new_pred_2brks, old_pred_2brks)
})
