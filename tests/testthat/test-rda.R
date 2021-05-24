library(testthat)
library(McMasterPandemic)
library(ggplot2)

test_that('identical_predict', {
skip_if_not_installed('anytime')
require(anytime)

new_pred_1 <- predict(ont_cal1)
new_pred_2brks <- predict(ont_cal_2brks)

load(system.file("testdata", 'ONcalib_2020Jun01.rda', package= "McMasterPandemic"))
load(system.file("testdata", 'ONcalib_2brks_2020Jun01.rda', package= "McMasterPandemic"))

old_pred_1 <- predict(ont_cal1)
old_pred_2brks <- predict(ont_cal_2brks)

pred_1_comparison <- (new_pred_1 ==old_pred_1)
expect(all(pred_1_comparison[!is.na(pred_1_comparison)]), 'New ont_cal1 predictions are not equal to old ont_cal1 predictions!')

pred_2_comparison <- (new_pred_2brks ==old_pred_2brks)
expect(all(pred_2_comparison[!is.na(pred_2_comparison)]), 'New ont_cal_2brks predictions are not equal to old ont_cal_2brks!')
})
