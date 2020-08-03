library(McMasterPandemic)
load(system.file("testdata","Ontario_basic.rda",package="McMasterPandemic"))
test_mle_pred <- predict(Ontario_fit)

str(test_mle_pred)

str(mle_prediction)
