options(stringsAsFactors=FALSE)
library(McMasterPandemic)
load(system.file("testdata","Ontario_basic.rda",package="McMasterPandemic"))
options(stringsAsFactors=FALSE)
test_mle_pred <- predict(Ontario_fit)

all.equal(test_mle_pred, mle_prediction)

str(test_mle_pred)

str(mle_prediction)
