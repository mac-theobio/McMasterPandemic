load("../ontario/mobility1.RData")
ont_cal_mob1
de <- attr(ont_cal_mob1,"de") ## DE component
names(de$member)
P <- de$member$pop
colnames(P) <- colnames(de$member$Sigma)
pairs(P,gap=0)
car::scatterplotMatrix(P
                     , gap = 0
                     , ellipse=TRUE
                     , regLine=FALSE
                     , smooth=FALSE)
library(ggplot2)
library(GGally)
ggpairs(as.data.frame(P))
