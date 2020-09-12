## An example file for computing credible regions from simulated datapoints
## Requires R version 3.6.3 (Holding the Windsock), and the R package gplots. 
## In this example, the datapoints are in file 'thixibetaBB_allres' and file 'thigenomegBB_allres'
## the output is stored in R objects specifying the points for the credible regions, 
## which can be loaded at a later time and  used for comparison with actual data.
## Here 'credible region' refers to a two-dimensional region containing the specified probability mass of
## of the given datapoints, as estimated by the function "ci2d". 
library(gplots)
m <- matrix(scan("./thixibetaBB_allres"), 31310000, 3, byrow=1)
mgenombetaBBest89 <- ci2d( m[,1:2], show='n', bandwidth=0.01, nbins=100,  ci.levels=0.89)
mgenombetaBBest99 <- ci2d( m[,1:2], show='n', bandwidth=0.01, nbins=100,  ci.levels=0.99)
rm(m)
m <- matrix(scan("./thigenomegBB_allres"), 10540000 , 3, byrow=1)
mgenomegBBest99 <- ci2d( m[,1:2], show='n', bandwidth=0.01, nbins=100,  ci.levels=0.99)
mgenomegBBest89 <- ci2d( m[,1:2], show='n', bandwidth=0.01, nbins=100,  ci.levels=0.89)
rm(m)
save.image(file="thiBBgenomiccredregionobjects.RData")
##save.image(file="thiABgenomiccredregionobjects.RData")
##in another R:
##load("credregionobjects.RData")
quit(save="yes")
