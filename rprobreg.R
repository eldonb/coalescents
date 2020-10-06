# requires the ks library 

gogn <- matrix ( scan("rnorm0"),  990000, 2,  byro=1)

##library(KernSmooth)
#library(MASS)
#library(ks, lib.loc='/home/bjarki/verk/mfn/Rpakkar/')
##library(fields)
d <- data.frame( x=gogn[,1], y=gogn[,2])
rm(gogn)
##k  <-  kde2d(r1, r2, n=250)
##points <- data.matrix( expand.grid( seq(0, 1, by=0.01), seq(0, 1, by=0.01) ))
k <- ks::kde( d, compute.cont=1)


pdf("p0.pdf")
                                        #filled.contour( f1, plot.axes= { contour(f1, levels=0.01, add=TRUE); axis(1); axis(2)} )
## showing contours containing 99% and 89% of distribution
contour( z=k$estimate, levels=c( as.numeric(k$cont[99]), as.numeric(k$cont[89])), xlim=c(0,1), ylim=c(0,1), labcex=.8, lwd=2, las=1, bty='l', labels=c('99%', '89%') )
gogn <- matrix ( scan("rnormK"), 8910000 , 2,  byro=1)
rm(d)
d <- data.frame( x=gogn[,1], y=gogn[,2])
rm(gogn,k)
k <- ks::kde( d, compute.cont=1)
contour( z=k$estimate, levels=c( as.numeric(k$cont[99]), as.numeric(k$cont[89])),  labcex=.8, lwd=2, add=1, lty=2, col='grey', labels=c('99%', '89%') )
rm(k,d)
graphics.off()
quit(save='yes')
