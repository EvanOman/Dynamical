# Creates an orbit diagram for the map x^2 + c + beta/x^2

library(parallel)
library(lattice)
logistic.map <- function(c, x, N, M){
    ## c: bifurcation parameter
    ## x: initial value
    ## N: number of iteration
    ## M: number of iteration points to be returned
    z <- 1:N
    z[1] <- x
    for(i in c(1:(N-1))){
    	z[i+1] <- z[i]^2 + c + .001/z[i]^2	
    }
    
    ## Return the last M iterations 
    z[c((N-M):N)]
}

## Set scanning range for bifurcation parameter r
my.c <- seq(-.21, -.09, length= 20000)
system.time(Orbit <- sapply(my.c, logistic.map,  x=.001^(1/4)
                            , N=500, M=100))


Orbit <- as.vector(Orbit)
c <- sort(rep(my.c, 101))

crit = .001^(.25)

bitmap("pertperdub.tiff", height = 10, width = 15, units = 'in', type="tifflzw", res=500)
xyplot(Orbit ~ c, col = "black", pch = ".", cex = 1, ylim = c(-crit - .1,1), panel=function(...) {
    panel.abline(h=c(crit,0, -crit))
    panel.xyplot(...)
})
dev.off()


