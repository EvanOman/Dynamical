points <- seq(-2,2, length=20000)
p <- points
sing <- function(x, c, beta){return(x^2+c+beta/(x^2))}
c <- -.062456
bet <- .001

png(filename="./bounded.png",  width=1000, height=30000)
par(mfrow=c(30,1))
for (i in 1:30)
{
    p <- sing(p, c, bet)
    plotable <- points[p < 5]
    plot(x=plotable, y=rep(0, length(plotable)), pch = ".", xlab="point", ylab="")
}

dev.off()

