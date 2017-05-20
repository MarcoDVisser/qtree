##' quadtree
##'
##' 
##'
##'
##'
##'
##'
##'
##'
##'
##'
quadtree <- function(xy, k=1) {
  d <- dim(xy)[2]
  quad <- function(xy, i, id=1) {
    if (length(xy) < 2*k*d) {
      rv = list(id=id, value=xy)
      class(rv) <- "quadtree.leaf"
    }
    else {
      q0 <- (1 + runif(1,min=-1/2,max=1/2)/dim(xy)[1])/2 # Random quantile near the median
      x0 <- quantile(xy[,i], q0)
      ## Modulus: 1 becomes 2, 2 becomes 1 when d=2
      j <- i %% d + 1 # (Works for octrees, too...) 
      rv <- list(index=i, threshold=x0, range=range(xy[,i]),
                 lower=quad(xy[xy[,i] <= x0, ], j, id*2), 
                 upper=quad(xy[xy[,i] > x0, ], j, id*2+1))
      class(rv) <- "quadtree"
    }
    return(rv)
  }
  quad(xy, 1)
}

points.quadtree <- function(q, ...) {
  points(q$lower, ...); points(q$upper, ...)
}

points.quadtree.leaf <- function(q, ...) {
    points(q$value, col=hcl(q$id), ...)
}


lines.quadtree <- function(q, xylim, ...) {
  i <- q$index
  j <- 3 - q$index
  clip <- function(xylim.clip, i, upper) {
    if (upper) xylim.clip[1, i] <- max(q$threshold, xylim.clip[1,i]) else 
      xylim.clip[2,i] <- min(q$threshold, xylim.clip[2,i])
    xylim.clip
  } 
  if(q$threshold > xylim[1,i]) lines(q$lower, clip(xylim, i, FALSE), ...)
  if(q$threshold < xylim[2,i]) lines(q$upper, clip(xylim, i, TRUE), ...)
  xlim <- xylim[, j]
  xy <- cbind(c(q$threshold, q$threshold), xlim)
  lines(xy[, order(i:j)],  ...)
}

lines.quadtree.leaf <- function(q, xylim, ...) {} # Nothing to do at leaves!

##'extract records
##' @param q previously build quadtree
##' @param x extaction coordinates in the first dimension
##' @param y extaction coordinates in the first dimension
##' @param r extraction radius

### SUBSET based on range
qextract <- function(q,x,y,r=1){

    refxy <- c(x,y)
    uppxy <- refxy+r
    lowxy <- refxy-r
    
    climbTree <- function(q,xy, i) {

        d <- dim(xy)
        ## Modulus: 1 becomes 2, 2 becomes 1 when d=2
        j <- i %% d + 1 # (Works for octrees, too...)

        if(class(q)=="quadtree.leaf"){
        return(q$value)
        } else {
            
            ## check overlap in ranges
            ## (StartA <= EndB)  and  (EndA >= StartB)
            if((lowxy[j]<=q$lower$range[2])&(uppxy[j]>=q$lower$range[1])){
                climbTree(q$lower,xy,j)
#                cat("$lower")
            } else if((lowxy[j]<=q$upper$range[2])&(uppxy[j]>=q$upper$range[1])) {
              climbTree(q$upper,xy,j)
#                 cat("$upper")
            }
        }
        }

  
     climbTree(q,xy=refxy,1)

}



n <- 2500       # Points per cluster
n.centers <- 20  # Number of cluster centers
sd <- 1/2        # Standard deviation of each cluster
set.seed(17)
centers <- matrix(runif(n.centers*2, min=c(-90, 30), max=c(-75, 40)), ncol=2, byrow=TRUE)
xy <- matrix(apply(centers, 1, function(x) rnorm(n*2, mean=x, sd=sd)), ncol=2, byrow=TRUE)
k <- 5
system.time(qt <- quadtree(xy, k))
#
# Set up to map the full extent of the quadtree.
#
xylim <- cbind(x=c(min(xy[,1]), max(xy[,1])), y=c(min(xy[,2]), max(xy[,2])))
plot(xylim, type="n", xlab="x", ylab="y", main="Quadtree")
#
# This is all the code needed for the plot!
#
points(qt, pch=16)
lines(qt, xylim, col="gray")
