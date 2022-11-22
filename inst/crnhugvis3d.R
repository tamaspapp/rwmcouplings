##############
# Visualisation of contractive behaviour of Hug kernels under a CRN coupling
##############
# (Mostly by hand.)
# Hacked to shreds... but it works.

pdf(file = "crnhugvis.pdf",  width = 8, height = 8) 
par(cex = 2,pty="s",mar = rep(0, 4))

##' Orthographic projection
##' @param r Latitude-longitude coordinates in a matrix with columns
##' labelled \code{phi} (latitude) and \code{lambda} (longitude)
##' @param proj.centre Location of centre of projection as matrix with
##' column names \code{phi} (elevation) and \code{lambda} (longitude).
##' @param ... Arguments not used by this projection.n
##' @return Two-column matrix with columns labelled \code{x} and
##' \code{y} of locations of projection of coordinates on plane
##' @references \url{http://en.wikipedia.org/wiki/Map_projection},
##' \url{http://mathworld.wolfram.com/OrthographicProjection.html}
##' @author David Sterratt
##' @export
orthographic <- function(r,
                         proj.centre=cbind(phi=0, lambda=0),
                         ...) {
  if (is.character(r) & identical(r, "boundary")) {
    return(circle(360))
  }
  lambda0 <- proj.centre[1, "lambda"]
  phi0    <- proj.centre[1, "phi"]
  
  ## First translate to Cartesian coordinates with only the rotation
  ## about the polar axis so that (0, lambda0) is at the centre of the
  ## projection.
  P <- cbind(x=cos(r[,"phi"])*sin(r[,"lambda"] - lambda0),
             y=sin(r[,"phi"]),
             z=cos(r[,"phi"])*cos(r[,"lambda"] - lambda0))
  
  ## The rotate about the x axis so that (phi0, lambda0) is at the
  ## centre
  P <- P %*% rbind(c(1, 0, 0),
                   c(0,  cos(phi0), sin(phi0)),
                   c(0, -sin(phi0), cos(phi0)))
  
  ## Projection onto the x-y plane
  rc <- cbind(x=P[, 1], y=P[, 2])
  
  ## Anything with negative z is obsured
  rc[P[, 3] < 0,] <- NA
  return(rc)
}

# We'll need radians for the orthographic projection
deg2rad <- function(deg) {(deg * pi) / (180)}
rad2deg <- function(rad) {(rad * 180) / (pi)}

library(geosphere) # For plotting great circles

projectedGreatCircleInRadians <- function(p1,p2,proj.centre, n = 1e3) {
  proj.centre <- deg2rad(proj.centre)
  gc <- deg2rad(geosphere::greatCircle(p1,p2,n))
  gc <- data.frame(gc); names(gc) <- c("lambda", "phi")
  out <- orthographic(gc,proj.centre)
  return(out)
}

projectedGreatCircleIntermediateInRadians <- function(p1,p2,proj.centre, n = 1e3) {
  proj.centre <- deg2rad(proj.centre)
  gc <- deg2rad(geosphere::gcIntermediate(p1,p2,n, addStartEnd=T))
  gc <- data.frame(gc); names(gc) <- c("lambda", "phi")
  out <- orthographic(gc,proj.centre)
  return(out)
}

projectedPointInRadians <- function(p,proj.centre) { 
  p <- deg2rad(p)
  proj.centre <- deg2rad(proj.centre)
  out <- orthographic(p,proj.centre)
  return(out)
}

# We'll work with degrees longitude and latitude

# (longitude first, latitude second)
proj.centre <- cbind(lambda=0, phi = 15)
x0_deg <- cbind(lambda = -45, phi = -10) 
y0_deg <- cbind(lambda = 45, phi = -10) 
A_deg <- cbind(lambda = 0, phi = 80) 

# Define x1 and y1 as intersections of appropriate grate circles. Makes plotting a segment easier
point1.x1y1gc <- cbind(lambda = -45, phi = 40)
point2.x1y1gc <- cbind(lambda = 45, phi = 40)
x1_deg <- gcIntersect(x0_deg,A_deg,point1.x1y1gc,point2.x1y1gc)[1,1:2]
x1_deg <- cbind(lambda = x1_deg[1], phi = x1_deg[2])
y1_deg <- gcIntersect(y0_deg,A_deg,point1.x1y1gc,point2.x1y1gc)[1,1:2]
y1_deg <- cbind(lambda = y1_deg[1], phi = y1_deg[2])

# Project down
x0 <- projectedPointInRadians(x0_deg, proj.centre)
y0 <- projectedPointInRadians(y0_deg, proj.centre)
A <- projectedPointInRadians(A_deg, proj.centre)
x1 <- projectedPointInRadians(x1_deg, proj.centre)
y1 <- projectedPointInRadians(y1_deg, proj.centre)

gc3 <- na.omit(projectedGreatCircleInRadians(x0_deg,A_deg,proj.centre))
gc1 <- na.omit(projectedGreatCircleInRadians(y0_deg,A_deg,proj.centre))
gc2 <- na.omit(projectedGreatCircleIntermediateInRadians(x0_deg,y0_deg,proj.centre))


#Draw Sphere outline
plot.limits<- c(-1,1)
plot(sin(2*pi*seq(0,1,0.01)), cos(2*pi*seq(0,1,0.01)), xlim = plot.limits, ylim = plot.limits, col = "blue", type = "l", 
     main = NULL, xlab = NA, ylab = NA, axes=F, frame.plot=F)

#Draw great circles behind the spheere
lines(-gc1, lty = 2)
lines(-gc3, lty = 2)

#Draw great circles (and great circle segments) on surface

#Draw line x0 -> y0 in red
lines(gc2, lty = 6)

#For coloring, separate out segment between x0->x1
gc3.1 <- gc3[gc3[,2] <= x0[2],] #segment below x0
gc3.2 <- gc3[gc3[,2] >= x1[2],] #segment above x1
gc3.3 <- na.omit(projectedGreatCircleIntermediateInRadians(x0_deg,x1_deg,proj.centre)) #segment between x0 and x1
lines(gc3.1)
lines(gc3.2)
lines(gc3.3, col = "red")

#For coloring, separate out segment between y0->y1
gc1.1 <- gc1[gc1[,2] <= y0[2],] 
gc1.2 <- gc1[gc1[,2] >= y1[2],]
gc1.3 <- na.omit(projectedGreatCircleIntermediateInRadians(y0_deg,y1_deg,proj.centre))#gc1[(gc1[,2] >= y0[2]) & (gc1[,2] <= y1[2]),]
lines(gc1.1)
lines(gc1.2)
lines(gc1.3, col = "red")

#Draw line x1 -> y1 in red
gc4<-na.omit(projectedGreatCircleIntermediateInRadians(x1_deg,y1_deg,proj.centre))
lines(gc4, lty = 6, col = "red")

#Draw points x0,y0,A
points(x0, pch = 19)
text(x0, expression(X[t]), pos = 2)
points(y0, pch = 19)
text(y0, expression(Y[t]), pos = 4)
points(A, pch = 19)
text(A, "A", pos = 1)

#Draw points x1,y1
points(y1, pch = 19, col = "red")
text(y1, expression(Y[t+1]), pos = 4, col = "red")
points(x1, pch = 19, col = "red")
text(x1, expression(X[t+1]), pos = 2, col = "red")

#Draw velocity at x0 and y0
v_offset <- c(0,0.35)
v0.x <- rbind(x0,x0+v_offset)
v0.y <- rbind(y0,y0+v_offset)

lines(v0.x, col = "red", lwd = 3)
points(x0+v_offset,col="red", pch = 17)
text(x0+v_offset, expression(v), pos = 2, col = "red")
lines(v0.y, col = "red", lwd = 3)
points(y0+v_offset,col="red", pch = 17)
text(y0+v_offset, expression(v), pos = 4, col = "red")

dev.off()
