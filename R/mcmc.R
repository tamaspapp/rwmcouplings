
#####
# Hug and Hop
#####

#' Hug + Hop 
#'@export
hughop <- function(x0,Time,Bounces,lam,kap,iter,logpi,gradlogpi){
  
  normalize <- function(x) {return((1/sqrt(sum(x^2)))*x)}
  mu = sqrt(lam * kap);
  
  d = length(x0);
  
  delta = Time / Bounces;
  
  x = x0;
  
  logpi_x = logpi(x);
  
  acc_x_hug = 0; acc_x_hop = 0;
  for(i in 1:iter) {
    # Hug
    # Propose: bounce B times
    z = rnorm(d); # X-velocity
    # Initial position half-step
    xp = x + 0.5 * delta * z;
    g_xp = normalize(gradlogpi(xp));
    
    for(b in 1:Bounces)
    {
      # Reflect velocity in gradient
      z = z -  2 * sum(z * g_xp) * g_xp;
      # Full step
      xp = xp + delta * z;
      g_xp = normalize(gradlogpi(xp));
    }
    # Went too far
    xp = xp - 0.5 * delta * z;
    
    # Accept-reject
    logpi_xp = logpi(xp);
    logHR_x = logpi_xp - logpi_x;
    
    log_u = log(runif(1));
    if (logHR_x > 0 || log_u < logHR_x)
    {
      x = xp;
      logpi_x = logpi_xp;
      acc_x_hug = acc_x_hug + 1;
    }
    
    # Hop
    g_x = gradlogpi(x); gnx = sum(g_x^2);
    
    z = rnorm(d);
    z1 = sum(z * g_x)/gnx *g_x;
    xp <- x + lam/sqrt(gnx) * z1 + mu/sqrt(gnx) *(z - z1)
    
    ###
    # Accept-reject
    ###
    g_xp = gradlogpi(xp); gnxp = sum(g_xp^2);
    logpi_xp = logpi(xp);
    logHR_x = logpi_xp - logpi_x + 0.5*d*log(gnxp/gnx) - 0.5/mu^2*sum((xp-x)^2)*(gnxp - gnx) - 0.5*(1/lam^2 - 1/mu^2)*(sum((xp - x)*g_xp)^2 - sum((x - xp)*g_x)^2);
    
    log_u <- log(runif(1)) # Common uniform acceptance random variable
    if (logHR_x > 0 || log_u < logHR_x) {
      x <- xp
      logpi_x <- logpi_xp
      acc_x_hop <- acc_x_hop + 1
    }
  }
  return(list("x" = x, 
              "acc_hop" = acc_x_hop/iter, 
              "acc_hug" = acc_x_hug/iter))
}

#' CRN Hug + GCRN Hop 
#' Two-scale coupling for Hop: 
#'   1. apply rejection-sampled maximal proposal coupling when fall under some threshold of squared distance
#'  2. apply GCRN when fall above the threshold
#'@export
cplhughop <- function(x0,y0,Time,Bounces,lam,kap,thresh,iter,logpi,gradlogpi) {
  
  normalize <- function(x) {return((1/sqrt(sum(x^2)))*x)}
  loghop <- function(x,y,g_x,gnx) {
    0.5*d*log(gnx) - 0.5/(mu*mu)*sum((y-x)^2)*gnx - 0.5*(1/(lam*lam) - 1/(mu*mu))*sum((y-x)*g_x)^2
  }
  
  mu <- sqrt(lam * kap) 
  
  d = length(x0);
  
  delta = Time / Bounces;
  
  x = x0;
  y = y0;
  
  logpi_x = logpi(x);
  logpi_y = logpi(y);
  
  squaredist = rep(0, iter+1);
  squaredist[1] = sum((x-y)^2);
  
  acc_x_hug = 0; acc_x_hop = 0;
  tau = -1;
  for(i in 1:iter) {
    # Hug
    # Propose: bounce B times
    z = rnorm(d); # X-velocity
    z_y = z;  # Y-velocity
    
    # Initial position half-step
    xp = x + 0.5 * delta * z;
    yp = y + 0.5 * delta * z_y;
    g_xp = normalize(gradlogpi(xp));
    g_yp = normalize(gradlogpi(yp));
    
    for(b in 1:Bounces)
    {
      # Reflect velocity in gradient
      z = z -  2 * sum(z * g_xp) * g_xp;
      z_y = z_y -  2 * sum(z_y * g_yp) * g_yp;
      
      # Full step
      xp = xp + delta * z;
      yp = yp + delta * z_y;
      g_xp = normalize(gradlogpi(xp));
      g_yp = normalize(gradlogpi(yp));
    }
    # Went too far
    xp = xp - 0.5 * delta * z;
    yp = yp - 0.5 * delta * z_y;
    
    # Accept-reject
    logpi_xp = logpi(xp);
    logpi_yp = logpi(yp);
    
    logHR_x = logpi_xp - logpi_x;
    logHR_y = logpi_yp - logpi_y;
    
    log_u = log(runif(1));
    if (logHR_x > 0 || log_u < logHR_x)
    {
      x = xp;
      logpi_x = logpi_xp;
      acc_x_hug = acc_x_hug + 1;
    }
    if (logHR_y > 0 || log_u < logHR_y)
    {
      y = yp;
      logpi_y = logpi_yp;
    }
    
    # Hop
    if(sum((x-y)^2) < thresh) {
      ###
      # Generate proposals
      ###
      g_x = gradlogpi(x); gnx = sum(g_x^2);
      g_y = gradlogpi(y); gny = sum(g_y^2);
      
      z = rnorm(d);
      z1 = sum(z * g_x)/gnx *g_x;
      xp = x + lam/sqrt(gnx)*z1 + mu/sqrt(gnx)*(z - z1);
      g_xp = gradlogpi(xp); gnxp = sum(g_xp^2);
      logpi_xp = logpi(xp);
      
      # Rejection-sample a maximal proposal coupling
      logacc = loghop(y,xp,g_y,gny) - loghop(x,xp,g_x,gnx);
      if(logacc >= 0 || log(runif(1)) < logacc){
        yp = xp;
        logpi_yp = logpi_xp;
        g_yp = g_xp;
        gnyp = gnxp;
        cpl = TRUE;
      }else{
        cpl = FALSE;
        repeat{
          z = rnorm(d);
          z1 = sum(z * g_y)/gny *g_y;
          yp = y + lam/gny *z1 + mu/sqrt(gny) *(z - z1);
          logacc = loghop(y,yp,g_y,gny) - loghop(x,yp,g_x,gnx);
          if(logacc >= 0 || log(runif(1)) < logacc){
            logpi_yp = logpi(yp);
            g_yp = gradlogpi(yp); gnyp = sum(g_yp^2);
            break;
          }
        }
      }
      ###
      # Accept-reject
      ###
      logHR_x = logpi_xp - logpi_x + 0.5*d*log(gnxp/gnx) - 0.5/mu^2*sum((xp-x)^2)*(gnxp - gnx) - 0.5*(1/lam^2 - 1/mu^2)*(sum((xp - x)*g_xp)^2 - sum((x - xp)*g_x)^2);
      logHR_y = logpi_yp - logpi_y + 0.5*d*log(gnyp/gny) - 0.5/mu^2*sum((yp-y)^2)*(gnyp - gny) - 0.5*(1/lam^2 - 1/mu^2)*(sum((yp - y)*g_yp)^2 - sum((y - yp)*g_y)^2);
      
      log_u <- log(runif(1)) # Common uniform acceptance random variable
      if (logHR_x > 0 || log_u < logHR_x) {
        x <- xp
        logpi_x <- logpi_xp
        acc_x_hop <- acc_x_hop + 1
        acc_x = TRUE;
      } else {acc_x = FALSE;}
      if (logHR_y > 0 || log_u < logHR_y) {
        y <- yp
        logpi_y <- logpi_yp
        acc_y = TRUE;
      } else {acc_y = FALSE;}
      
      # Check if we've coupled, and stop
      if(acc_x && acc_y && cpl){tau = i; break;}
      squaredist[i+1] = sum((x-y)^2);
      
    } else { # GCRN
      ###
      # Generate proposals
      ###
      g_x = gradlogpi(x); gnx = sum(g_x^2);
      g_y = gradlogpi(y); gny = sum(g_y^2);
      
      z = rnorm(d);
      z1 = rnorm(1);
      xp <- x + lam/gnx *z1 *g_x + mu/sqrt(gnx) *(z - sum(z * g_x)/gnx *g_x)
      yp <- y + lam/gny *z1 *g_y + mu/sqrt(gny) *(z - sum(z * g_y)/gny *g_y)
      
      ###
      # Accept-reject
      ###
      g_xp = gradlogpi(xp); gnxp = sum(g_xp^2);
      g_yp = gradlogpi(yp); gnyp = sum(g_yp^2);
      
      logpi_xp = logpi(xp);
      logpi_yp = logpi(yp);
      
      logHR_x = logpi_xp - logpi_x + 0.5*d*log(gnxp/gnx) - 0.5* (1/mu^2) * sum((xp-x)^2) * (gnxp - gnx) - 0.5 *  (1/lam^2 - 1/mu^2) * (sum((x - xp)*g_xp)^2 - sum((xp - x)*g_x)^2);
      logHR_y = logpi_yp - logpi_y + 0.5*d*log(gnyp/gny) - 0.5* (1/mu^2) * sum((yp-y)^2) * (gnyp - gny) - 0.5 *  (1/lam^2 - 1/mu^2) * (sum((y - yp)*g_yp)^2 - sum((yp - y)*g_y)^2);
      
      log_u <- log(runif(1)) # Common uniform acceptance random variable
      if (logHR_x > 0 || log_u < logHR_x) {
        x <- xp
        logpi_x <- logpi_xp
        acc_x_hop <- acc_x_hop + 1
      }
      if (logHR_y > 0 || log_u < logHR_y) {
        y <- yp
        logpi_y <- logpi_yp
      }
      
      squaredist[i+1] = sum((x-y)^2);
    }
  }
  
  return(list("squaredist" = squaredist,
              #"acc_x_hug" = acc_x_hug / i,
              #"acc_x_hop" = acc_x_hop / i,
              "tau" = tau))
}


#' CRN Hug + CRN Hop 
#' Two-scale coupling for Hop: 
#'   1. apply rejection-sampled maximal proposal coupling when fall under some threshold of squared distance
#'  2. apply GCRN when fall above the threshold
#'@export
crnhughop <- function(x0,y0,Time,Bounces,lam,kap,thresh,iter,logpi,gradlogpi) {
  
  normalize <- function(x) {return((1/sqrt(sum(x^2)))*x)}
  loghop <- function(x,y,g_x,gnx) {
    0.5*d*log(gnx) - 0.5/(mu*mu)*sum((y-x)^2)*gnx - 0.5*(1/(lam*lam) - 1/(mu*mu))*sum((y-x)*g_x)^2
  }
  
  mu <- sqrt(lam * kap) 
  
  d = length(x0);
  
  delta = Time / Bounces;
  
  x = x0;
  y = y0;
  
  logpi_x = logpi(x);
  logpi_y = logpi(y);
  
  squaredist = rep(0, iter+1);
  squaredist[1] = sum((x-y)^2);
  
  acc_x_hug = 0; acc_x_hop = 0;
  tau = -1;
  for(i in 1:iter) {
    # Hug
    # Propose: bounce B times
    z = rnorm(d); # X-velocity
    z_y = z;  # Y-velocity
    
    # Initial position half-step
    xp = x + 0.5 * delta * z;
    yp = y + 0.5 * delta * z_y;
    g_xp = normalize(gradlogpi(xp));
    g_yp = normalize(gradlogpi(yp));
    
    for(b in 1:Bounces)
    {
      # Reflect velocity in gradient
      z = z -  2 * sum(z * g_xp) * g_xp;
      z_y = z_y -  2 * sum(z_y * g_yp) * g_yp;
      
      # Full step
      xp = xp + delta * z;
      yp = yp + delta * z_y;
      g_xp = normalize(gradlogpi(xp));
      g_yp = normalize(gradlogpi(yp));
    }
    # Went too far
    xp = xp - 0.5 * delta * z;
    yp = yp - 0.5 * delta * z_y;
    
    # Accept-reject
    logpi_xp = logpi(xp);
    logpi_yp = logpi(yp);
    
    logHR_x = logpi_xp - logpi_x;
    logHR_y = logpi_yp - logpi_y;
    
    log_u = log(runif(1));
    if (logHR_x > 0 || log_u < logHR_x)
    {
      x = xp;
      logpi_x = logpi_xp;
      acc_x_hug = acc_x_hug + 1;
    }
    if (logHR_y > 0 || log_u < logHR_y)
    {
      y = yp;
      logpi_y = logpi_yp;
    }
    
    # Hop
    if(sum((x-y)^2) < thresh) {
      ###
      # Generate proposals
      ###
      g_x = gradlogpi(x); gnx = sum(g_x^2);
      g_y = gradlogpi(y); gny = sum(g_y^2);
      
      z = rnorm(d);
      z1 = sum(z * g_x)/gnx *g_x;
      xp = x + lam/sqrt(gnx)*z1 + mu/sqrt(gnx)*(z - z1);
      g_xp = gradlogpi(xp); gnxp = sum(g_xp^2);
      logpi_xp = logpi(xp);
      
      # Rejection-sample a maximal proposal coupling
      logacc = loghop(y,xp,g_y,gny) - loghop(x,xp,g_x,gnx);
      if(logacc >= 0 || log(runif(1)) < logacc){
        yp = xp;
        logpi_yp = logpi_xp;
        g_yp = g_xp;
        gnyp = gnxp;
        cpl = TRUE;
      }else{
        cpl = FALSE;
        repeat{
          z = rnorm(d);
          z1 = sum(z * g_y)/gny *g_y;
          yp = y + lam/gny *z1 + mu/sqrt(gny) *(z - z1);
          logacc = loghop(y,yp,g_y,gny) - loghop(x,yp,g_x,gnx);
          if(logacc >= 0 || log(runif(1)) < logacc){
            logpi_yp = logpi(yp);
            g_yp = gradlogpi(yp); gnyp = sum(g_yp^2);
            break;
          }
        }
      }
      ###
      # Accept-reject
      ###
      logHR_x = logpi_xp - logpi_x + 0.5*d*log(gnxp/gnx) - 0.5/mu^2*sum((xp-x)^2)*(gnxp - gnx) - 0.5*(1/lam^2 - 1/mu^2)*(sum((xp - x)*g_xp)^2 - sum((x - xp)*g_x)^2);
      logHR_y = logpi_yp - logpi_y + 0.5*d*log(gnyp/gny) - 0.5/mu^2*sum((yp-y)^2)*(gnyp - gny) - 0.5*(1/lam^2 - 1/mu^2)*(sum((yp - y)*g_yp)^2 - sum((y - yp)*g_y)^2);
      
      log_u <- log(runif(1)) # Common uniform acceptance random variable
      if (logHR_x > 0 || log_u < logHR_x) {
        x <- xp
        logpi_x <- logpi_xp
        acc_x_hop <- acc_x_hop + 1
        acc_x = TRUE;
      } else {acc_x = FALSE;}
      if (logHR_y > 0 || log_u < logHR_y) {
        y <- yp
        logpi_y <- logpi_yp
        acc_y = TRUE;
      } else {acc_y = FALSE;}
      
      # Check if we've coupled, and stop
      if(acc_x && acc_y && cpl){tau = i; break;}
      squaredist[i+1] = sum((x-y)^2);
      
    } else { # CRN
      ###
      # Generate proposals
      ###
      g_x = gradlogpi(x); gnx = sum(g_x^2);
      g_y = gradlogpi(y); gny = sum(g_y^2);
      
      z = rnorm(d);
      z1x = sum(z * g_x)/gnx*g_x;
      z1y = sum(z * g_y)/gny*g_y;
      xp <- x + lam/sqrt(gnx)*z1x + mu/sqrt(gnx)*(z - z1x)
      yp <- y + lam/sqrt(gny)*z1y + mu/sqrt(gny)*(z - z1y)
      
      ###
      # Accept-reject
      ###
      g_xp = gradlogpi(xp); gnxp = sum(g_xp^2);
      g_yp = gradlogpi(yp); gnyp = sum(g_yp^2);
      
      logpi_xp = logpi(xp);
      logpi_yp = logpi(yp);
      
      logHR_x = logpi_xp - logpi_x + 0.5*d*log(gnxp/gnx) - 0.5* (1/mu^2) * sum((xp-x)^2) * (gnxp - gnx) - 0.5 *  (1/lam^2 - 1/mu^2) * (sum((x - xp)*g_xp)^2 - sum((xp - x)*g_x)^2);
      logHR_y = logpi_yp - logpi_y + 0.5*d*log(gnyp/gny) - 0.5* (1/mu^2) * sum((yp-y)^2) * (gnyp - gny) - 0.5 *  (1/lam^2 - 1/mu^2) * (sum((y - yp)*g_yp)^2 - sum((yp - y)*g_y)^2);
      
      log_u <- log(runif(1)) # Common uniform acceptance random variable
      if (logHR_x > 0 || log_u < logHR_x) {
        x <- xp
        logpi_x <- logpi_xp
        acc_x_hop <- acc_x_hop + 1
      }
      if (logHR_y > 0 || log_u < logHR_y) {
        y <- yp
        logpi_y <- logpi_yp
      }
      
      squaredist[i+1] = sum((x-y)^2);
    }
  }
  
  return(list("squaredist" = squaredist,
              #"acc_x_hug" = acc_x_hug / i,
              #"acc_x_hop" = acc_x_hop / i,
              "tau" = tau))
}


#####
# Coupled samplers, 1 target
#####

#' GCRN coupling
#'@export
RWM_gcrn <- function (x,y,h,iter,logpi,gradlogpi) {
  d <- length(x)
  
  logpi_x <- logpi(x)
  logpi_y <- logpi(y)
  grad_x <- gradlogpi(x)
  grad_y <- gradlogpi(y)
  
  squaredist <- rep(NA, iter + 1)
  squaredist[1] <- sum((x-y)^2)
  
  x_square <- rep(NA, iter + 1)
  x_square[1] <- sum(x^2)
  
  y_square <- rep(NA, iter + 1)
  y_square[1] <- sum(y^2)
  
  accepts <- 0
  
  for(i in 1:iter) {
    ###
    # Proposal noise
    ###
    e_grad_x <- grad_x / sqrt(sum(grad_x^2))
    e_grad_y <- grad_y / sqrt(sum(grad_y^2))
    
    z1 <- rnorm(1)
    z  <- rnorm(d)
    
    xdot <- z + (z1 - sum(z * e_grad_x)) * e_grad_x
    ydot <- z + (z1 - sum(z * e_grad_y)) * e_grad_y
    
    ###
    # Generate proposals
    ###
    xp <- x + h * xdot
    yp <- y + h * ydot
    
    logpi_xp <- logpi(xp)
    logpi_yp <- logpi(yp)
    gradlogpi_xp <- gradlogpi(xp)
    gradlogpi_yp <- gradlogpi(yp)
    
    ###
    # Accept-reject
    ###
    log_u <- log(runif(1)) # Common uniform acceptance random variable
    
    logHR_x <- logpi_xp - logpi_x
    logHR_y <- logpi_yp - logpi_y
    
    if (logHR_x > 0 || log_u < logHR_x) {
      x <- xp
      logpi_x <- logpi_xp
      grad_x <- gradlogpi_xp
      accepts <- accepts + 1
    }
    if (logHR_y > 0 || log_u < logHR_y) {
      y <- yp
      logpi_y <- logpi_yp
      grad_y <- gradlogpi_yp
    }
    
    ###
    # Store
    ###
    x_square[i + 1] <- sum(x^2)
    y_square[i + 1] <- sum(y^2)
    squaredist[i + 1] <- sum((x-y)^2)
  }
  
  return(list("xsq" = x_square,
              "ysq" = y_square,
              "squaredist" = squaredist,
              "acc_rate" = accepts / iter))
}

#' CRN coupling
#'@export
RWM_crn <- function (x,y,h,iter,logpi,gradlogpi = NULL) {
  d <- length(x)
  
  logpi_x <- logpi(x)
  logpi_y <- logpi(y)
  
  squaredist <- rep(NA, iter + 1)
  squaredist[1] <- sum((x-y)^2)
  
  x_square <- rep(NA, iter + 1)
  x_square[1] <- sum(x^2)
  
  y_square <- rep(NA, iter + 1)
  y_square[1] <- sum(y^2)
  
  accepts <- 0
  
  for(i in 1:iter) {
    ###
    # Proposal noise
    ###
    z  <- rnorm(d)
    
    xdot <- z
    ydot <- z
    
    ###
    # Generate proposals
    ###
    xp <- x + h * xdot
    yp <- y + h * ydot
    
    logpi_xp <- logpi(xp)
    logpi_yp <- logpi(yp)
    
    ###
    # Accept-reject
    ###
    log_u <- log(runif(1)) # Common uniform acceptance random variable
    
    logHR_x <- logpi_xp - logpi_x
    logHR_y <- logpi_yp - logpi_y
    
    if (logHR_x > 0 || log_u < logHR_x) {
      x <- xp
      logpi_x <- logpi_xp
      accepts <- accepts + 1
    }
    if (logHR_y > 0 || log_u < logHR_y) {
      y <- yp
      logpi_y <- logpi_yp
    }
    
    ###
    # Store
    ###
    x_square[i + 1] <- sum(x^2)
    y_square[i + 1] <- sum(y^2)
    squaredist[i + 1] <- sum((x-y)^2)
  }
  
  return(list("xsq" = x_square,
              "ysq" = y_square,
              "squaredist" = squaredist,
              "acc_rate" = accepts / iter))
}

#' Reflection coupling
#'@export
RWM_refl <- function (x,y,h,iter,logpi,gradlogpi = NULL) {
  d <- length(x)
  
  logpi_x <- logpi(x)
  logpi_y <- logpi(y)
  
  squaredist <- rep(NA, iter + 1)
  diff <- x - y
  squaredist[1] <- sum(diff^2)
  
  x_square <- rep(NA, iter + 1)
  x_square[1] <- sum(x^2)
  
  y_square <- rep(NA, iter + 1)
  y_square[1] <- sum(y^2)
  
  accepts <- 0
  
  for(i in 1:iter) {
    ###
    # Proposal noise
    ###
    e <- diff / sqrt(squaredist[i])
    
    z  <- rnorm(d)
    
    xdot <- z
    ydot <- z - (2*sum(e*z)) * e 
    
    ###
    # Generate proposals
    ###
    xp <- x + h * xdot
    yp <- y + h * ydot
    
    logpi_xp <- logpi(xp)
    logpi_yp <- logpi(yp)
    
    ###
    # Accept-reject
    ###
    log_u <- log(runif(1)) # Common uniform acceptance random variable
    
    logHR_x <- logpi_xp - logpi_x
    logHR_y <- logpi_yp - logpi_y
    
    if (logHR_x > 0 || log_u < logHR_x) {
      x <- xp
      logpi_x <- logpi_xp
      accepts <- accepts + 1
    }
    if (logHR_y > 0 || log_u < logHR_y) {
      y <- yp
      logpi_y <- logpi_yp
    }
    ###
    # Store
    ###
    x_square[i + 1] <- sum(x^2)
    y_square[i + 1] <- sum(y^2)
    diff <- x - y
    squaredist[i + 1] <- sum(diff^2)
  }
  
  return(list("xsq" = x_square,
              "ysq" = y_square,
              "squaredist" = squaredist,
              "acc_rate" = accepts / iter))
}



#####
# Coupled samplers, 2 targets
#####

#' GCRN coupling
#'@export
RWM2targets_gcrn <- function (x,y,h,iter,logpix,gradlogpix,
                              logpiy,gradlogpiy) {
  d <- length(x)
  
  logpi_x <- logpix(x)
  logpi_y <- logpiy(y)
  grad_x <- gradlogpix(x)
  grad_y <- gradlogpiy(y)
  
  squaredist <- rep(NA, iter + 1)
  squaredist[1] <- sum((x-y)^2)
  
  x_square <- rep(NA, iter + 1)
  x_square[1] <- sum(x^2)
  
  y_square <- rep(NA, iter + 1)
  y_square[1] <- sum(y^2)
  
  accepts <- 0
  
  for(i in 1:iter) {
    ###
    # Proposal noise
    ###
    e_grad_x <- grad_x / sqrt(sum(grad_x^2))
    e_grad_y <- grad_y / sqrt(sum(grad_y^2))
    
    z1 <- rnorm(1)
    z  <- rnorm(d)
    
    xdot <- z + (z1 - sum(z * e_grad_x)) * e_grad_x
    ydot <- z + (z1 - sum(z * e_grad_y)) * e_grad_y
    
    ###
    # Generate proposals
    ###
    xp <- x + h * xdot
    yp <- y + h * ydot
    
    logpi_xp <- logpix(xp)
    logpi_yp <- logpiy(yp)
    gradlogpi_xp <- gradlogpix(xp)
    gradlogpi_yp <- gradlogpiy(yp)
    
    ###
    # Accept-reject
    ###
    log_u <- log(runif(1)) # Common uniform acceptance random variable
    
    logHR_x <- logpi_xp - logpi_x
    logHR_y <- logpi_yp - logpi_y
  
    if (logHR_x > 0 || log_u < logHR_x) {
      x <- xp
      logpi_x <- logpi_xp
      grad_x <- gradlogpi_xp
      accepts <- accepts + 1
    }
    if (logHR_y > 0 || log_u < logHR_y) {
      y <- yp
      logpi_y <- logpi_yp
      grad_y <- gradlogpi_yp
    }
    
    ###
    # Store
    ###
    x_square[i + 1] <- sum(x^2)
    y_square[i + 1] <- sum(y^2)
    squaredist[i + 1] <- sum((x-y)^2)
  }
  
  return(list("xsq" = x_square,
              "ysq" = y_square,
              "squaredist" = squaredist,
              "acc_rate" = accepts / iter))
}
#' CRN coupling
#'@export
RWM2targets_crn <- function (x,y,h,iter,logpix,gradlogpix = NULL,
                             logpiy,gradlogpiy = NULL) {
  d <- length(x)
  
  logpi_x <- logpix(x)
  logpi_y <- logpiy(y)
  
  squaredist <- rep(NA, iter + 1)
  squaredist[1] <- sum((x-y)^2)
  
  x_square <- rep(NA, iter + 1)
  x_square[1] <- sum(x^2)
  
  y_square <- rep(NA, iter + 1)
  y_square[1] <- sum(y^2)
  
  accepts <- 0
  
  for(i in 1:iter) {
    ###
    # Proposal noise
    ###
    z  <- rnorm(d)
    
    xdot <- z
    ydot <- z
    
    ###
    # Generate proposals
    ###
    xp <- x + h * xdot
    yp <- y + h * ydot
    
    logpi_xp <- logpix(xp)
    logpi_yp <- logpiy(yp)
    
    ###
    # Accept-reject
    ###
    log_u <- log(runif(1)) # Common uniform acceptance random variable
    
    logHR_x <- logpi_xp - logpi_x
    logHR_y <- logpi_yp - logpi_y
    
    if (logHR_x > 0 || log_u < logHR_x) {
      x <- xp
      logpi_x <- logpi_xp
      accepts <- accepts + 1
    }
    if (logHR_y > 0 || log_u < logHR_y) {
      y <- yp
      logpi_y <- logpi_yp
    }
    
    ###
    # Store
    ###
    x_square[i + 1] <- sum(x^2)
    y_square[i + 1] <- sum(y^2)
    squaredist[i + 1] <- sum((x-y)^2)
  }
  
  return(list("xsq" = x_square,
              "ysq" = y_square,
              "squaredist" = squaredist,
              "acc_rate" = accepts / iter))
}

#' Reflection coupling
#'@export
RWM2targets_refl <- function (x,y,h,iter,logpix,gradlogpix = NULL,
                              logpiy,gradlogpiy = NULL) {
  d <- length(x)
  
  logpi_x <- logpix(x)
  logpi_y <- logpiy(y)
  
  squaredist <- rep(NA, iter + 1)
  diff <- x - y
  squaredist[1] <- sum(diff^2)
  
  x_square <- rep(NA, iter + 1)
  x_square[1] <- sum(x^2)
  
  y_square <- rep(NA, iter + 1)
  y_square[1] <- sum(y^2)
  
  accepts <- 0
  
  for(i in 1:iter) {
    ###
    # Proposal noise
    ###
    e <- diff / sqrt(squaredist[i])
    
    z  <- rnorm(d)
    
    xdot <- z
    ydot <- z - (2*sum(e*z)) * e 
    
    ###
    # Generate proposals
    ###
    xp <- x + h * xdot
    yp <- y + h * ydot
    
    logpi_xp <- logpix(xp)
    logpi_yp <- logpiy(yp)
    
    ###
    # Accept-reject
    ###
    log_u <- log(runif(1)) # Common uniform acceptance random variable
    
    logHR_x <- logpi_xp - logpi_x
    logHR_y <- logpi_yp - logpi_y
    
    if (logHR_x > 0 || log_u < logHR_x) {
      x <- xp
      logpi_x <- logpi_xp
      accepts <- accepts + 1
    }
    if (logHR_y > 0 || log_u < logHR_y) {
      y <- yp
      logpi_y <- logpi_yp
    }
    ###
    # Store
    ###
    x_square[i + 1] <- sum(x^2)
    y_square[i + 1] <- sum(y^2)
    diff <- x - y
    squaredist[i + 1] <- sum(diff^2)
  }
  
  return(list("xsq" = x_square,
              "ysq" = y_square,
              "squaredist" = squaredist,
              "acc_rate" = accepts / iter))
}










#####
# Target log-density and its gradient
#####

#' Spherical Gaussian target log-density
#'@export
spherical <- function(x) {
  (-0.5) * sum(x^2)
}
#' Spherical Gaussian target log-density gradient
#'@export
spherical_grad <- function(x) {
  -x
}

#' #' Elliptical Gaussian target log-density
#' #'@export
#' elliptical <- function(x, Sigma_inv_chol) {
#'   (-0.5) * sum(as.vector(Sigma_inv_chol %*% x)^2)
#' }
#' #' Elliptical Gaussian target log-density gradient
#' #'@export
#' elliptical_grad <- function(x, Sigma_inv) {
#'   -as.vector(Sigma_inv %*% x)
#' }
