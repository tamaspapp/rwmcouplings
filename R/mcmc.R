
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
    
    log_u <- log(runif(1))
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
#'  1. apply rejection-sampled maximal proposal coupling when fall under some threshold of squared distance
#'  2. apply GCRN when fall above the threshold
#'@export
cplhughop <- function(x0,y0,Time,Bounces,lam,kap,thresh,iter,logpi,gradlogpi) {
  
  mu <- sqrt(lam * kap)
  
  normalize <- function(x) {return((1/sqrt(sum(x^2)))*x)}
  loghop <- function(x,y,g_x,gnx) {
    0.5*d*log(gnx) - 0.5/(mu*mu)*sum((y-x)^2)*gnx - 0.5*(1/(lam*lam) - 1/(mu*mu))*sum((y-x)*g_x)^2
  }
  
  
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


#' RWM + HMC, mixed
#'@export
rwmhmc <- function(x0,
                   Time,L, # HMC params
                   gamma = 0.05, # Probability of RWM
                   h, # RWM params
                   iter,
                   logpi,gradlogpi){
  d = length(x0);
  delta = Time / L;
  
  x = x0;
  logpi_x = logpi(x);
  g_x = gradlogpi(x);
  rwm = F;
  
  acc_x_rwm = acc_x_hmc = 0;
  iter_rwm = iter_hmc = 0;
  for(i in 1:iter) {
    if(runif(1) < gamma) {
      iter_rwm = iter_rwm + 1;

      xp = x + h*rnorm(d);
      logpi_xp = logpi(xp);
      logHR_x  = logpi_xp - logpi_x;
      log_u = log(runif(1));
      if (logHR_x > 0 || log_u < logHR_x)
      {
        x = xp;
        g_x = gradlogpi(xp);
        logpi_x = logpi_xp;
        acc_x_rwm = acc_x_rwm + 1;
      }
      
    } else {
      iter_hmc = iter_hmc + 1;
      
      # Propose: L leapfrog steps
      xp = x;
      g_xp = g_x;
      z = rnorm(d); # X-velocity
      
      zp = z + 0.5 * delta * g_xp # momentum half-step
      for (lfs in 1:L) {
        xp = xp + delta * zp;# position full-step
        
        g_xp = gradlogpi(xp); 
        zp = zp + ifelse(lfs < L, 1, 0.5) * delta * g_xp; # momentum full-step, except the half-step at end
      }
      # zp = -zp # Momentum is symmetric, no need to negate at end
      
      # Accept-reject:
      logpi_xp = logpi(xp);
      logHR_x = logpi_xp - logpi_x - 0.5 * sum(zp^2) + 0.5 * sum(z^2);
      
      log_u = log(runif(1));
      if (logHR_x > 0 || log_u < logHR_x)
      {
        x = xp;
        g_x = g_xp;
        logpi_x = logpi_xp;
        acc_x_hmc = acc_x_hmc + 1;
      }
    }
  }
  return(list("x" = x,
              "acc_x_rwm" = acc_x_rwm / iter_rwm,
              "acc_x_hmc" = acc_x_hmc / iter_hmc))
}

#' two-scale GNRN-RWM + CRN-HMC, mixed
#'@export
twoscalehmc <- function(x0,y0,
                        Time,L, # HMC params
                        gamma = 0.05, # Probability of RWM
                        h,thresh = 0.1, # RWM params
                        iter,
                        logpi,gradlogpi){
  
  normalize <- function(z){z/sqrt(sum(z^2))}
  reflect <- function(x,z){x - (2/sum(z^2)*sum(x*z))*z}
  
  d = length(x0);
  delta = Time / L;
  
  x = x0;
  logpi_x = logpi(x);
  g_x = gradlogpi(x);
  
  y = y0;
  logpi_y = logpi(y);
  g_y = gradlogpi(y);
  
  acc_x_hmc = acc_x_rwm = 0;
  acc_y_hmc = acc_y_rwm = 0;
  iter_rwm = 0;
  iter_hmc = 0;
  
  squaredist = rep(0, iter + 1);
  squaredist[1] = sum((x-y)^2);
  tau = -1;
  rwm = F;
  
  for(i in 1:iter) {
    coupled <- F
    if(runif(1) < gamma) {
      # Two-scale GCRN
      if(squaredist[i] < thresh) { # Refl-max
        iter_rwm = iter_rwm + 1; rwm = T;
        x_dot = rnorm(d)
        xp = x+h*x_dot;
        logpi_xp = logpi(xp);
        
        z = (x - y) / h;
        logcpl = -0.5 * sum((x_dot + z)^2) + 0.5 * sum(x_dot^2);
        
        ucpl = runif(1);
        if (logcpl > 0 || log(ucpl) < logcpl) # Maximal coupling of proposals
        {
          coupled = T;
          yp = xp;
          logpi_yp = logpi_xp;
        }
        else # Reflect
        {
          yp = y+h*reflect(x_dot, z);
          logpi_yp = logpi(yp);
        }
      } else { #GCRN
        e_grad_x <- normalize(g_x)
        e_grad_y <- normalize(g_y)
        
        z1 <- rnorm(1)
        z  <- rnorm(d)
        
        x_dot <- z + (z1 - sum(z * e_grad_x)) * e_grad_x
        y_dot <- z + (z1 - sum(z * e_grad_y)) * e_grad_y
        
        xp <- x + h * x_dot
        yp <- y + h * y_dot
        
        logpi_xp <- logpi(xp)
        logpi_yp <- logpi(yp)
      }
      
      logHR_x <- logpi_xp - logpi_x
      logHR_y <- logpi_yp - logpi_y
      
      g_xp = gradlogpi(xp);
      g_yp = gradlogpi(yp);
      
    } else { # CRN - HMC
      iter_hmc = iter_hmc + 1; rwm = F;
      # Propose: L leapfrog steps
      xp = x;
      g_xp = g_x;
      
      yp = y;
      g_yp = g_y;
      
      z = rnorm(d); # velocity
      
      # momentum half-step
      zpx = z + 0.5 * delta * g_xp 
      zpy = z + 0.5 * delta * g_yp
      for (lfs in 1:L) {
        # position full-step
        xp = xp + delta * zpx;
        yp = yp + delta * zpy;
        
        # momentum full-step, except the half-step at end
        g_xp = gradlogpi(xp);
        g_yp = gradlogpi(yp);
        zpx = zpx + ifelse(lfs < L, 1, 0.5) * delta * g_xp; 
        zpy = zpy + ifelse(lfs < L, 1, 0.5) * delta * g_yp;
      }
      
      # Accept-reject:
      logpi_xp = logpi(xp);
      logHR_x = logpi_xp - logpi_x - 0.5 * sum(zpx^2) + 0.5 * sum(z^2);
      
      logpi_yp = logpi(yp);
      logHR_y = logpi_yp - logpi_y - 0.5 * sum(zpy^2) + 0.5 * sum(z^2);
    }
    
    log_u = log(runif(1));
    if (logHR_x > 0 || log_u < logHR_x)
    {
      x = xp;
      g_x = g_xp;
      logpi_x = logpi_xp;
      if(rwm){acc_x_rwm = acc_x_rwm + 1;} else {acc_x_hmc = acc_x_hmc + 1;}
    }
    if (logHR_y > 0 || log_u < logHR_y)
    {
      y = yp;
      g_y = g_yp;
      logpi_y = logpi_yp;
      if(rwm){acc_y_rwm = acc_y_rwm + 1;} else {acc_y_hmc = acc_y_hmc + 1;}
    }
    
    squaredist[i + 1] = sum((x-y)^2);
    if(isTRUE(all.equal(x,y)) && coupled) {
      tau = i;
      break;
    }
    
  }
  
  return(list(#"x" = x, "y" = y,
              #"acc_x_rwm" = acc_x_rwm / iter_rwm, 
              #"acc_y_rwm" = acc_y_rwm / iter_rwm,
              #"acc_x_hmc" = acc_x_hmc / iter_hmc,
              #"acc_y_hmc" = acc_y_hmc / iter_hmc,
              "tau" = tau,
              "squaredist" = squaredist))
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

#' Optimal coupling
#' #'@export
RWM2targets_optimal <- function(x,y,h,iter,
                                logpix,gradlogpix,
                                logpiy,gradlogpiy) {
  # NB: this is slow
  qemg <- function(p, mu, sigma, lambda) {
    uniroot(function(y){qeng::pemg(y, mu, sigma, lambda) - p}, c(-1e11,1e11), extendInt = "upX", maxiter = 2e9L, check.conv = T, tol = 1e-9)$root
  }
  
  d <- length(x)
  
  logpi_x <- logpix(x)
  logpi_y <- logpiy(y)
  
  grad_x <- gradlogpix(x)
  norm_gx <- sqrt(sum(grad_x^2))
  e_grad_x <- grad_x / norm_gx
  
  grad_y <- gradlogpiy(y)
  norm_gy <- sqrt(sum(grad_y^2))
  e_grad_y <- grad_y / norm_gy
  
  squaredist <- rep(NA, iter + 1)
  squaredist[1] <- sum((x-y)^2)
  
  accepts <- 0
  
  for(i in 1:iter) {
    
    u1 <- runif(1)
    d1 <- emg::qemg(u1,mu = 0, sigma = h * norm_gx, lambda = 1) 
    d2 <- emg::qemg(u1,mu = 0, sigma = h * norm_gy, lambda = 1)
    
    u2 <- runif(1)
    z_gx <- TruncatedNormal::qtnorm(u2, mu = (h * norm_gx)^2, sd = h * norm_gx, lb = -Inf, ub = d1) # will need to divide by (h * norm_gx) to get down to a std Gaussian
    z_gy <- TruncatedNormal::qtnorm(u2, mu = (h * norm_gy)^2, sd = h * norm_gy, lb = -Inf, ub = d2)
    log_ux <- -(d1 - z_gx)
    log_uy <- -(d2 - z_gy)
    z_gx <- z_gx / (h * norm_gx)
    z_gy <- z_gy / (h * norm_gy)
    
    ###
    # Coupled proposals
    ###
    z  <- rnorm(d)
    
    xdot <- z + (z_gx - sum(z * grad_x)/norm_gx)/norm_gx * grad_x
    ydot <- z + (z_gy - sum(z * grad_y)/norm_gy)/norm_gy * grad_y
    
    xp <- x + h * xdot
    yp <- y + h * ydot
    
    logpi_xp <- logpix(xp)
    logpi_yp <- logpiy(yp)
    ###
    # Coupled accept-reject
    ###
    logHR_x <- logpi_xp - logpi_x
    logHR_y <- logpi_yp - logpi_y
    
    if (logHR_x > 0 || log_ux < logHR_x) {
      accepts <- accepts + 1
      
      x <- xp
      logpi_x <- logpi_xp
      
      # Compute gradients only when needed
      grad_x <- gradlogpix(x)
      norm_gx <- sqrt(sum(grad_x^2))
      e_grad_x <- grad_x / norm_gx
    }
    if (logHR_y > 0 || log_uy < logHR_y) {
      y <- yp
      logpi_y <- logpi_yp
      
      grad_y <- gradlogpiy(y)
      norm_gy <- sqrt(sum(grad_y^2))
      e_grad_y <- grad_y / norm_gy
    }
    
    ###
    # Store
    ###
    squaredist[i + 1] <- sum((x-y)^2)
  }
  
  return(list("squaredist" = squaredist,
              "acc_rate" = accepts / iter))
}












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
