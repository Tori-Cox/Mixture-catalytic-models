#######################################
## functions for fitting the mixture ##
#######################################



### calculate weibull distribution parameters

    mweibull <- function(mean, stdev)
      # Returns the Method-of-Moments shape and scale parameter estimates for the Weibull distribution, given a mean and variance.
      # Ref:  Garcia, O.  NZ Journal of Forestry Science
    {
      cv <- stdev / mean
      # if (cv > 1.2 || cv <= 0)
      #  stop("Out of range")
      a <- 1 / (cv * (((((0.007454537 * cv * cv - 0.08354348) * cv + 0.153109251) *
                          cv - 0.001946641) * cv - 0.22000991) * (1 - cv) * (1 - cv) + 1))
      list(shape = a, scale = mean / gamma(1 + 1 / a))
    }


   
### fit distributions
    
    # N.B. the following code is adapted from the mixdist package
    # Ref: Macdonald P, Du J. mixdist: Finite Mixture Distribution. 2018
    
    grp_prob <- function (mixdat, mixpar, dist1, dist2) { 
      m <- nrow(mixdat)
      mu <- mixpar[, 2]
      sigma <- mixpar[, 3]
      Gdist1 <- as.character(dist1)
      Gdist2 <- as.character(dist2)
      
      
      # distribution 1
      if(Gdist1 == "norm"){
        par1 <- mu[1]
        par2 <- sigma[1]
        mixcdf <- t(sapply(mixdat[-m, 1], pnorm, par1, par2))
      }else if (Gdist1 == "gamma"){
        par1 <- (mu[1]/sigma[1])^2
        par2 <- mu[1]/(sigma[1]^2)
        mixcdf <- t(sapply(mixdat[-m, 1], pgamma, par1, rate = par2))
      }else if (Gdist1 == "weibull") {
        par <- mweibull(mu[1], sigma[1])
        par1 <- par$shape
        par2 <- par$scale
        mixcdf <- t(sapply(mixdat[-m, 1], pweibull, shape=par1, scale=par2))
      }
      
      
      # distribution 2
      if(Gdist2 == "norm"){
        par1 <- mu[2]
        par2 <- sigma[2]
        mixcdf2 <- t(sapply(mixdat[-m, 1], pnorm, par1, par2))
      }else if (Gdist2 == "gamma"){
        par1 <- (mu[2]/sigma[2])^2
        par2 <- mu[2]/(sigma[2]^2)
        mixcdf2 <- t(sapply(mixdat[-m, 1], pgamma, par1, rate = par2))
      }else if (Gdist2 == "weibull") {
        par <- mweibull(mu[2], sigma[2])
        par1 <- par$shape
        par2 <- par$scale
        mixcdf2 <- t(sapply(mixdat[-m, 1], pweibull, shape=par1, scale=par2))
      }
      
      mixcdf <- cbind(t(mixcdf),t(mixcdf2))
      mixcdf <- rbind(mixcdf, 1) - rbind(0, mixcdf)
    }


    mix <- function(mixdat, mixpar, dist1, dist2, constr,
                emsteps = emsteps, 
                exptol = 5e-06, ...){
  
                # check distributions are typed in correct format
                if (is.na(match(dist1, c("norm", "gamma", "weibull")))) 
                  stop(paste("Unknown / unsupported distribution 1: ", dist1, ". Try norm, gamma or weibull", sep = ""))
                if (is.na(match(dist2, c("norm", "gamma", "weibull")))) 
                  stop(paste("Unknown / unsupported distribution 2: ", dist2, ". Try norm, gamma or weibull", sep = ""))
  
   
    estep <- function(emixdat, emixpar, edist1, edist2, econstr){
      
      prj.i <- grp_prob(emixdat, emixpar, edist1, edist2)
      prji <- sweep(prj.i, 2, emixpar[, 1], "*")
      pri.j <- sweep(prji, 1, apply(prji, 1, sum), "/")
      pseudodat <- sweep(pri.j, 1, emixdat[, 2], "*")
      nplusi <- apply(pseudodat, 2, sum)
      pi <- nplusi/sum(nplusi)
      
      list(pi = pi, pseudomixdat = data.frame(x = emixdat[, 1], pseudodat = pseudodat))
    }
  
  
  
    mstep <- function(mspseudomixdat, msmixpar, msdist1, msdist2, msconstr, ...) {
    
    
    
    multineg2llg <- function(mulpseudomixdat, mulmixpar, 
                             muldist1, muldist2, mulconstr, multheta){
      
      muly <- as.matrix(mulpseudomixdat[, -1])
      muln <- apply(muly, 2, sum)
      mulpara <- mixtheta2par(multheta, mulmixpar, mulconstr,mixprop = FALSE)
      upmixpar <- data.frame(mulpi = mulmixpar[, 1], mulmu = mulpara[, 2], mulsigma = mulpara[, 3])
      
      mulparvalid <- testpar(upmixpar, muldist1, mulconstr)
      mulparvalid2 <- testpar(upmixpar, muldist2, mulconstr)
      
      
      ## to ensure weibull params don't exceed limits
      if(muldist1=="weibull" & upmixpar$mulsigma[1]/upmixpar$mulmu[1]>1.2){
        mulparvalid <- FALSE
      }
      if(muldist2=="weibull" & upmixpar$mulsigma[2]/upmixpar$mulmu[2]>1.2){
        mulparvalid2 <- FALSE
      }
      
      
      if (mulparvalid & mulparvalid2) {
        mulp <- grp_prob(mulpseudomixdat, upmixpar, muldist1, muldist2)
        mulnp <- sweep(mulp, 2, muln, "*")
        mulyu <- as.vector(muly)
        mulnpu <- as.vector(mulnp)[mulyu > 0]
        mulyu <- mulyu[mulyu > 0]
        min(-2 * sum(mulyu * log(mulnpu/mulyu)), 1e+16)
      }
      else 1e+16
    }
    
    
    
    theta0 <- mixpar2theta(msmixpar, msconstr, mixprop = FALSE) 
    
    msnlmobj <- nlm(multineg2llg, mulpseudomixdat = mspseudomixdat, 
                    mulmixpar = msmixpar, muldist1 = msdist1, muldist2 = msdist2, mulconstr = msconstr, 
                    p = theta0)#, ...)
    
      
    mspara <- mixtheta2par(msnlmobj$est, msmixpar, msconstr,mixprop = FALSE)
    msmu <- mspara[, 2]
    mssigma <- mspara[, 3]
    
    data.frame(pi = msmixpar[, 1], mu = msmu, sigma = mssigma)
    
  }
  
  
  fitpar <- list()
  fitpar[[1]] <- mixpar
  for (j in 1:emsteps) {
      
    fitpi <- estep(emixdat=mixdat, emixpar=fitpar[[j]], edist1=dist1, edist2=dist2, econstr=constr)
    fitpar[[j]][, 1] <- fitpi$pi
    
    fitpar[[j+1]] <- mstep(mspseudomixdat = fitpi$pseudomixdat, msmixpar=fitpar[[j]], msdist1=dist1, msdist2=dist2, msconstr=constr)
   
  }
  fitpar <- fitpar[[emsteps+1]]
  print(c(dist1, dist2))
  print(c("Params from EM step", fitpar))
  

  mixtheta0 <- mixpar2theta(fitpar, constr)
  
  mixlike <- function(lmixdat, lmixpar, ldist1, ldist2, lconstr, lmixtheta, lexptol) {
    
    lk <- nrow(lmixpar)
    lm <- nrow(lmixdat)
    uppar <- mixtheta2par(lmixtheta, lmixpar, lconstr)
    parvalid <- testpar(uppar, ldist1, lconstr)
    if (!parvalid) {
      like <- 1e+16
      ldf <- lm - 1
    }
    parvalid <- testpar(uppar, ldist2, lconstr)
    if (!parvalid) {
      like <- 1e+16
      ldf <- lm - 1
    }
    else {
      pmat <- grp_prob(lmixdat, uppar, ldist1, ldist2)
      ln <- sum(lmixdat[, 2])
      joint <- sweep(pmat, 2, uppar[, 1], "*")
      mixed <- apply(joint, 1, sum)
      ly <- as.vector(as.matrix(lmixdat[, 2]))
      
      pu <- as.vector(mixed)[ly > 0]
      npu <- ln * pu
      yu <- ly[ly > 0]
      
        like <- min(-2 * sum(yu * log(npu/yu)), 1e+16)
        
        ldf <- sum(mixed > lexptol) - 1 - length(lmixtheta)
      
    }
    attr(like, "df") <- ldf
    like
  }
  
  nlmobj <- nlm(mixlike, lmixdat = mixdat, lmixpar = fitpar, 
                ldist1 = dist1, ldist2 = dist2, lconstr = constr, 
                lexptol = exptol, p = mixtheta0, hessian = TRUE)
  
  param <- mixtheta2par(nlmobj$est, mixpar, constr)
 
  
  ## for CIs on mu and sd
          covmat <- function(covmixpar, covconstr, covhessian) {
            invmat <- try(solve(covhessian/2))
            if (inherits(invmat, "try-error")) 
              invmat <- matrix(NA, nrow = nrow(covhessian), ncol = ncol(covhessian))
            dr <- nrow(invmat)
            covk <- nrow(covmixpar)
            if (covconstr$conpi == "NONE") 
              lpi <- covk - 1
            else if (covconstr$conpi == "PFX" & sum(covconstr$fixpi) < 
                     covk - 1) {
              pi.e <- covmixpar[!covconstr$fixpi, 1]
              lpi <- length(pi.e) - 1
            }
            else lpi <- 0
            if (covconstr$conmu == "NONE") 
              mu.e <- covmixpar[, 2]
            lmu <- length(mu.e)
            if (covconstr$consigma == "NONE") 
              sigma.e <- covmixpar[, 3]
            lsigma <- length(sigma.e)
            sigmat <- diag(c(rep(1, dr - lsigma), sigma.e))
            vmat <- sigmat %*% invmat %*% sigmat
            se <- sqrt(diag(vmat))
            pi.se <- rep(NA, covk)
            if (lpi > 0) {
              pi.se[1:covk] <- c(se[1:lpi], sqrt(sum(vmat[1:lpi, 
                                                          1:lpi])))
            }
            mu.se <- rep(NA, covk)
            if (lmu > 0) {
              mu.se[1:lmu] <- se[(lpi + 1):(lpi + lmu)]
            }
            sigma.se <- rep(NA, covk)
            if (lsigma > 0) {
              sigma.se[1:lsigma] <- se[(lpi + lmu + 1):(lpi + 
                                                          lmu + lsigma)]
            }
            list(vmat = vmat, mixse = data.frame(pi.se = pi.se, mu.se = mu.se, 
                                                 sigma.se = sigma.se))
          }
          
          vmat <- covmat(param, constr, nlmobj$hessian)
  
  
  
  chisq <- nlmobj$minimum
  df <- attr(mixlike(mixdat, fitpar, dist1, dist2, constr, nlmobj$est, 
                     exptol), "df")
  P <- 1 - pchisq(chisq, df)
  mixfit <- list(parameters = param, distribution = list(dist1, dist2), se = vmat$mixse, 
                 chisq = chisq, P = P,df=df, vmat = vmat$vmat, 
                 mixdata = mixdat)
  
  print(c("Params from second stage", mixfit$parameters))
  return(mixfit)
  
}




