#!/usr/bin/env Rscript

## ==============================================================================
## author          :James S. Clark
## email           :jimclark@duke.edu
## web             :https://sites.nicholas.duke.edu/clarklab/people/clark/
## license         :GPLv3
## ==============================================================================

# INTRACO project :https://www.fondationbiodiversite.fr/la-frb-en-action/programmes-et-projets/le-cesab/intraco/
# Complexity and diversity explained by the E × S interaction

# Loop for simulations
simLoop <- function(DES, ninit, nsite, S, nt, nstart, C, z, cs, ss, tau,
                    D, xb){
  
  n <- matrix(ninit, nsite, S)
  disp <- rep(0, S)
  nsum <- n*0         # average over nt - nstart time steps
  
  slab <- rep( 1:S, each = ninit )
  n0   <- length(slab)
  slabel <- matrix( 0, nsite, ninit*C )
  rownames(slabel) <- rownames(x)
  slabel[,1:n0] <- matrix( slab, nsite, n0, byrow=T)
  glabel <- slabel*0
  
  for(t in 1:nt){  # time iteration
    
    ns <- rowSums(n)         # no. individuals by site
    ntot[t,] <- colSums(n)   # time by species
    disp <- disp*0           # global dispersal pool
    
    for(k in 1:nsite){  # all sites
      
      nk <- rep(0, ns[k])          # individuals at k
      if(length(nk) == 0)next
      nn <- length(nk)
      mk <- rep(1:S, n[k,])
      names(nk) <- mk              # species label for nn individuals
      ik <- as.matrix( expand.grid(mk,mk) )
      
      cvk <- matrix( cs[ ik ], nn, nn )  # individual by indiviaual
      c0  <- sum( diag(cvk) )/nn                 # total variance divided by population size
      
      if(DES == 'cov'){
      
        if(nn > 1){
          gi <- as.vector(rmvnArma(1, rep(1, nn), cvk + diag( tau, nn ) )  )           # individual
        }else{
          gi <- rnorm(1, 1, sqrt(cvk[1] + tau) )
        }
        
        # inherit RE from previous individual
        now2last <- match( mk, slabel[k,] ) # align current spec label to previous
        labk <- glabel[k,]
        wf   <- which( is.finite(now2last) )
        if(t == 1){
          labk[1:nn] <- gi
        }
        gi[ wf ] <- labk[ now2last[wf] ]
        
        glabel[k,1:nn] <- gi
        slabel[k,1:nn] <- mk
        
        if( nn < ncol(slabel) ){
          nc <- (nn+1):ncol(slabel)
          slabel[k,nc] <- glabel[k,nc] <- NA
        }
      }
      if( DES == 'beta' ){
        gi <- rnorm( nn, 1 + xb[k,mk], sqrt(tau) )
      }
      if( DES == 'unstr' ){
        gi <- rnorm( nn, 1, sqrt( c0 + tau ) )
      }
      ge <- rnorm( nn, 0, sqrt( ss ) )
      g <- gi + ge
      g[ g < 0 ] <- 0
      
      ts   <- tapply( g, mk, sum )   # sum individuals
      it   <- as.numeric(names(ts))
      
      #  survival
      ps <-  1 - pnorm( (sum(ts) - C)/z )  
      if(ps < .07)ps <- .07
      
      n[k,it]  <-  (1 - D)*ts*ps
      disp[it] <- disp[it] + D*ts*ps
      n[ is.na(n) | n < 0 ] <- 0
    }
    
    di <- matrix(disp/nsite, nsite, S, byrow=T)
    n  <- round( n +  di )
    
    if(t > nstart) nsum <- nsum + n
    
  } # time loop
  
  nmean <- nsum/(nt - nstart)
  list( nmean = nmean, ntot = ntot, nlast = n )
}

# Shannon diversity
diversity <- function( siteBySpec ) {
  
  pmat   <- sweep( siteBySpec, 1, rowSums(siteBySpec), '/' )
  pmat[ pmat == 0 ] <- NA
  bySite <- -rowSums( pmat*log(pmat), na.rm=T )
  
  byTot  <- colSums( round(siteBySpec) )     # at least one individual
  rich   <- length( which(byTot > 1) )
  byTot  <- byTot/sum(byTot)
  bt     <- byTot[ byTot > 0 ]
  byTot  <- -sum( bt*log(bt) )
  if( is.na(byTot) )byTot <- 0
  
  list( bySite = bySite, byTot = byTot, richness = rich )
}

# plotTime
plotTime <- function( timeBySpec ){
  
  time <- 1:nrow(timeBySpec)
  S    <- ncol(timeBySpec)
  
  plot(time, timeBySpec[,1], type='l', ylim=c(1, 2*ninit*nsite), 
       log = 'y', xlab='', ylab='' )
  for(s in 2:S){
    lines(time, timeBySpec[,s], col=s)
  }
}

# getColor
.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}

# smooth.na
smooth.na <- function(x,y){   
  
  #remove missing values
  #x is the index
  #y is a matrix with rows indexed by x
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  
  wy <- which(!is.finite(y),arr.ind   = TRUE)
  if(length(wy) == 0)return(cbind(x,y))
  wy <- unique(wy[,1])
  ynew <- y[-wy,]
  xnew <- x[-wy]
  
  return(cbind(xnew,ynew))
}

# shadeInterval
shadeInterval <- function(xvalues,loHi,col='grey',PLOT  = TRUE, add  = TRUE,
                          xlab=' ',ylab=' ', xlim = NULL, ylim = NULL, 
                          LOG = FALSE, trans = .5){
  
  #draw shaded interval
  
  loHi <- as.matrix(loHi)
  
  tmp <- smooth.na(xvalues,loHi)
  
  xvalues <- tmp[,1]
  loHi <- tmp[,-1]
  
  xbound <- c(xvalues,rev(xvalues))
  ybound <- c(loHi[,1],rev(loHi[,2]))
  if (is.null(ylim)) ylim <- range(as.numeric(loHi))
  if (is.null(xlim)) xlim <- range(xvalues)
  
  if (!add) {
    if (!LOG) plot(NULL, xlim = xlim, ylim=ylim, 
                   xlab=xlab, ylab=ylab)
    if (LOG) suppressWarnings(plot(NULL, xlim = xlim, ylim=ylim, 
                              xlab=xlab, ylab=ylab, log='y'))
  }
  
  if (PLOT) polygon(xbound, ybound, border=NA, col=.getColor(col, trans))
  invisible(cbind(xbound, ybound))
  
}

# End Of File