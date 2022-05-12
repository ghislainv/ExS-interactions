
# experiments

# tradeoffs involving C (carrying capacity), z (sensitivity), D (dispersal)
# shifting source of variance ss (noise), vs varBeta/varX
# E (how complex) vs S (no. species) and K (no. sites)

# how does the effect of complexity E depend on local interactions?
#   small K, big C means non-local competition

#library('Rcpp')
#library('RcppArmadillo')
#sourceCpp("Rcpp_functions.cpp")
library("mvnfast")

nsite <- 50
ninit <- 5
nt    <- 100
Q <- 200         # dimensionality E (columns in X)
S <- 50          # no. species
varBeta <- .1    # determines magnitude of variation in XB
varX    <- .1

#C  <- ninit*10    # 'carrying capacity for a local site in survival probit

C  <- 5*10    # 'carrying capacity for a local site in survival probit
z  <- C/2         #  standard deviation in survival probit
D  <- .05          # dispersal fraction
ss <- .001          # residual variance
tau <- 1e-3       # nugget variance
beta <- matrix( rnorm(Q*S, 0, sqrt(varBeta) ), Q, S )     # X, B, V
vx   <- cov( rmvn(10, rep(0, Q), varX*diag(Q) ) )
x    <- rmvn( nsite, rep(0, Q), vx )
vx   <- cov(x)
cs <- t(beta)%*%vx%*%beta                                 # C = B'VB
xb <- x%*%beta                                            # XB
colnames(x) <- rownames(beta) <- paste('q',1:nrow(beta), sep='')

words <- paste('There are ', ninit, ' initial individuals for each of $S$ = ', S, ' species on $K$ = ', nsite, 
               ' sites. Parameter values are $E$ = ', Q, ', $c$ = ', C, ', $q$ = ' , D, 
               ', $z$ = ', z, ', $sigma^2$ = ', ss, ', $tau^2$ = ', tau,
               ', Var[vec($\bB$)] = ', varBeta, '$\bI_{ES}$, Var[vec($\bX$)] = ', varX, '$\bI_{nE}$', sep='')
words



par(mfrow=c(2,3), bty='n', mar=c(4,4,1,1), oma=c(4,3,2,4))


DESIGN <- 'beta'    # XB (full knowledge)

DESIGN <- 'cov'     # spp differences summarized by C

DESIGN <- 'unstr'   # noisy individuals

n <- matrix(ninit, nsite, S)
ntot <- matrix(0, nt, S)

colnames(beta) <- colnames(n) <- colnames(ntot) <- paste('s',1:S, sep='')
rownames(x) <- rownames(n) <- paste('site',1:nsite,sep='_')

nsim <- 25
nstart <- nt/2

div <- matrix(0, nsim, 3)
colnames(div) <- c('diversity','richness','sorting') # landscape diversity

# three designs for 5 levels of E

eseq <- round( c(4, 8, 12, 24, 32, 44, 55, 65, 80, 100, 120, 150) )

#eseq <- round( c(4, 12, 32, 65, 80, 120,  200) )


#eseq <- round( c(4, 12) )

eseq <- rev(eseq)                        # must be descending
dseq <- c('unstr','cov','beta')
design <- as.matrix( expand.grid(eseq,dseq) )
colnames(design) <- c('E','DESIGN')

nexp <- nrow(design)

effects <- matrix(NA, nexp, 9)
ci <- c( 'Mu', paste( c('-','+'), 'SD', sep='' ) )
colnames(effects) <- as.vector( t(outer( colnames(div), ci, paste, sep='' )) )

par(mfrow=c(2,4), bty='n', mar=c(4,4,1,1), oma=c(4,3,2,4))

for(m in 1:nexp){           # number of experiments
  
 # set.seed(999)
  
  nave  <- ntot*0            # average over simulations
  betam <- beta
  DES   <- DESIGN
  
  if('DESIGN' %in% colnames(design)){
    DES <- design[m,'DESIGN']
  }
  
  print( c(m, design[m,]) )
  
  for(i in 1:nsim){  # repeat simulation
    
    beta <- matrix( rnorm(Q*S, 0, sqrt(varBeta) ), Q, S )     # X, B, V
    vx   <- cov( rmvn(10, rep(0, Q), varX*diag(Q) ) )
    x    <- rmvn(nsite, rep(0, Q), vx )
    vx   <- cov(x)
    cs <- t(beta)%*%vx%*%beta                                 # C = B'VB
    xb <- x%*%beta                                            # XB
    colnames(x) <- rownames(beta) <- paste('q',1:nrow(beta), sep='')
    
    
    if( 'E' %in% colnames(design) ){
      q2    <- as.numeric( design[m,'E'] )
      
      q2 <- 200
      
      betam <- beta[1:q2,]
      xb    <- x[,1:q2]%*%betam
      vx    <- cov(x[,1:q2])
      cs    <- t(betam)%*%vx%*%betam  # matrix C
    }
    
    print(i)
    
    tmp <- simLoop( DES, ninit, nsite, S, nt, nstart, C, z, cs, ss, tau,
                        D, xb )
    nmean <- tmp$nmean
    ntot  <- tmp$ntot
    
  #  plotTime( ntot )
    
    # sorting: site suitability vs mean abundance
    sorting <- cor( as.vector(xb), as.vector(nmean) )
    
    # species diversity
    tmp        <- diversity( nmean )
    div[i,] <- c( exp(tmp$byTot), tmp$richness, sorting )
    
    nave <- nave + ntot
    
  }       # replicate simlations
  
  nave <- nave/nsim     # by-time average over simulations
  
  dnow <-  apply( div, 2, quantile, pnorm(c(0, -1, 1) ), na.rm=T )
  dnow[ is.na(dnow) ] <- 0
  
  effects[m,] <- signif( dnow, 4 )
}  # design loop

effects <- data.frame(design, effects)
effects$E <- as.numeric(effects$E)


greens <- c('#ffffcc','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32')
dcolor <- c('#d53e4f','#d53e4f','#fc8d59','#fee08b','#ffffbf','#e6f598','#99d594','#3288bd','#3288bd')
qcolor <- c('#7fc97f','#beaed4','#fdc086','#386cb0')

cramp <- colorRampPalette(qcolor)
cols <- cramp(3)


par(bty='n', mfrow=c(1,3), bty='n')

vars <- c('diversity','richness', 'sorting')
xlab <- ''

for(k in 1:3){
  xlab <- ''
  if( k == 2 )xlab <- 'E'
  mmm  <- paste(vars[k], 'Mu',sep='')
  sdds <- paste(vars[k], c('.SD','.SD.1'), sep='')
  plot(NA, xlim = range(effects$E), ylim = range(effects[,grep(vars[k], colnames(effects))]), 
       xlab = xlab, ylab = vars[k])
  d <- unique(effects$DESIGN)
  for(i in 1:length(d)){
    wi <- which(effects$DESIGN == d[i])
    lines(effects$E[wi], effects[wi,mmm], col = cols[i], lwd=2)
    .shadeInterval( xvalues= effects$E[wi],
                    loHi = effects[wi, sdds],
                    col= cols[i], add  = TRUE, trans = .4)
  }
}




title( DESIGN )

mtext('Time', side=1, line=.4, outer  = TRUE, cex=1.2)
mtext('Abundance', side=2, line=.4, outer  = TRUE, cex=1.2)
mtext('E = 100', side=3, line=.4, outer  = TRUE, cex=1.2)
mtext('E = 2', side=3, line=-17, outer  = TRUE, cex=1.2)





# species covariance

par(bty='n', mfrow=c(1,2))

cramp <- colorRampPalette( rev(dcolor) )
cramp <- cramp(40)
cmm  <- max( abs(cs) )
cseq <- seq(-cmm, cmm, length=40)

#cc <- cov2cor(cs)

cols <- findInterval( cs, cseq, all.inside = T)
cols <- cramp[ cols ]
ix   <- expand.grid( 1:S, 1:S )


plot( ix[,1], ix[,2], col = cols, xaxt='n', yaxt='n',
      xlab = '', ylab = '', pch=15, cex=.6, asp=1)
kmax <- apply(cs, 1, which.max )
ik   <- cbind( 1:S, kmax )
points(ik[,2], ik[,1], pch=0, cex=.3, col='brown')

title('E = 2')


# outcomes map
cramp <- colorRampPalette(greens)
cramp <- cramp(40)
cseq <- seq(min(xb), max(xb), length=40)
cols <- findInterval( xb, cseq, all.inside = T)
cols <- cramp[ cols ]
ix   <- expand.grid( 1:S, 1:nsite )

par(bty='n', mfrow=c(1,2))
plot( ix[,1], ix[,2], col = cols, xaxt='n', yaxt='n',
      xlab = 'Species', ylab = '', pch=15, cex=.5, asp=1)
kmax <- apply(xb, 1, which.max )
ik   <- cbind( 1:nsite, kmax )
points(ik[,2], ik[,1], pch=0, cex=.5)

plot( xb, nsum, xlab='Suitability', ylab='Mean abundance', log = 'y' )
points( xb[ ik ], nsum[ ik ], col = 2, pch=16)

mtext('Site', side=2, line=-5, outer  = TRUE, cex=1.2)

##################### functions

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
          gi <- as.vector(rmvn(1, rep(1, nn), cvk + diag( tau, nn ) )  )           # individual
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

diversity <- function( siteBySpec ){ # shannon diversity 
  
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


plotTime <- function( timeBySpec ){
  
  time <- 1:nrow(timeBySpec)
  S    <- ncol(timeBySpec)
  
  plot(time, timeBySpec[,1], type='l', ylim=c(1, 2*ninit*nsite), 
       log = 'y', xlab='', ylab='' )
  for(s in 2:S){
    lines(time, timeBySpec[,s], col=s)
  }
}

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}


.shadeInterval <- function(xvalues,loHi,col='grey',PLOT  = TRUE, add  = TRUE,
                           xlab=' ',ylab=' ', xlim = NULL, ylim = NULL, 
                           LOG = FALSE, trans = .5){
  
  #draw shaded interval
  
  loHi <- as.matrix(loHi)
  
  tmp <- smooth.na(xvalues,loHi)
  
  xvalues <- tmp[,1]
  loHi <- tmp[,-1]
  
  xbound <- c(xvalues,rev(xvalues))
  ybound <- c(loHi[,1],rev(loHi[,2]))
  if(is.null(ylim))ylim <- range(as.numeric(loHi))
  if(is.null(xlim))xlim <- range(xvalues)
  
  if(!add){
    if(!LOG)plot(NULL, xlim = xlim, ylim=ylim, 
                 xlab=xlab, ylab=ylab)
    if(LOG)suppressWarnings( plot(NULL,  xlim = xlim, ylim=ylim, 
                                  xlab=xlab, ylab=ylab, log='y') )
  }
  
  
  if(PLOT)polygon(xbound,ybound, border=NA,col=.getColor(col, trans))
  
  invisible(cbind(xbound,ybound))
  
}

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


.tnorm <- function(n,lo,hi,mu,sig, tiny=0){   
  
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  if(length(lo) != length(mu)){
    print(length(lo))
    print(length(mu))
    stop()
  }
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] + tiny
  z
}


.tnormMVNmatrix <- function(avec, muvec, smat, 
                            lo=matrix(-10000,nrow(muvec),ncol(muvec)), 
                            hi=matrix(10000,nrow(muvec),ncol(muvec)),
                            whichSample = c(1:nrow(smat))){
  
  # lo, hi must be same dimensions as muvec,avec
  # each sample is a row
  
  lo[lo < -10000] <- -10000
  hi[hi > 10000]  <- 10000
  
  if(max(whichSample) > length(muvec))
    stop('\nwhichSample outside length(muvec)\n')
  
  whichSample <- sample(whichSample) # randomize order
  
  nd <- dim(avec)
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, 
                       lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) ) 
  r[,whichSample] <- a[,whichSample]
  r
}
