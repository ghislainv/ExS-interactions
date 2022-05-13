#!/usr/bin/env Rscript

## ==============================================================================
## author          :James S. Clark
## email           :jimclark@duke.edu
## web             :https://sites.nicholas.duke.edu/clarklab/people/clark/
## license         :GPLv3
## ==============================================================================

# INTRACO project :https://www.fondationbiodiversite.fr/la-frb-en-action/programmes-et-projets/le-cesab/intraco/
# Complexity and diversity explained by the E × S interaction

# TO BE DONE:
# tradeoffs involving C (carrying capacity), z (sensitivity), D (dispersal)
# shifting source of variance ss (noise), vs varBeta/varX
# E (how complex) vs S (no. species) and K (no. sites)
# how does the effect of complexity E depend on local interactions?
# small K, big C means non-local competition

# Load libraries
library("Rcpp")
library("RcppArmadillo")
library("here")

# Source C++ code for function rmvnArma
sourceCpp(here("src/rmvnArma.cpp"))

# Source R functions
source(here("R/R_functions.R"))

# Variables
nsite <- 50      # Number of sites
ninit <- 5       # Initial number of individuals per species
nt <- 100        # Time steps
Q <- 200         # Dimensionality of E (columns in X)
S <- 10 #50          # Number of species
varBeta <- .1    # Determines magnitude of variation in B
varX <- .1       # Determines magnitude of variation in X
#C  <- ninit*10   # Carrying capacity for a local site in survival probit
C <- 5 * 10        # Carrying capacity for a local site in survival probit
z <- C / 2         # Standard deviation in survival probit
D <- .05         # Dispersal fraction
ss <- .001       # Residual variance
tau <- 1e-3      # Nugget variance

# Species parameters
beta <- matrix(rnorm(Q * S, 0, sqrt(varBeta)), Q, S)     # X, B, V
# Get covariance matrix vx for environmental dimensions X
# Must be symmetric, positive semi-definite (add small values on the diagonal)
vx   <- cov(rmvnArma(10, rep(0, Q), varX * diag(Q))) + 1e-10 * diag(Q)
# Draw environmental dimensions for each site
x    <- rmvnArma(nsite, rep(0, Q), vx)
vx   <- cov(x)
cs <- t(beta) %*% vx %*% beta                             # C = B'VB
xb <- x %*% beta                                          # XB
colnames(x) <- rownames(beta) <- paste('q',1:nrow(beta), sep='')

words <- paste('There are ', ninit, ' initial individuals for each of $S$ = ', S, ' species on $K$ = ', nsite, 
               ' sites. Parameter values are $E$ = ', Q, ', $c$ = ', C, ', $q$ = ' , D, 
               ', $z$ = ', z, ', $sigma^2$ = ', ss, ', $tau^2$ = ', tau,
               ', Var[vec($\bB$)] = ', varBeta, '$\bI_{ES}$, Var[vec($\bX$)] = ', varX, '$\bI_{nE}$', sep='')
words

# DESIGN <- 'unstr'   # noisy individuals
# DESIGN <- 'cov'     # spp differences summarized by C
# DESIGN <- 'beta'    # XB (full knowledge)

n <- matrix(ninit, nsite, S)
ntot <- matrix(0, nt, S)

colnames(beta) <- colnames(n) <- colnames(ntot) <- paste('s',1:S, sep='')
rownames(x) <- rownames(n) <- paste('site',1:nsite,sep='_')

nsim <- 2 # 25 # Number of repetitions
nstart <- nt/2

div <- matrix(0, nsim, 3)
colnames(div) <- c('diversity', 'richness', 'sorting') # landscape diversity

# three designs for 5 levels of E
#eseq <- round( c(4, 8, 12, 24, 32, 44, 55, 65, 80, 100, 120, 150) )
#eseq <- round( c(4, 12, 32, 65, 80, 120,  200) )
eseq <- round( c(2, 100) )
eseq <- rev(eseq)                        # must be descending
dseq <- c('unstr', 'cov', 'beta')
design <- as.matrix(expand.grid(eseq, dseq))
colnames(design) <- c('E', 'DESIGN')
nexp <- nrow(design)
effects <- matrix(NA, nexp, 9)
ci <- c( 'Mu', paste(c('-','+'), 'SD', sep='' ))
colnames(effects) <- as.vector( t(outer( colnames(div), ci, paste, sep='' )) )

par(mfcol=c(2, 3), bty='n', mar=c(4, 4, 2, 1), oma=c(4, 3, 2, 4))

for (m in 1:nexp) { # number of experiments
  
  # set.seed(999)
  
  nave  <- ntot * 0 # average over simulations
  betam <- beta
  DES <- design[m, 'DESIGN']
  print(c(m, design[m,]))
  
  for (i in 1:nsim) { # repeat simulation
    
    beta <- matrix( rnorm(Q*S, 0, sqrt(varBeta) ), Q, S )     # X, B, V
    vx   <- cov( rmvnArma(10, rep(0, Q), varX*diag(Q) ) ) + 1e-10 * diag(Q)
    x    <- rmvnArma(nsite, rep(0, Q), vx )
    vx   <- cov(x)
    cs <- t(beta) %*% vx %*% beta                               # C = B'VB
    xb <- x %*% beta                                            # XB
    colnames(x) <- rownames(beta) <- paste('q',1:nrow(beta), sep='')
    
    # GV: This part needs to be checked...
    # if( 'E' %in% colnames(design) ){
    #   q2    <- as.numeric( design[m,'E'] )
    #   
    #   q2 <- 200
    #   
    #   betam <- beta[1:q2,]
    #   xb    <- x[,1:q2]%*%betam
    #   vx    <- cov(x[,1:q2])
    #   cs    <- t(betam)%*%vx%*%betam  # matrix C
    # }
    # 
    print(i)
    
    tmp <- simLoop( DES, ninit, nsite, S, nt, nstart, C, z, cs, ss, tau,
                        D, xb )
    nmean <- tmp$nmean
    ntot  <- tmp$ntot
    
    # Plot species abundances
    if (i==1 & design[m, "E"] %in% c("100", "  2")) {
      plotTime(ntot)
      title(DES)
    }
    
    # sorting: site suitability vs mean abundance
    sorting <- cor( as.vector(xb), as.vector(nmean) )
    
    # species diversity
    tmp        <- diversity( nmean )
    div[i,] <- c( exp(tmp$byTot), tmp$richness, sorting )
    
    nave <- nave + ntot
    
  } # replicate simulations
  
  nave <- nave/nsim     # by-time average over simulations
  
  dnow <-  apply( div, 2, quantile, pnorm(c(0, -1, 1) ), na.rm=T )
  dnow[ is.na(dnow) ] <- 0
  
  effects[m,] <- signif( dnow, 4 )
} # design loop

effects <- data.frame(design, effects)
effects$E <- as.numeric(effects$E)

mtext('Time', side=1, line=.4, outer=TRUE, cex=1.2)
mtext('Abundance', side=2, line=.4, outer=TRUE, cex=1.2)
mtext('E = 100', side=3, line=.4, outer=TRUE, cex=1.2)
mtext('E = 2', side=3, line=-27, outer=TRUE, cex=1.2)

# Save R objects
#save(list=ls(), file=here("outputs", "backup.RData"))
load("outputs/backup.RData")

# ==================================
# Plots
# ==================================

# Plot width
plotwidth <- 12.12537 # Corresponding to default \textwidth for LaTex in cm for a4paper

# Colors
greens <- c('#ffffcc','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32')
dcolor <- c('#d53e4f','#d53e4f','#fc8d59','#fee08b','#ffffbf','#e6f598','#99d594','#3288bd','#3288bd')
qcolor <- c('#7fc97f','#beaed4','#fdc086','#386cb0')
cramp <- colorRampPalette(qcolor)

# =======================================
# Effect of dimensionality on coexistence
cols <- cramp(3)
vars <- c('diversity', 'richness', 'sorting')
# Plot
png(here("outputs", "dimensionality-coexistence.png"), width=plotwidth, height=plotwidth*0.5, res=800, units="cm")
par(bty='n', mfrow=c(1, 3), mar=c(4, 4, 1, 1))
for (k in 1:3) {
  xlab <- ''
  if (k==2) {xlab <- 'E'}
  mmm  <- paste(vars[k], 'Mu', sep='')
  sdds <- paste(vars[k], c('.SD','.SD.1'), sep='')
  plot(NA, xlim = range(effects$E), ylim = range(effects[,grep(vars[k], colnames(effects))]), 
       xlab = xlab, ylab = vars[k])
  d <- unique(effects$DESIGN)
  for(i in 1:length(d)){
    wi <- which(effects$DESIGN == d[i])
    lines(effects$E[wi], effects[wi,mmm], col = cols[i], lwd=2)
    shadeInterval(xvalues= effects$E[wi],
                   loHi = effects[wi, sdds],
                   col = cols[i], add = TRUE, trans = .4)
  }
}
dev.off()

# ==================
# species covariance
par(bty='n', mfrow=c(1, 2))

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

# ============
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

plot( xb, nmean, xlab='Suitability', ylab='Mean abundance', log = 'y' )
points( xb[ ik ], nmean[ ik ], col = 2, pch=16)

mtext('Site', side=2, line=-5, outer  = TRUE, cex=1.2)

# End Of File