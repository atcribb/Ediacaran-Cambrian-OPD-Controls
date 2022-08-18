#Authors: Alison Cribb, Sebastiaan van de Velde
#Code to simulate the relative flux of O2 sourced from bioirrigation vs diffusion
#need Model_Function.R loaded in to run analyses 

# #====== PACKAGES ========#
library(ReacTran)
library(rootSolve)
library(marelac)
library(writexl)

#====== BIOTURBATION PARAMETERS FOR EACH TIME INTERVAL ======#
bioturbation.parms <- as.data.frame(matrix(NA, nrow=2, ncol=3))
colnames(bioturbation.parms) <- c('Ediacaran', 'Terreneuvian', 'Modern')
rownames(bioturbation.parms) <- c('Db', 'irr')
bioturbation.parms['Db',] <- c(0.01, 0.98, 10)
bioturbation.parms['irr',] <- c(3.65, 35.77, 365)

#====== FUNCTIONS ======#
# function to calculate tortuiosity:
f.tort <- function(por) 1-2*log(por)
tort   <- function(D.grid, por.grid){
  D.grid$mid <- D.grid$mid/f.tort(por.grid$mid)
  D.grid$int <- D.grid$int/f.tort(por.grid$int)
  return(D.grid)
}

# make sure model function from Model_Function is in environment!

#===========================================================
# Parameters list
#===========================================================
PL <- list() #create empty parameters list 

# species model keeps track of 
PL$N.var <- 9
PL$var.names <- c("CH2O.f","CH2O.s","O2","SO4","FeS","HCO3","NH4","NO3","HS")

# reaction rates
PL$N.rate <- 6
PL$rate.names <- c("flux.up","flux.down","reac","irr","ddt","deficit")
PL$rate.summary <- data.frame(matrix(nrow=PL$N.var,ncol=PL$N.rate),row.names=PL$var.names)
names(PL$rate.summary) <- PL$rate.names

# 7 reactions the model simulates 
PL$N.reac <- 7
PL$reac.names <- c("R1","R2","R3","R4","R5","R6","R7")
PL$reaction.summary <- data.frame (matrix(nrow=PL$N.reac,ncol=1),row.names=PL$reac.names)
names(PL$reaction.summary) <- c("rate")

#======== MODEL DOMAIN AND GRID DEFINITION =========#
PL$L <- 15   # depth of sediment domain (cm)
PL$N <- 400  # number of grid layers
PL$grid <- setup.grid.1D(x.up=0, x.down=PL$L, N=PL$N, dx.1=PL$L/2000, p.dx.1=1.1)
PL$Depth <- PL$grid$x.mid

#======= PARAMETERS =========#
PL$S       <- 30       # salinity
PL$TC      <- 25       # temperature (deg C)
PL$P       <- 1.013    # pressure (bar)
PL$pH      <- 7.5

#porosity profile:
PL$por.0     <- 0.8    # porosity at the sediment-water interface
PL$por.inf   <- 0.8    # asymptotic porosity at depth
PL$por.x.att <- 1.0    # attenuation depth (cm)

PL$por.grid <- setup.prop.1D(func=p.exp, grid=PL$grid, y.0=PL$por.0,     y.inf=PL$por.inf,     x.L=0, x.att=PL$por.x.att)
PL$svf.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=(1-PL$por.0),y.inf=(1-PL$por.inf),x.L=0,x.att=PL$por.x.att)

#transport parameters:
PL$rho.sed <- 2.6     # density solid sediment (g cm-3)
PL$v.0     <- 0.2     # sedimentation velocity (cm yr-1)
PL$v.inf   <- 0.2     # sedimentation velocity (cm yr-1)
PL$velocity.info <- setup.compaction.1D(v.0 = PL$v.0, por.0=PL$por.0, por.inf=PL$por.inf, por.grid=PL$por.grid)
PL$v.grid <- PL$velocity.info$v
PL$u.grid <- PL$velocity.info$u

#Molecular diffusion (O2, SO4, HCO3, HS)
c.fac <- 10000*(3600*24*365.25) # unit conversion

Dmol.all <- c.fac*diffcoeff(S=PL$S, t=PL$TC, P=PL$P, species=c('O2', 'SO4', 'HCO3', 'NH4', 'NO3', 'HS'))

PL$D.O2.grid <- setup.prop.1D(value=Dmol.all$O2, grid=PL$grid)
PL$D.O2.grid <- tort(PL$D.O2.grid, PL$por.grid)

PL$D.SO4.grid <- setup.prop.1D(value=Dmol.all$SO4, grid=PL$grid)
PL$D.SO4.grid <- tort(PL$D.SO4.grid, PL$por.grid)

PL$D.HCO3.grid <- setup.prop.1D(value=Dmol.all$HCO3, grid=PL$grid)
PL$D.HCO3.grid <- tort(PL$D.HCO3.grid, PL$por.grid)

PL$D.NH4.grid <- setup.prop.1D(value=Dmol.all$NH4, grid=PL$grid)
PL$D.NH4.grid <- tort(PL$D.NH4.grid, PL$por.grid)

PL$D.NO3.grid <- setup.prop.1D(value=Dmol.all$NO3, grid=PL$grid)
PL$D.NO3.grid <- tort(PL$D.NO3.grid, PL$por.grid)

PL$D.HS.grid <- setup.prop.1D(value=Dmol.all$HS, grid=PL$grid)
PL$D.HS.grid <- tort(PL$D.HS.grid, PL$por.grid)

#========= REACTION PARAMETERS =========#
#organic matter decay
PL$k.f <- 10.0     # decay constant organic matter (year-1) - Katsev et al. 2006
PL$k.s <- 0.1      # decay constant organic matter (year-1) - Fossing et al. 2004

PL$K_O2  <- 0.001  # Monod constant O2 consumption (umol cm-3 or mM) - van Cappellen and Wang (1996), Meysman et al. (2015)
PL$K_NO3 <- 0.001  # Monod constant NO3 reduction  (umol cm-3 or mM) - van Cappellen and Wang (1996)
PL$K_SO4 <- 0.9    # Monod constant SO4 reduction  (umol cm-3 or mM) - Meysman et al. (2015)

PL$k_Sox <- 1E+06  # kinetic constant sulfide oxidation (umol-1 cm3 yr-1) - Meysman et al. (2015)
PL$k_NH4 <- 1E+06  # kinetic constant NH4 oxidation (umol-1 cm3 yr-1)  - van Cappellen and Wang (1996)
PL$k_Sni <- 1E+06  # kinetic sulfide oxidation with nitrate (umol-1 cm3 yr-1) - Meysman et al. (2015)

PL$CNratio <- 106.0/16.0 # C to N ratio organic matter
PL$f.FeS   <- 0.1        # fraction of sulphide precipitating as FeS  

#==========================================#
#======= LOW END MEMBER SIMULATIONS =======#
#==========================================#

#======= BOUNDARY CONDITIONS ========#
#flux boundary conditions:
CH2O.tot     <- c(150, 300, 450)
frac.CH2O    <- 0.5 #ratio of labile:refractory organic matter 
F.CH2O.f.seq <- frac.CH2O*CH2O.tot        # organic matter deposition (umol cm-2 yr-1)
F.CH2O.s.seq <- (1-frac.CH2O)*CH2O.tot    # organic matter deposition (umol cm-2 yr-1)
PL$F.FeS     <- 0.0                       # FeS deposition (umol cm-2 yr-1)

#modern concentrations for normalization
O2.ow.max   <- 0.28  #O2 concentration bottom water [umol cm-3 or mM]
SO4.ow.max  <- 28.8  # SO4 concentration bottom water [umol cm-3 or mM]
NO3.ow.max  <- 0.1 # NO3 concentration bottom water [umol cm-3 or mM]
HCO3.ow.max <- 2.2   # HCO3 concentration bottom water [umol cm-3 or mM]
NH4.ow.max  <- 0.0   # NH4 concentration bottom water [umol cm-3 or mM]
HS.ow.max   <- 0.0   # HS concentration bottom water [umol cm-3 or mM]

#upper boundary concentrations - 25% PAL atmO2
PL$O2.ow   <- 0.25*O2.ow.max   # O2 concentration bottom water (25% PAL)     (umol cm-3 or mM)
PL$SO4.ow  <- 0.25*SO4.ow.max  # SO4 concentration bottom water  (umol cm-3 or mM)
PL$HCO3.ow <- 2.2              # HCO3 concentration bottom water (umol cm-3 or mM)
PL$NH4.ow  <- 0.0              # NH4 concentration bottom water  (umol cm-3 or mM)
PL$HS.ow   <- 0.0              # Fe concentration bottom water   (umol cm-3 or mM)
PL$NO3.ow  <- 0.25*NO3.ow.max   # NO3 concentration bottom water  (umol cm-3 or mM)

#==================================#
#============MODEL RUN=============#
#==================================#
#Set up orgniac carbon flux sequences
CH2O.fast.parms <- c(rep(F.CH2O.f.seq[1],(ncol(bioturbation.parms)+1)), rep(F.CH2O.f.seq[2],(ncol(bioturbation.parms)+1)), rep(F.CH2O.f.seq[3],(ncol(bioturbation.parms)+1)))
CH2O.slow.parms <- c(rep(F.CH2O.s.seq[1],(ncol(bioturbation.parms)+1)), rep(F.CH2O.s.seq[2],(ncol(bioturbation.parms)+1)), rep(F.CH2O.s.seq[3],(ncol(bioturbation.parms)+1)))
CH2O.tot.parms  <- CH2O.fast.parms + CH2O.slow.parms

#Set up bioturbation sequences
Db.0.parms <- rep(c(0, bioturbation.parms['Db','Ediacaran'], bioturbation.parms['Db','Terreneuvian'], bioturbation.parms['Db','Modern']),length(CH2O.tot))
irr.0.parms <-rep(c(0, bioturbation.parms['irr','Ediacaran'], bioturbation.parms['irr','Terreneuvian'], bioturbation.parms['irr','Modern']),length(CH2O.tot))
bioturbation_labels <- rep(c('None', 'Ediacaran', 'Terreneuvian', 'Modern'),length(CH2O.tot))

#set up dataframe to save results
results_df <- as.data.frame(matrix(data=NA, nrow=length(CH2O.tot.parms), ncol=7))
colnames(results_df) <- c('O2', 'F.CH2O', 'Bioturbation', 'Db', 'irr', 'irr.O2', 'diff.O2')
results_df$O2 <- PL$O2.ow

#run model
for(i in 1:length(CH2O.tot.parms)){
  
  #=====CH2O and bioturbation paramaters=====#
  #organic matter flux
  PL$F.CH2O.f <- CH2O.fast.parms[i]
  PL$F.CH2O.s <- CH2O.slow.parms[i]
  
  #biodiffusion transport
  PL$Db.0     <- Db.0.parms[i]    # biodiffusion coefficient Db at SWI [cm2 yr-1]
  PL$L.mix    <- 1 + 9*(1-exp((-PL$Db.0/3))) # Mixed depth layer
  PL$Db.inf   <- 0.0      # deep Db [cm2 yr-1]
  PL$Db.x.att <- 2      # attenuation depth [cm]
  PL$Db.grid  <- setup.prop.1D(func=p.sig, grid=PL$grid,y.0=PL$Db.0,y.inf=PL$Db.inf,x.L=PL$L.mix,x.att=PL$Db.x.att)
  
  #bioirrigation transport
  PL$L.irr     <- 0.0  # irrigation depth [cm]
  PL$irr.0     <- irr.0.parms[i]  # irrigation rate at SWI [yr-1]
  PL$irr.inf   <- 0.0  # deep irrigation rate [yr-1]
  PL$irr.x.att <- 3.0  # irrigation attenuation coef [cm]
  PL$irr.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$L.irr,x.att=PL$irr.x.att)
  
  #====== Initial conditions ========#
  #Initialize state variable vectors
  CH2O.f  <- 0
  CH2O.s  <- 0
  O2    <- PL$O2.ow
  SO4   <- PL$SO4.ow
  FeS   <- 0
  HCO3  <- PL$HCO3.ow
  NH4   <- PL$NH4.ow
  NO3   <- PL$NO3.ow
  HS    <- PL$HS.ow
  
  N <- PL$N
  yini <- vector(length=PL$N.var*PL$N)
  yini[(0*N+1):(1*N)] <- CH2O.f
  yini[(1*N+1):(2*N)] <- CH2O.s
  yini[(2*N+1):(3*N)] <- O2
  yini[(3*N+1):(4*N)] <- SO4
  yini[(4*N+1):(5*N)] <- FeS
  yini[(5*N+1):(6*N)] <- HCO3
  yini[(6*N+1):(7*N)] <- NH4
  yini[(7*N+1):(8*N)] <- NO3
  
  #====== Run model to steady state =====#
  out <- steady.1D(y=yini, func=model, parms=PL, nspec=PL$N.var, pos=TRUE)
  print(paste('Completed run', i, 'of', length(CH2O.tot.parms)))
  
  #====== save output in results_df ======#
  results_df$F.CH2O[i]          <- PL$F.CH2O.f + PL$F.CH2O.s
  results_df$Bioturbation[i]    <- bioturbation_labels[i]
  results_df$Db[i]              <- PL$Db.0
  results_df$irr[i]             <- PL$irr.0
  results_df$irr.O2[i]          <- out$int.irr.O2
  results_df$diff.O2[i]         <- out$O2.diff.flux.up
  
}

lowO2_irr_results <- results_df
save(lowO2_irr_results, file='O2irrflux_lowO2_081022.RData')

#==========================================#
#======= HIGH END MEMBER SIMULATIONS =======#
#==========================================#

#======= BOUNDARY CONDITIONS ========#
#flux boundary conditions:
CH2O.tot     <- c(150, 300, 450)
frac.CH2O    <- 0.5 #ratio of labile:refractory organic matter 
F.CH2O.f.seq <- frac.CH2O*CH2O.tot        # organic matter deposition (umol cm-2 yr-1)
F.CH2O.s.seq <- (1-frac.CH2O)*CH2O.tot    # organic matter deposition (umol cm-2 yr-1)
PL$F.FeS     <- 0.0                       # FeS deposition (umol cm-2 yr-1)

#modern concentrations for normalization
O2.ow.max   <- 0.28  #O2 concentration bottom water [umol cm-3 or mM]
SO4.ow.max  <- 28.8  # SO4 concentration bottom water [umol cm-3 or mM]
NO3.ow.max  <- 0.1 # NO3 concentration bottom water [umol cm-3 or mM]
HCO3.ow.max <- 2.2   # HCO3 concentration bottom water [umol cm-3 or mM]
NH4.ow.max  <- 0.0   # NH4 concentration bottom water [umol cm-3 or mM]
HS.ow.max   <- 0.0   # HS concentration bottom water [umol cm-3 or mM]

#upper boundary concentrations - 25% PAL atmO2
PL$O2.ow   <- 0.50*O2.ow.max   # O2 concentration bottom water (25% PAL)     (umol cm-3 or mM)
PL$SO4.ow  <- 0.50*SO4.ow.max  # SO4 concentration bottom water  (umol cm-3 or mM)
PL$HCO3.ow <- 2.2              # HCO3 concentration bottom water (umol cm-3 or mM)
PL$NH4.ow  <- 0.0              # NH4 concentration bottom water  (umol cm-3 or mM)
PL$HS.ow   <- 0.0              # Fe concentration bottom water   (umol cm-3 or mM)
PL$NO3.ow  <- 0.25*NO3.ow.max   # NO3 concentration bottom water  (umol cm-3 or mM)

#==================================#
#============MODEL RUN=============#
#==================================#
#Set up orgniac carbon flux sequences
CH2O.fast.parms <- c(rep(F.CH2O.f.seq[1],(ncol(bioturbation.parms)+1)), rep(F.CH2O.f.seq[2],(ncol(bioturbation.parms)+1)), rep(F.CH2O.f.seq[3],(ncol(bioturbation.parms)+1)))
CH2O.slow.parms <- c(rep(F.CH2O.s.seq[1],(ncol(bioturbation.parms)+1)), rep(F.CH2O.s.seq[2],(ncol(bioturbation.parms)+1)), rep(F.CH2O.s.seq[3],(ncol(bioturbation.parms)+1)))
CH2O.tot.parms  <- CH2O.fast.parms + CH2O.slow.parms

#Set up bioturbation sequences
Db.0.parms <- rep(c(0, bioturbation.parms['Db','Ediacaran'], bioturbation.parms['Db','Terreneuvian'], bioturbation.parms['Db','Modern']),length(CH2O.tot))
irr.0.parms <-rep(c(0, bioturbation.parms['irr','Ediacaran'], bioturbation.parms['irr','Terreneuvian'], bioturbation.parms['irr','Modern']),length(CH2O.tot))
bioturbation_labels <- rep(c('None', 'Ediacaran', 'Terreneuvian', 'Modern'),length(CH2O.tot))

#set up dataframe to save results
results_df <- as.data.frame(matrix(data=NA, nrow=length(CH2O.tot.parms), ncol=7))
colnames(results_df) <- c('O2', 'F.CH2O', 'Bioturbation', 'Db', 'irr', 'irr.O2', 'diff.O2')
results_df$O2 <- PL$O2.ow

#run model
for(i in 1:length(CH2O.tot.parms)){
  
  #=====CH2O and bioturbation paramaters=====#
  #organic matter flux
  PL$F.CH2O.f <- CH2O.fast.parms[i]
  PL$F.CH2O.s <- CH2O.slow.parms[i]
  
  #biodiffusion transport
  PL$Db.0     <- Db.0.parms[i]    # biodiffusion coefficient Db at SWI [cm2 yr-1]
  PL$L.mix    <- 1 + 9*(1-exp((-PL$Db.0/3))) # Mixed depth layer
  PL$Db.inf   <- 0.0      # deep Db [cm2 yr-1]
  PL$Db.x.att <- 2      # attenuation depth [cm]
  PL$Db.grid  <- setup.prop.1D(func=p.sig, grid=PL$grid,y.0=PL$Db.0,y.inf=PL$Db.inf,x.L=PL$L.mix,x.att=PL$Db.x.att)
  
  #bioirrigation transport
  PL$L.irr     <- 0.0  # irrigation depth [cm]
  PL$irr.0     <- irr.0.parms[i]  # irrigation rate at SWI [yr-1]
  PL$irr.inf   <- 0.0  # deep irrigation rate [yr-1]
  PL$irr.x.att <- 3.0  # irrigation attenuation coef [cm]
  PL$irr.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$L.irr,x.att=PL$irr.x.att)
  
  #====== Initial conditions ========#
  #Initialize state variable vectors
  CH2O.f  <- 0
  CH2O.s  <- 0
  O2    <- PL$O2.ow
  SO4   <- PL$SO4.ow
  FeS   <- 0
  HCO3  <- PL$HCO3.ow
  NH4   <- PL$NH4.ow
  NO3   <- PL$NO3.ow
  HS    <- PL$HS.ow
  
  N <- PL$N
  yini <- vector(length=PL$N.var*PL$N)
  yini[(0*N+1):(1*N)] <- CH2O.f
  yini[(1*N+1):(2*N)] <- CH2O.s
  yini[(2*N+1):(3*N)] <- O2
  yini[(3*N+1):(4*N)] <- SO4
  yini[(4*N+1):(5*N)] <- FeS
  yini[(5*N+1):(6*N)] <- HCO3
  yini[(6*N+1):(7*N)] <- NH4
  yini[(7*N+1):(8*N)] <- NO3
  
  #====== Run model to steady state =====#
  out <- steady.1D(y=yini, func=model, parms=PL, nspec=PL$N.var, pos=TRUE)
  print(paste('Completed run', i, 'of', length(CH2O.tot.parms)))
  
  #====== save output in results_df ======#
  results_df$F.CH2O[i]          <- PL$F.CH2O.f + PL$F.CH2O.s
  results_df$Bioturbation[i]    <- bioturbation_labels[i]
  results_df$Db[i]              <- PL$Db.0
  results_df$irr[i]             <- PL$irr.0
  results_df$irr.O2[i]          <- out$int.irr.O2
  results_df$diff.O2[i]         <- out$O2.diff.flux.up
  
}

highO2_irr_results <- results_df
save(highO2_irr_results, file='O2irrflux_highO2_081022.RData')


