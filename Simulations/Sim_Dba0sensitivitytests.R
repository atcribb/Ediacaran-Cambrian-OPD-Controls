#Authors: Alison Cribb, Sebastiaan van de Velde
#Code to simulate biodiffusion and bioirrigation coefficient sensitivity analyses in Figure 1 and Figures S1, respectively
#need Model_Function.R loaded in to run anlayses 

#====== PACKAGES ========#
library(ReacTran)
library(rootSolve)
library(marelac)
library(writexl)


#==== SENSITIVITY ANALYSIS FOR BIODIFFUSION COEFFICIENTS =====#
#Generate sensitivity analysis curve by increasing Db from 0 to 3

#OM influx level (low    = 150 umol cm-2 yr-1
#                 medium = 300 umol cm-2 yr-1 
#                 high   = 700 umol cm-2 yr-1 )

CH2O.tot <- 150 #low
#CH2O.tot <- 300 #medium
#CH2O.tot <- 700 #high

#====== FUNCTIONS ======#
# function to calculate tortuiosity:
f.tort <- function(por) 1-2*log(por)
tort   <- function(D.grid, por.grid){
  D.grid$mid <- D.grid$mid/f.tort(por.grid$mid)
  D.grid$int <- D.grid$int/f.tort(por.grid$int)
  return(D.grid)
}

#==== PARAMETER LIST  ====#
PL <- list()

PL$N.var <- 9
PL$var.names <- c("CH2O.f","CH2O.s","O2","SO4","FeS","HCO3","NH4","NO3","HS")

PL$N.rate <- 6
PL$rate.names <- c("flux.up","flux.down","reac","irr","ddt","deficit")
PL$rate.summary <- data.frame(matrix(nrow=PL$N.var,ncol=PL$N.rate),row.names=PL$var.names)
names(PL$rate.summary) <- PL$rate.names

PL$N.reac <- 7
PL$reac.names <- c("R1","R2","R3","R4","R5","R6","R7")
PL$reaction.summary <- data.frame (matrix(nrow=PL$N.reac,ncol=1),row.names=PL$reac.names)
names(PL$reaction.summary) <- c("rate")

#======== MODEL DOMAIN AND GRID DEFINITION =========#
PL$L <- 15   # depth of sediment domain [cm]
PL$N <- 400  # number of grid layers
PL$grid <- setup.grid.1D(x.up=0, x.down=PL$L, N=PL$N, dx.1=PL$L/2000, p.dx.1=1.1)
PL$Depth <- PL$grid$x.mid

#======= PARAMETERS =========#
PL$S       <- 30    # salinity
PL$TC      <- 25    # temperature [deg C]
PL$P       <- 1.013 # pressure [bar]
PL$pH      <- 7.5

#porosity profile:
PL$por.0     <- 0.8     # porosity at the sediment-water interface
PL$por.inf   <- 0.8   # asymptotic porosity at depth
PL$por.x.att <- 1.0    # attenuation depth [cm]

PL$por.grid <- setup.prop.1D(func=p.exp, grid=PL$grid, y.0=PL$por.0,     y.inf=PL$por.inf,     x.L=0, x.att=PL$por.x.att)
PL$svf.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=(1-PL$por.0),y.inf=(1-PL$por.inf),x.L=0,x.att=PL$por.x.att)

#transport parameters:
PL$rho.sed <- 2.6     # density solid sediment [g cm-3]
PL$v.0     <- 0.2 # sedimentation velocity [cm yr-1]
PL$v.inf   <- 0.2 # sedimentation velocity [cm yr-1]
PL$velocity.info <- setup.compaction.1D(v.0 = PL$v.0, por.0=PL$por.0, por.inf=PL$por.inf, por.grid=PL$por.grid)
PL$v.grid <- PL$velocity.info$v
PL$u.grid <- PL$velocity.info$u

#Molecular diffusion (O2, SO4, HCO3, HS)
c.fac <- 10000*(3600*24*365.25)

Dmol.all <- c.fac*diffcoeff(S=PL$S, t=PL$TC, P=PL$P, species=c('O2', 'SO4', 'HCO3', 'NH4', 'NO3', 'HS', 'H', 'Fe'))

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
PL$k.f <- 10.0    # decay constant organic matter [yr-1]
PL$k.s <- 0.1     # decay constant organic matter [yr-1]

PL$K_O2  <- 0.008  # Monod constant O2 consumption  [umol cm-3 or mM]
PL$K_NO3 <- 0.008  # Monod constant NO3 reduction [umol cm-3 or mM]
PL$K_SO4 <- 0.9   # Monod constant SO4 reduction [umol cm-3 or mM]

PL$k_Sox <- 1E+06   # kinetic constant sulfide oxidation (umol-1 cm3 yr-1)
PL$k_NH4 <- 1E+06   # kinetic constant NH4 oxidation (umol-1 cm3 yr-1) -> nitrification?
PL$k_Sni <- 1E+06   # kinetic sulfide oxidation with nitrate (umol-1 cm3 yr-1)

PL$CNratio <- 106.0/16.0 # C to N ratio organic matter
PL$f.FeS   <- 0.1        # fraction of sulphide precipitating as FeS  

#======= BOUNDARY CONDITION ========#
#flux boundary conditions
frac.CH2O <- 0.5
PL$F.CH2O.f <- CH2O.tot*frac.CH2O       # fast degrading organic matter depostion (umol cm-2 yr-1)
PL$F.CH2O.s <- CH2O.tot*(1-frac.CH2O)   # slow degrading organic matter deposition (umol cm-2 yr-1)
PL$F.FeS    <- 0.0                      # FeS deposition [umol cm-2 yr-1]

#upper boundary concentrations (O2, SO4, and NO3 scale) 
O2.ow.max   <- 0.28  #O2 concentration bottom water [umol cm-3 or mM]
SO4.ow.max  <- 28.8  # SO4 concentration bottom water [umol cm-3 or mM]
NO3.ow.max  <- 0.1 # NO3 concentration bottom water [umol cm-3 or mM]
HCO3.ow.max <- 2.2   # HCO3 concentration bottom water [umol cm-3 or mM]
NH4.ow.max  <- 0.0   # NH4 concentration bottom water [umol cm-3 or mM]
HS.ow.max   <- 0.0   # HS concentration bottom water [umol cm-3 or mM]

O2.ow.min   <- 0.014  #O2 concentration bottom water [umol cm-3 or mM]     (1% PAL)
SO4.ow.min  <- 1      # SO4 concentration bottom water [umol cm-3 or mM]   (Canfield)
NO3.ow.min  <- 0.012    # NO3 concentration bottom water [umol cm-3 or mM] (FOAM)
HCO3.ow.min <- 2.2    # HCO3 concentration bottom water [umol cm-3 or mM]  (FOAM)
NH4.ow.min  <- 0.0    # NH4 concentration bottom water [umol cm-3 or mM]   (FOAM)
HS.ow.min   <- 0.0    # HS concentration bottom water [umol cm-3 or mM]    (FOAM)

iterations <- 4
#solute boundary conditions
O2.ow.seq    <- seq(from=O2.ow.min,    to=O2.ow.max,    length.out=iterations)
SO4.ow.seq   <- seq(from=SO4.ow.min,   to=SO4.ow.max,   length.out=iterations)
NO3.ow.seq   <- seq(from=NO3.ow.min,   to=NO3.ow.max,   length.out=iterations)
HCO3.ow.seq  <- seq(from=HCO3.ow.min,  to=HCO3.ow.max,  length.out=iterations)
NH4.ow.seq   <- seq(from=NH4.ow.min,   to=NH4.ow.max,   length.out=iterations)
HS.ow.seq    <- seq(from=HS.ow.min,    to=HS.ow.max,    length.out=iterations)

#======= BIOTURBATION PARAMETERS =======# 
#range Db from 0 to 3 cm2 yr-1

Db.seq <- seq(from=0, to=3, by=0.1)
#set up the rest in the loop

#Bioirrigation grid
PL$L.irr     <- 0  # irrigation depth [cm]
PL$irr.0     <- 0 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
PL$irr.inf   <- 0  # deep irrigation rate [yr-1]
PL$irr.x.att <- 3.0  # irrigation attenuation coef [cm]
PL$irr.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$L.irr,x.att=PL$irr.x.att)

#========= DATAFRAME SETUP FOR OUTPUT ========#
results_df <- as.data.frame(matrix(nrow=length(O2.ow.seq)*length(Db.seq), ncol=3))
colnames(results_df) <- c('Db', 'OPD', 'O2_bw')

#========= RUN ANALYSIS ========#
for(j in 1:length(O2.ow.seq)){ 
  
  #change boundary conditions here and set yini for initial conditions 
  PL$O2.ow <- O2.ow.seq[j]
  PL$NO3.ow <- NO3.ow.seq[j]
  PL$SO4.ow <- SO4.ow.seq[j]
  PL$HCO3.ow <- HCO3.ow.seq[j]
  PL$NH4.ow <- NH4.ow.seq[j]
  PL$HS.ow  <- HS.ow.seq[j]
  
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
  
  for(i in 1:length(Db.seq)){
    
    #set up biodiffusion grid
    PL$Db.0 <- Db.seq[i]
    PL$L.mix    <- 1 + 9*(1-exp((-PL$Db.0/3))) # Mixed depth layer
    PL$Db.inf   <- 0.0      # deep Db [cm2 yr-1]
    PL$Db.x.att <- 2      # attenuation depth [cm]
    PL$Db.grid  <- setup.prop.1D(func=p.sig, grid=PL$grid,y.0=PL$Db.0,y.inf=PL$Db.inf,x.L=PL$L.mix,x.att=PL$Db.x.att)
    
    #====== Run model to steady state =====#
    out <- steady.1D(y=yini, func=model, parms=PL, nspec=PL$N.var, pos=TRUE)
    
    #====== Measure OPD =====#
    O2.out <- out$y[(2*PL$N+1):(3*PL$N)] #get output of O2
    O2.search <- max(which(O2.out>1.0e-3))
    O2.penetrationdepth <- PL$Depth[O2.search]
    
    #====== save output in results_df ======#
    #put into dataframe
    idx <- match(NA, results_df$O2_bw)
    results_df$Db[idx]     <- PL$Db.0
    results_df$OPD[idx]    <- O2.penetrationdepth
    results_df$O2_bw[idx]  <- round(PL$O2.ow, digits=3)
    
    
  }
}

#save as appropriate output file:
save(results_df, file='BiodiffusionSensitivityTestOutput_low.RData')
#save(results_df, file='BiodiffusionSensitivityTestOutput_med.RData')
#save(results_df, file='BiodiffusionSensitivityTestOutput_high.RData')

#==== SENSITIVITY ANALYSIS FOR BIOIRRIGATION COEFFICIENTS =====#
#Generate sensitivity analysis curve by increasing irr.0 from 0 to 300
#Parameters and boundary conditions the same as above, so we only need to change the Db.0 and irr.0 sequences

#======= BIOTURBATION PARAMETERS =======# 
#Db.0 come from Ediacaran_Db and Terreneuvian_Db
Ediacaran_Db <- 0.1
Terreneuvian_Db <- 0.98

# #Bioirrigation grid equence
irr.0.seq <- seq(from=0, to=300, by=10)

#========= DATAFRAME SETUP FOR OUTPUT ========#
#Two dataframes, one for Ediacaran_Db and the other for Terreneuvian_Db
#record OPD, irr.0, and changing bottom water O2 

Ediacaran_out <- as.data.frame(matrix(NA, nrow=length(O2.ow.seq)*length(irr.0.seq), ncol=3))
colnames(Ediacaran_out) <- c('irr', 'O2_bw', 'OPD')

Terreneuvian_out <- as.data.frame(matrix(NA, nrow=length(O2.ow.seq)*length(irr.0.seq), ncol=3))
colnames(Terreneuvian_out) <- c('irr', 'O2_bw', 'OPD')

#========= RUN ANALYSIS - EDIACARAN ========#
#set up biodiffusion grid
PL$Db.0 <- Ediacaran_Db
PL$L.mix    <- 1 + 9*(1-exp((-PL$Db.0/3))) # Mixed depth layer
PL$Db.inf   <- 0.0      # deep Db [cm2 yr-1]
PL$Db.x.att <- 2      # attenuation depth [cm]
PL$Db.grid  <- setup.prop.1D(func=p.sig, grid=PL$grid,y.0=PL$Db.0,y.inf=PL$Db.inf,x.L=PL$L.mix,x.att=PL$Db.x.att)

for(j in 1:length(O2.ow.seq)){ 
  
  #change boundary conditions here and set yini for initial conditions 
  PL$O2.ow <- O2.ow.seq[j]
  PL$NO3.ow <- NO3.ow.seq[j]
  PL$SO4.ow <- SO4.ow.seq[j]
  PL$HCO3.ow <- HCO3.ow.seq[j]
  PL$NH4.ow <- NH4.ow.seq[j]
  PL$HS.ow  <- HS.ow.seq[j]
  
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
  
  for(i in 1:length(irr.0.seq)){
    
    #Bioirrigation grid
    PL$L.irr     <- 0  # irrigation depth [cm]
    PL$irr.0     <- irr.0.seq[i] # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
    PL$irr.inf   <- 0  # deep irrigation rate [yr-1]
    PL$irr.x.att <- 3.0  # irrigation attenuation coef [cm]
    PL$irr.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$L.irr,x.att=PL$irr.x.att)
    
    #====== Run model to steady state =====#
    out <- steady.1D(y=yini, func=model, parms=PL, nspec=PL$N.var, pos=TRUE)
    
    #====== Measure OPD =====#
    O2.out <- out$y[(2*PL$N+1):(3*PL$N)] #get output of O2
    O2.search <- max(which(O2.out>1.0e-3))
    O2.penetrationdepth <- PL$Depth[O2.search]
    
    #====== save output in results_df ======#
    #put into dataframe
    idx <- match(NA, Ediacaran_out$O2_bw)
    Ediacaran_out$irr[idx]     <- PL$irr.0
    Ediacaran_out$OPD[idx]    <- O2.penetrationdepth
    Ediacaran_out$O2_bw[idx]  <- round(PL$O2.ow, digits=3)
    
    
  }
}

#save as appropriate file
save(Ediacaran_out, file='Ediacaran_BioirrigationSensitivityTest_low.RData')
#save(Ediacaran_out, file='Ediacaran_BioirrigationSensitivityTest_med.RData')
#save(Ediacaran_out, file='Ediacaran_BioirrigationSensitivityTest_high.RData')

#==== SAME ANALYSIS - TERRENEUVIAN ====#
#set up biodiffusion grid
PL$Db.0 <- Terreneuvian_Db
PL$L.mix    <- 1 + 9*(1-exp((-PL$Db.0/3))) # Mixed depth layer
PL$Db.inf   <- 0.0      # deep Db [cm2 yr-1]
PL$Db.x.att <- 2      # attenuation depth [cm]
PL$Db.grid  <- setup.prop.1D(func=p.sig, grid=PL$grid,y.0=PL$Db.0,y.inf=PL$Db.inf,x.L=PL$L.mix,x.att=PL$Db.x.att)

for(j in 1:length(O2.ow.seq)){ 
  
  #change boundary conditions here and set yini for initial conditions 
  PL$O2.ow <- O2.ow.seq[j]
  PL$NO3.ow <- NO3.ow.seq[j]
  PL$SO4.ow <- SO4.ow.seq[j]
  PL$HCO3.ow <- HCO3.ow.seq[j]
  PL$NH4.ow <- NH4.ow.seq[j]
  PL$HS.ow  <- HS.ow.seq[j]
  
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
  
  for(i in 1:length(irr.0.seq)){
    
    #Bioirrigation grid
    PL$L.irr     <- 0  # irrigation depth [cm]
    PL$irr.0     <- irr.0.seq[i] # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
    PL$irr.inf   <- 0  # deep irrigation rate [yr-1]
    PL$irr.x.att <- 3.0  # irrigation attenuation coef [cm]
    PL$irr.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$L.irr,x.att=PL$irr.x.att)
    
    #====== Run model to steady state =====#
    out <- steady.1D(y=yini, func=model, parms=PL, nspec=PL$N.var, pos=TRUE)
    
    #====== Measure OPD =====#
    O2.out <- out$y[(2*PL$N+1):(3*PL$N)] #get output of O2
    O2.search <- max(which(O2.out>1.0e-3))
    O2.penetrationdepth <- PL$Depth[O2.search]
    
    #====== save output in results_df ======#
    #put into dataframe
    idx <- match(NA, Terreneuvian_out$O2_bw)
    Terreneuvian_out$irr[idx]     <- PL$irr.0
    Terreneuvian_out$OPD[idx]    <- O2.penetrationdepth
    Terreneuvian_out$O2_bw[idx]  <- round(PL$O2.ow, digits=3)
    
    
  }
}

#save as appropriate output file
save(Terreneuvian_out, file='Terreneuvian_BioirrigationSensitivityTest_low.RData')
#save(Terreneuvian_out, file='Terreneuvian_BioirrigationSensitivityTest_med.RData')
#save(Terreneuvian_out, file='Terreneuvian_BioirrigationSensitivityTest_high.RData')
