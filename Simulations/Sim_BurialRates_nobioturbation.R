#Authors: ALison Cribb, Sebastiaan van de Velde
#Code to simulate burial rates of organic carbon and reduced compounds -- no bioturbation
#Need Model_Function.R loaded in to run analyses 

# #====== PACKAGES ========#
library(ReacTran)
library(rootSolve)
library(marelac)
library(writexl)

#====== BIOTURBATION PARAMETERS ======#
Db.0.input <- 0
a0.0.input <- 0

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

#======== BIOTURBATION GRID =======#
#biodiffusion transport
PL$Db.0     <- Db.0.input    # biodiffusion coefficient Db at SWI [cm2 yr-1]
PL$L.mix    <- 1 + 9*(1-exp((-PL$Db.0/3))) # Mixed depth layer
PL$Db.inf   <- 0.0      # deep Db [cm2 yr-1]
PL$Db.x.att <- 2      # attenuation depth [cm]
PL$Db.grid  <- setup.prop.1D(func=p.sig, grid=PL$grid,y.0=PL$Db.0,y.inf=PL$Db.inf,x.L=PL$L.mix,x.att=PL$Db.x.att)

#bioirrigation transport
PL$L.irr     <- 0.0  # irrigation depth [cm]
PL$irr.0     <- a0.0.input  # irrigation rate at SWI [yr-1]
PL$irr.inf   <- 0.0  # deep irrigation rate [yr-1]
PL$irr.x.att <- 3.0  # irrigation attenuation coef [cm]
PL$irr.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$L.irr,x.att=PL$irr.x.att)

#======== BOUNDARY CONDITIONS - SEQUENCES =========#
#flux boundary conditions
CH2O.tot     <- c(150, 300, 450, 700)
frac.CH2O    <- 0.5 #ratio of labile:refractory organic matter 
F.CH2O.f.seq <- frac.CH2O*CH2O.tot        # organic matter deposition (umol cm-2 yr-1)
F.CH2O.s.seq <- (1-frac.CH2O)*CH2O.tot    # organic matter deposition (umol cm-2 yr-1)
PL$F.FeS     <- 0.0                       # FeS deposition (umol cm-2 yr-1)

#upper boundary concentrations (O2, SO4, and NO3 scale) 
O2.ow.max   <- 0.28  #O2 concentration bottom water [umol cm-3 or mM]
SO4.ow.max  <- 28.8  # SO4 concentration bottom water [umol cm-3 or mM]
NO3.ow.max  <- 0.1 # NO3 concentration bottom water [umol cm-3 or mM]
HCO3.ow.max <- 2.2   # HCO3 concentration bottom water [umol cm-3 or mM]
NH4.ow.max  <- 0.0   # NH4 concentration bottom water [umol cm-3 or mM]
HS.ow.max   <- 0.0   # HS concentration bottom water [umol cm-3 or mM]

#sequence to scale boundary conditions - solutes scale linearly with PAL
PALscaleseq <- c(0.05, 0.1, 0.25, 0.50, 0.70, 1) 

#solute boundary conditions
O2.ow.seq    <- O2.ow.max*PALscaleseq
SO4.ow.seq   <- SO4.ow.max*PALscaleseq
NO3.ow.seq   <- NO3.ow.max*PALscaleseq
HCO3.ow.seq  <- rep(HCO3.ow.max, length(PALscaleseq)) #keep HCO3 the same
NH4.ow.seq   <- NH4.ow.max*PALscaleseq
HS.ow.seq    <- HS.ow.max*PALscaleseq

#========= DATAFRAME SETUP FOR OUTPUT ========#
results_df <- as.data.frame(matrix(nrow=length(O2.ow.seq)*length(F.CH2O.f.seq), ncol=9))
colnames(results_df) <- c('Db', 'irr', 'CH2O_tot', 'O2_bw', 
                          'burial_CH2O.f', 'burial_CH2O.s', 'burial_FeS', 'burial_NH4', 'burial_HS')
results_df$Db <- Db.0.input
results_df$irr <- a0.0.input 

#=================================#
#             MODEL RUN           #
#=================================#

#full boundary condition sequences
O2.seq.full <- rep(O2.ow.seq, length(F.CH2O.f.seq))
SO4.seq.full <- rep(SO4.ow.seq, length(F.CH2O.f.seq))
NO3.seq.full <- rep(NO3.ow.seq, length(F.CH2O.f.seq))
HCO3.seq.full <- rep(HCO3.ow.seq, length(F.CH2O.f.seq))
NH4.seq.full <- rep(NH4.ow.seq, length(F.CH2O.f.seq))
HS.seq.full <- rep(HS.ow.seq, length(F.CH2O.f.seq))

CH2O.f.seq.full <- c(rep(F.CH2O.f.seq[1], length(O2.ow.seq)),
                     rep(F.CH2O.f.seq[2], length(O2.ow.seq)),
                     rep(F.CH2O.f.seq[3], length(O2.ow.seq)),
                     rep(F.CH2O.f.seq[4], length(O2.ow.seq)))
CH2O.s.seq.full <- c(rep(F.CH2O.s.seq[1], length(O2.ow.seq)),
                     rep(F.CH2O.s.seq[2], length(O2.ow.seq)),
                     rep(F.CH2O.s.seq[3], length(O2.ow.seq)),
                     rep(F.CH2O.s.seq[4], length(O2.ow.seq)))

results_df$CH2O_tot <- CH2O.f.seq.full + CH2O.s.seq.full
results_df$O2_bw <- O2.seq.full

for(i in 1:nrow(results_df)){  

    #solute boundary conditions  
    PL$O2.ow <- O2.seq.full[i]
    PL$SO4.ow <- SO4.seq.full[i]
    PL$NO3.ow <- NO3.seq.full[i]
    PL$HCO3.ow <- HCO3.seq.full[i]
    PL$NH4.ow <- NH4.seq.full[i]
    PL$HS.ow <- HS.seq.full[i]
    
    PL$F.CH2O.f <- CH2O.f.seq.full[i]
    PL$F.CH2O.s <- CH2O.f.seq.full[i]
  
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
    
    out <- steady.1D(y=yini, func=model, parms=PL, nspec=PL$N.var, pos=TRUE)
    
    #save outputs in results_df 
    results_df$burial_CH2O.f[i] <- out$CH2O.f.burial
    results_df$burial_CH2O.s[i] <- out$CH2O.s.burial
    results_df$burial_FeS[i]    <- out$FeS.burial
    results_df$burial_NH4[i]    <- out$NH4.burial
    results_df$burial_HS[i]     <- out$HS.burial 
    
}

View(results_df)
save(results_df, file='nobio_burial_df.RData')