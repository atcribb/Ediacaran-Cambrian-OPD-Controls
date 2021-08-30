#Authors: Sebastiaan J. van de Velde, Alison T. Cribb
#Simulation: Weak biomixing + weak bioirrigation 
#Use biogeochemical model function from Model_Function.R script 

#====== PACKAGES ========#
library(ReacTran)
library(rootSolve)
library(marelac)
library(writexl)

#====== NO DATA INPUT ======#

#====== BIOTURBATION OPTIONS =======#
#modify bioturbation settings here
biomixing <- 1.0 #biodiffusion coefficient (cm2 yr-1)
bioirrigation <- 36.5 #bioirrigation coefficient (yr-1)

#====== FUNCTIONS ======#
# function to calculate tortuiosity:
f.tort <- function(por) 1-2*log(por)
tort   <- function(D.grid, por.grid){
  D.grid$mid <- D.grid$mid/f.tort(por.grid$mid)
  D.grid$int <- D.grid$int/f.tort(por.grid$int)
  return(D.grid)
}

# make sure model function from Model_Function is in environment 

#====== PARAMETERS LIST =====#
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

#======= BOUNDARY CONDITIONS=========#
umolcmyr <- (365.25/10) #conver to umol cm-2 y-1

#MAX BOUNDARY CONDITIONS (Modern)
#flux boundary conditions
F.CH2O.tot.max    <- 547.875                    #total organic matter deposition (umol cm-2 yr-1)
f.CH2O.fast.max   <- 0.50                       #fraction of organic matter that is fast degrading
F.CH2O.f.max   <- F.CH2O.tot.max*(f.CH2O.fast.max)   #fast-remineralizing organic matter deposition (umol cm-2 yr-1)
F.CH2O.s.max   <- F.CH2O.tot.max*(1-f.CH2O.fast.max) #slow-remineralizing organic matter deposition (umol cm-2 yr-1)
F.FeS.max      <- 0.0                        # FeS deposition (umol cm-2 yr-1)

#solute boundary conditions 
O2.ow.max   <- 0.28  #O2 concentration bottom water (umol cm-3 or mM) 
SO4.ow.max  <- 28.8  # SO4 concentration bottom water (umol cm-3 or mM) 
NO3.ow.max  <- 0.1 # NO3 concentration bottom water (umol cm-3 or mM) 
HCO3.ow.max <- 2.2   # HCO3 concentration bottom water (umol cm-3 or mM) 
NH4.ow.max  <- 0.0   # NH4 concentration bottom water (umol cm-3 or mM) 
HS.ow.max   <- 0.0   # HS concentration bottom water (umol cm-3 or mM) 

#MIN BOUNDARY CONDITIONS (Neoproterozoic)
#flux boundary conditions
F.CH2O.tot.min    <- 27.4                       #total organic matter deposition (umol cm-2 yr-1)
f.CH2O.fast.min   <- 0.50                       #fraction of organic matter that is fast degrading 
F.CH2O.f.min      <- F.CH2O.tot.min*(f.CH2O.fast.min)   #fast-remineralizing organic matter deposition (umol cm-2 yr-1)
F.CH2O.s.min      <- F.CH2O.tot.min*(1-f.CH2O.fast.min) #slow-remineralizing organic matter deposition (umol cm-2 yr-1)
F.FeS.min         <- 0.0                        # FeS deposition (umol cm-2 yr-1)

#solute boundary conditions 
O2.ow.min   <- 0.014  #O2 concentration bottom water (umol cm-3 or mM)     (1% PAL)
SO4.ow.min  <- 1      # SO4 concentration bottom water (umol cm-3 or mM)   - Canfield and Farquhar (2009)
NO3.ow.min  <- 0.012    # NO3 concentration bottom water (umol cm-3 or mM) - FOAM from Zhao et al. (2020)
HCO3.ow.min <- 2.2    # HCO3 concentration bottom water (umol cm-3 or mM)  - FOAM from Zhao et al. (2020)
NH4.ow.min  <- 0.0    # NH4 concentration bottom water (umol cm-3 or mM)   - FOAM from Zhao et al. (2020)
HS.ow.min   <- 0.0    # HS concentration bottom water (umol cm-3 or mM)    - FOAM from Zhao et al. (2020)

#ESTALBISH BOUNDARY CONDITION SEQUENCES
iterations <- 100 #length of each sequence vector -- interations^2 = number of model runs 
#flux boundary conditions
F.CH2O.f.seq <- seq(from=F.CH2O.f.min, to=F.CH2O.f.max, length.out=iterations)
F.CH2O.s.seq <- seq(from=F.CH2O.s.min, to=F.CH2O.s.max, length.out=iterations)
F.FeS.seq    <- seq(from=F.FeS.min,    to=F.FeS.max,    length.out=iterations)

#solute boundary conditions
O2.ow.seq    <- seq(from=O2.ow.min,    to=O2.ow.max,    length.out=iterations)
SO4.ow.seq   <- seq(from=SO4.ow.min,   to=SO4.ow.max,   length.out=iterations)
NO3.ow.seq   <- seq(from=NO3.ow.min,   to=NO3.ow.max,   length.out=iterations)
HCO3.ow.seq  <- seq(from=HCO3.ow.min,  to=HCO3.ow.max,  length.out=iterations)
NH4.ow.seq   <- seq(from=NH4.ow.min,   to=NH4.ow.max,   length.out=iterations)
HS.ow.seq    <- seq(from=HS.ow.min,    to=HS.ow.max,    length.out=iterations)

#set PL$Boundary Conditions in the sensitivty analysis loop

#====== BIOTURBATION PARAMETERS (no biorbation on this analysis) =======#
#takes biomixing and bioirrigation variables from top of the script (lines 14-15)
#or they can be manually changed here for PL$Db.0 and PL$irr.0

#Biodiffusion grid
PL$Db.0     <- biomixing # biodiffusion coefficient Db at SWI 
PL$L.mix    <- 1 + 9*(1-exp((-PL$Db.0/3))) # Mixed depth layer - Eq from van de Velde and Meysman (2016)
PL$Db.inf   <- 0.0      # deep Db (cm2 yr-1)
PL$Db.x.att <- 2.0      # attenuation depth (cm)
PL$Db.grid  <- setup.prop.1D(func=p.sig, grid=PL$grid,y.0=PL$Db.0,y.inf=PL$Db.inf,x.L=PL$L.mix,x.att=PL$Db.x.att)

#Bioirrigation grid
PL$L.irr     <- 0    # irrigation depth (cm)
PL$irr.0     <- bioirrigation # irrigation rate at SWI (yr-1)
PL$irr.inf   <- 0    # deep irrigation rate (yr-1)
PL$irr.x.att <- 3.0  # irrigation attenuation coefficient ([cm])cm)
PL$irr.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$L.irr,x.att=PL$irr.x.att)

#===================================================#
#======== MODEL RUN: SENSITIVITY ANALYSIS ===========#
#===================================================#
#results dataframe setup:
results.df <- as.data.frame(matrix(data=NA, nrow=length(O2.ow.seq)*length(F.CH2O.f.seq), ncol=4))
colnames(results.df) <- c('Corg_fast', 'Corg_slow', 'O2_bw', 'O2_depth')

for(i in 1:length(F.CH2O.f.seq)){ #loop through Corg flux
  
  print(paste('CH2O SEQUENCE NUMBER', i, 'OUT OF', length(F.CH2O.f.seq))) #counter for model runs - 100 100 CH2O sequences (10,000 total)
  
  PL$F.CH2O.f <- F.CH2O.f.seq[i]
  PL$F.CH2O.s <- F.CH2O.s.seq[i]
  
  for(j in 1:length(O2.ow.seq)){ #loop through 100 oxygen concentrationss
    
    #change oxygen boundary condition 
    PL$O2.ow <- O2.ow.seq[j]
    
    #get rest of boundary conditions
    PL$F.FeS   <- F.FeS.seq[j]
    PL$SO4.ow  <- SO4.ow.seq[j]
    PL$NO3.ow  <- NO3.ow.seq[j]
    PL$HCO3.ow <- HCO3.ow.seq[j]
    PL$NH4.ow  <- NH4.ow.seq[j]
    PL$HS.ow   <- HS.ow.seq[j]
    
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
    yini[(8*N+1):(9*N)] <- HS
    
    
    #run steady state solution
    out <- steady.1D(y=yini, func=model, parms=PL, nspec=PL$N.var, pos=TRUE)
    
    #only need to pull the O2 out of the results
    O2.out <- out$y[(2*PL$N+1):(3*PL$N)]
    
    #find where oxygen goes below 1.0e-3
    O2.search <- max(which(O2.out>1.0e-3))
    O2.penetrationdepth <- PL$Depth[O2.search]
    
    #put into dataframe
    idx <- match(NA, results.df$O2_bw)
    
    results.df$Corg_fast[idx] <- PL$F.CH2O.f
    results.df$Corg_slow[idx] <- PL$F.CH2O.s
    results.df$O2_bw[idx] <- PL$O2.ow
    results.df$O2_depth[idx] <- O2.penetrationdepth
    
  }
  
}

#save results 
save(results.df, file="Sim4_WeakBioturbation.RData")
write_xlsx(results.df, "Output_Sim4_WeakBioturbation.xlsx")
print('FINISHED SIMULATION 4')

