#Authors: Alison Cribb, Sebastiaan van de Velde
#Plots oxygen consumption profiles, produced data for Supplementary Material 6
#Use with Model_Function.R from Model directory 

#====== PACKAGES ========#
require(ReacTran)
require(marelac)
require(AquaEnv)

#====== NO DATA INPUT ========#

#===== FUNCTIONS =======#
# Calculate tortuosity
f.tort <- function(por) 1-2*log(por)
tortuosity.correction <- function (D.grid,por.grid)
{
  D.grid$mid <- D.grid$mid/f.tort(por.grid$mid)
  D.grid$int <- D.grid$int/f.tort(por.grid$int)
  return(D.grid)
}

#function fo check mass balances
check.balance <- function(tran,reac,irr,ddt,VF)
{
  r <- vector(length=6)
  r[1] <- tran$flux.up
  r[2] <- tran$flux.down
  r[3] <- sum(VF$mid*reac*grid$dx)
  r[4] <- sum(VF$mid*irr*grid$dx)
  r[5] <- sum(VF$mid*ddt*grid$dx)
  r[6] <- r[1] - r[2] + r[3] + r[4] - r[5]
  return(r)
}

#======= PARAMETERS LIST =====#
parameters <- c()

N.var <- 9
var.names <- c("CH2O.f","CH2O.s","O2","SO4","FeS","HCO3","NH4","NO3","HS")

N.rate <- 6
rate.names <- c("flux.up","flux.down","reac","irr","ddt","deficit")
rate.summary <- data.frame(matrix(nrow=N.var,ncol=N.rate),row.names=var.names)
names(rate.summary) <- rate.names

N.reac <- 7
reac.names <- c("R1","R2","R3","R4","R5","R6","R7")
reaction.summary <- data.frame (matrix(nrow=N.reac,ncol=1),row.names=reac.names)
names(reaction.summary) <- c("rate")

#======== MODEL DOMAIN AND GRID DEFINITION =========#

L <- 15   # depth of sediment domain [cm]
N <- 400  # number of grid layers
grid <- setup.grid.1D(x.up=0, x.down=L, N=N, dx.1=L/2000, p.dx.1=1.1)
Depth <- grid$x.mid

#======= PARAMETERS =========#
# Environmental parameters
S       <- 30    # salinity
TC      <- 25    # temperature [deg C]
P       <- 1.013 # pressure [bar]
pH      <- 7.5

# Porosity profile=
por.0     <- 0.8     # porosity at the sediment-water interface
por.inf   <- 0.8   # asymptotic porosity at depth
por.x.att <- 1.0    # attenuation depth [cm]

por.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=por.0,y.inf=por.inf,x.L=0,x.att=por.x.att)
svf.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=(1-por.0),y.inf=(1-por.inf),x.L=0,x.att=por.x.att)

#transport parameters:
rho.sed <- 2.6     # density solid sediment [g cm-3]
v.0     <- 0.2 # sedimentation velocity [cm yr-1]
v.inf   <- 0.2 # sedimentation velocity [cm yr-1]
velocity.info <- setup.compaction.1D(v.0 = v.0, por.0=por.0, por.inf=por.inf, por.grid=por.grid)
v.grid <- velocity.info$v
u.grid <- velocity.info$u

#Molecular diffusion (O2, SO4, HCO3, HS)
c.fac <- 10000*(3600*24*365.25)

Dmol.O2   <- c.fac*diffcoeff(S=S,t=TC,P=P,species="O2")$O2
D.O2.grid <- setup.prop.1D(value=Dmol.O2,grid=grid)
D.O2.grid <- tortuosity.correction(D.O2.grid,por.grid)

Dmol.NO3   <- c.fac*diffcoeff(S=S,t=TC,P=P,species="NO3")$NO3
D.NO3.grid <- setup.prop.1D(value=Dmol.NO3,grid=grid)
D.NO3.grid <- tortuosity.correction(D.NO3.grid,por.grid)

Dmol.SO4   <- c.fac*diffcoeff(S=S,t=TC,P=P,species="SO4")$SO4
D.SO4.grid <- setup.prop.1D(value=Dmol.SO4,grid=grid)
D.SO4.grid <- tortuosity.correction(D.SO4.grid,por.grid)

Dmol.HCO3   <- c.fac*diffcoeff(S=S,t=TC,P=P,species="HCO3")$HCO3
D.HCO3.grid <- setup.prop.1D(value=Dmol.HCO3,grid=grid)
D.HCO3.grid <- tortuosity.correction(D.HCO3.grid,por.grid)

Dmol.NH4   <- c.fac*diffcoeff(S=S,t=TC,P=P,species="NH4")$NH4
D.NH4.grid <- setup.prop.1D(value=Dmol.NH4,grid=grid)
D.NH4.grid <- tortuosity.correction(D.NH4.grid,por.grid)

Dmol.HS    <- c.fac*diffcoeff(S=S,t=TC,P=P,species="HS")$HS
D.HS.grid  <- setup.prop.1D(value=Dmol.HS,grid=grid)
D.HS.grid  <- tortuosity.correction(D.HS.grid,por.grid)


#========= REACTION PARAMETERS =========#
#organic matter decay
k.f <- 10.0    # decay constant organic matter [yr-1]
k.s <- 0.1     # decay constant organic matter [yr-1]

K_O2  <- 0.008  # Monod constant O2 consumption  [umol cm-3 or mM]
K_NO3 <- 0.008  # Monod constant NO3 reduction [umol cm-3 or mM]
K_SO4 <- 0.9   # Monod constant SO4 reduction [umol cm-3 or mM]

k_Sox <- 1E+06   # kinetic constant sulfide oxidation (umol-1 cm3 yr-1)
k_NH4 <- 1E+06   # kinetic constant NH4 oxidation (umol-1 cm3 yr-1) -> nitrification?
k_Sni <- 1E+06   # kinetic sulfide oxidation with nitrate (umol-1 cm3 yr-1)

CNratio <- 106.0/16.0 # C to N ratio organic matter
f.FeS   <- 0.1        # fraction of sulphide precipitating as FeS

#======= BOUNDARY CONDITIONS=========#
umolcmyr <- (365.25/10) #conver to umol cm-2 y-1

#==UNCHANGING BOUNDARY CONDITIONS==#
#flux boundary conditions
F.CH2O.tot <- 450
f.CH2O.fast <- 0.5
F.CH2O.f <- F.CH2O.tot*f.CH2O.fast
F.CH2O.s <- F.CH2O.tot*(1-f.CH2O.fast)
F.FeS    <- 0.0         # FeS deposition [umol cm-2 yr-1]

#upper boundary concentrations
O2.ow   <- 0.14    # O2 concentration bottom water [umol cm-3 or mM]
SO4.ow  <- 22      # SO4 concentration bottom water [umol cm-3 or mM]
HCO3.ow <- 2.0     # HCO3 concentration bottom water [umol cm-3 or mM]
NH4.ow  <- 0.0     # NH4 concentration bottom water [umol cm-3 or mM]
HS.ow   <- 0.0     # Fe concentration bottom water [umol cm-3 or mM]
NO3.ow  <- 0.012   # NO3 concentration bottom water [umol cm-3 or mM]




#===================================================#
#====== MODEL RUN: OXYGEN CONSUMPTION PROFILES ======#
#===================================================#

#====== BIOTURBATION PARAMETERS: NO BIOTURBATION =======#
#Biodiffusion grid
Db.0   <- 0.0      # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0  # irrigation depth [cm]
irr.0     <- 0 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.nobio <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.nobio <- results.nobio$rate.matrix
reac.summary.nobio <- results.nobio$reac.matrix
results.nobio$y <- out$y
save(results.nobio, file='O2cons_nobioturbation.RData')

#====== BIOTURBATION PARAMETERS: EDIACARAN BIOMIXING =======#
#Biodiffusion grid
Db.0   <- 0.1    # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0  # irrigation depth [cm]
irr.0     <- 0.0 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.edibiomix <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.edibiomix <- results.edibiomix$rate.matrix
reac.summary.edibiomix <- results.edibiomix$reac.matrix
results.edibiomix$y <- out$y
save(results.edibiomix, file='O2cons_ediacaranbiomixing.RData')

#====== BIOTURBATION PARAMETERS: EDIACARAN BIOIRRIGATION =======#
#Biodiffusion grid
Db.0   <- 0    # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0.0  # irrigation depth [cm]
irr.0     <- 3.65 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.edibioirr <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.edibioirr <- results.edibioirr$rate.matrix
reac.summary.edibioirr <- results.edibioirr$reac.matrix
results.edibioirr$y <- out$y
save(results.edibioirr, file='O2cons_ediacaranbioirrigation.RData')


#====== BIOTURBATION PARAMETERS: TERRENEUVIAN BIOTURBATION =======#
#Biodiffusion grid
Db.0   <- 0.1   # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0  # irrigation depth [cm]
irr.0     <- 3.65 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.ediacaran <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.ediacaran <- results.ediacaran$rate.matrix
reac.summary.ediacaran <- results.ediacaran$reac.matrix
results.ediacaran$y <- out$y
save(results.ediacaran, file='O2cons_ediacaranbioturbation.RData')

#====== BIOTURBATION PARAMETERS: TERRENEUVIAN BIOMIXING =======#
#Biodiffusion grid
Db.0   <- 0.98    # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0  # irrigation depth [cm]
irr.0     <- 0.0 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.terrebiomix <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.terrebiomix <- results.terrebiomix$rate.matrix
reac.summary.terrebiomix <- results.terrebiomix$reac.matrix
results.terrebiomix$y <- out$y
save(results.terrebiomix, file='O2cons_terreneuvianbiomixing.RData')

#====== BIOTURBATION PARAMETERS: TERRENEUVIAN BIOIRRIGATION =======#
#Biodiffusion grid
Db.0   <- 0    # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0.0  # irrigation depth [cm]
irr.0     <- 35.77 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.terrebioirr <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.terrebioirr <- results.terrebioirr$rate.matrix
reac.summary.terrebioirr <- results.terrebioirr$reac.matrix
results.terrebioirr$y <- out$y
save(results.terrebioirr, file='O2cons_terreneuvianbioirrigation.RData')


#====== BIOTURBATION PARAMETERS: TERRENEUVIAN BIOTURBATION =======#
#Biodiffusion grid
Db.0   <- 0.98    # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0  # irrigation depth [cm]
irr.0     <- 35.77 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.terreneuvian <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.terreneuvian <- results.terreneuvian$rate.matrix
reac.summary.terreneuvian <- results.terreneuvian$reac.matrix
results.terreneuvian$y <- out$y
save(results.terreneuvian, file='O2cons_terreneuvianbioturbation.RData')

#====== BIOTURBATION PARAMETERS: MODERN BIOMIXING =======#
#Biodiffusion grid
Db.0   <- 10     # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0  # irrigation depth [cm]
irr.0     <- 0.0 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.modernbiomix <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.modernbiomix <- results.modernbiomix$rate.matrix
reac.summary.modernbiomix <- results.modernbiomix$reac.matrix
results.modernbiomix$y <- out$y
save(results.modernbiomix, file='O2cons_modernbiomixing.RData')

#====== BIOTURBATION PARAMETERS: MODERN BIOIRRIGATION =======#
#Biodiffusion grid
Db.0   <- 0.0     # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0  # irrigation depth [cm]
irr.0     <- 365 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.modernbioirr <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.modernbioirr <- results.modernbioirr$rate.matrix
reac.summary.modernbioirr <- results.modernbioirr$reac.matrix
results.modernbioirr$y <- out$y
save(results.modernbioirr, file='O2cons_modernbioirrigation.RData')

#====== BIOTURBATION PARAMETERS: MODERN BIOTURBTAION =======#
#Biodiffusion grid
Db.0   <- 10.0     # biodiffusion coefficient Db at SWI [cm2 yr-1]
L.mix  <- 1 + 9*(1-exp((-Db.0/3))) # Mixed depth layer
Db.inf <- 0.0      # asymptotic Db at depth cm2 yr-1
Db.x.att <- 2      # attenuation depth [cm]
Db.grid <- setup.prop.1D(func=p.sig,grid=grid,y.0=Db.0,y.inf=Db.inf,x.L=L.mix,x.att=Db.x.att)

#Bioirrigation grid
L.irr     <- 0  # irrigation depth [cm]
irr.0     <- 365 # irrigation rate at SWI [yr-1] #weak=36.5, strong=365
irr.inf   <- 0  # deep irrigation rate [yr-1]
irr.x.att <- 3.0  # irrigation attenuation coef [cm]
irr.grid <- setup.prop.1D(func=p.exp,grid=grid,y.0=irr.0,y.inf=irr.inf,x.L=L.irr,x.att=irr.x.att)

#Initialize state variable vector
CH2O.f  <- 0
CH2O.s  <- 0
O2    <- O2.ow
SO4   <- SO4.ow
FeS   <- 0
HCO3  <- HCO3.ow
NH4   <- NH4.ow
NO3   <- NO3.ow
HS    <- HS.ow

yini <- vector(length=N.var*N)
yini[(0*N+1):(1*N)] <- CH2O.f
yini[(1*N+1):(2*N)] <- CH2O.s
yini[(2*N+1):(3*N)] <- O2
yini[(3*N+1):(4*N)] <- SO4
yini[(4*N+1):(5*N)] <- FeS
yini[(5*N+1):(6*N)] <- HCO3
yini[(6*N+1):(7*N)] <- NH4
yini[(7*N+1):(8*N)] <- NO3
yini[(8*N+1):(9*N)] <- HS

#run steady state 
out <- steady.1D(y=yini, func=model, parms=parameters, nspec=N.var, pos=TRUE)

results.modern <- model(t=0,state=out$y,parameters=parameters,full.output=TRUE)
rate.summary.modern<- results.modern$rate.matrix
reac.summary.modern <- results.modern$reac.matrix
results.modern$y <- out$y
save(results.modern, file='O2cons_modernbioturbation.RData')
