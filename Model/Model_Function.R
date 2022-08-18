#Model function simulating carbon, oxygen, nitrogen, and sulfur reactions
#Authors: Sebastiaan van de Velde, Alison Cribb

#===== REACTIONS ========#
#Primary reactions - organic matter remineralization 
  #R1: CH2O + O2 -> HCO3 + H
  #R2: CH2O + (4/5)NO3 -> HCO3 + (2/5)N2 + (2/5)H2O + (1/5)H
  #R3: CH2O + (1/5)SO4 -> HCO3 + (1/2)HS + (1/2)H

#Secondary redox reactions
  #R4: HS + 2O2 -> SO4 + H
  #R5: NH4 + 2O2 -> NO3 + H2O + 2H
  #R6: FeS + (9/4)O2 + (3/2)H2O -> FeOOH + SO4 + 2H
  #R7: HS + (8/5)NO3 + (3/5)H -> (4/5)N2 + SO4 + (4/5)H2O

#====== PACKAGES ========#
library(ReacTran)

#====== NO DATA INPUT =======#

#====== MODEL FUNCTION ======#
model <- function (t,state,parameters,full.output=FALSE){
  with(as.list(c(parameters)),{
    
    # Initialisation of state variables
    
    CH2O.f <- state[1:N] #N = number of grid cells
    CH2O.s <- state[(1*N+1):(2*N)]
    O2     <- state[(2*N+1):(3*N)]
    SO4    <- state[(3*N+1):(4*N)]
    FeS    <- state[(4*N+1):(5*N)]
    HCO3   <- state[(5*N+1):(6*N)]
    NH4    <- state[(6*N+1):(7*N)]
    NO3    <- state[(7*N+1):(8*N)]
    HS     <- state[(8*N+1):(9*N)]
    
    # Transport terms: Solids [umol cm-3 solid yr-1]
    
    tran.CH2O.f <- tran.1D(C=CH2O.f,flux.up=F.CH2O.f,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)
    tran.CH2O.s <- tran.1D(C=CH2O.s,flux.up=F.CH2O.s,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)
    tran.FeS    <- tran.1D(C=FeS,   flux.up=F.FeS,   v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)
    
    # Transport terms: Solutes [umol cm-3 pore water yr-1]
    
    tran.O2   <- tran.1D(C=O2,  C.up=O2.ow,  v=u.grid,D=D.O2.grid,  VF=por.grid,dx=grid)
    tran.SO4  <- tran.1D(C=SO4, C.up=SO4.ow, v=u.grid,D=D.SO4.grid, VF=por.grid,dx=grid)
    tran.HCO3 <- tran.1D(C=HCO3,C.up=HCO3.ow,v=u.grid,D=D.HCO3.grid,VF=por.grid,dx=grid)
    tran.NH4  <- tran.1D(C=NH4, C.up=NH4.ow, v=u.grid,D=D.NH4.grid, VF=por.grid,dx=grid)
    tran.HS   <- tran.1D(C=HS,  C.up=HS.ow,  v=u.grid,D=D.HS.grid,  VF=por.grid,dx=grid)
    tran.NO3  <- tran.1D(C=NO3, C.up=NO3.ow, v=u.grid,D=D.NO3.grid, VF=por.grid,dx=grid)
    
    # Reaction terms: Expressed per volume of bulk sediment [umol cm-3 yr-1]
    
    Rmin.f <- (1.-por.grid$mid)*k.f*CH2O.f # Organic matter mineralization fast
    Rmin.s <- (1.-por.grid$mid)*k.s*CH2O.s # Organic matter mineralization slow
    Rmin   <- Rmin.f + Rmin.s
    
    O2.lim    <- O2/(O2+K_O2)
    NO3.lim   <- (NO3/(NO3+K_NO3))*(K_O2/(O2+K_O2))
    SO4.lim   <- (SO4/(SO4+K_SO4))*(K_NO3/(NO3+K_NO3))*(K_O2/(O2+K_O2))
    
    f.O2   <- O2.lim/(O2.lim+NO3.lim+SO4.lim)
    f.NO3  <- NO3.lim/(O2.lim+NO3.lim+SO4.lim)
    f.SO4  <- SO4.lim/(O2.lim+NO3.lim+SO4.lim)
    
    R1 <- f.O2*Rmin    # Aerobic respiration
    R2 <- f.NO3*Rmin   # Nitrate reduction
    R3 <- f.SO4*Rmin   # Sulfate reduction 
    
    R4 <- por.grid$mid*k_Sox*O2*HS*(O2>0.0)*(HS>0.0)        # CSO
    R5 <- por.grid$mid*k_NH4*O2*NH4*(O2>0.0)*(NH4>0.0)      # NH4+ oxidation
    R6 <- (1.-por.grid$mid)*k_Sox*O2*FeS*(FeS>0.0)*(O2>0.0) # FeS oxidation
    
    R7 <- por.grid$mid*k_Sni*NO3*HS*(NO3>0.0)*(HS>0.0)      # Sulfide oxidation with nitrate
    
    # Consumption rates
    
    reac.CH2O.f <- -Rmin.f/svf.grid$mid
    reac.CH2O.s <- -Rmin.s/svf.grid$mid
    
    reac.HCO3 <- (R1 + R2 + R3)/(por.grid$mid)
    reac.NH4 <- ((1/CNratio)*(R1 + R2 + R3) - R5)/(por.grid$mid)
    
    reac.O2  <- (- R1 - 2*R4 - 2*R5 - 9/4*R6)/(por.grid$mid)
    reac.NO3 <- (-0.8*R2 + R5 - 8/5*R7)/(por.grid$mid)
    
    reac.FeS  <- (1/2*f.FeS*R3                 - R6)/svf.grid$mid
    reac.SO4  <- (- 1/2*R3           + R4 + R7 + R6)/(por.grid$mid)
    reac.HS   <- (+ 1/2*(1-f.FeS)*R3 - R4 - R7)/(por.grid$mid)
    
    # Irrigation rates
    
    irr.CH2O.f <- 0
    irr.CH2O.s <- 0
    irr.FeS    <- 0
    
    irr.O2   <- irr.grid$mid*(O2.ow-O2)
    irr.SO4  <- irr.grid$mid*(SO4.ow-SO4)
    irr.HCO3 <- irr.grid$mid*(HCO3.ow-HCO3)
    irr.NH4  <- irr.grid$mid*(NH4.ow-NH4)
    irr.NO3  <- irr.grid$mid*(NO3.ow-NO3)
    irr.HS   <- irr.grid$mid*(HS.ow-HS)
    
    # Rate of changes of all species
    
    dCH2O.f <- tran.CH2O.f$dC + reac.CH2O.f
    dCH2O.s <- tran.CH2O.s$dC + reac.CH2O.s
    
    dFeS   <- tran.FeS$dC + reac.FeS
    
    dO2   <- tran.O2$dC   + reac.O2   + irr.O2
    dSO4  <- tran.SO4$dC  + reac.SO4  + irr.SO4
    dHCO3 <- tran.HCO3$dC + reac.HCO3 + irr.HCO3
    dNH4  <- tran.NH4$dC  + reac.NH4  + irr.NH4
    dNO3  <- tran.NO3$dC  + reac.NO3  + irr.NO3
    dHS   <- tran.HS$dC   + reac.HS   + irr.HS
    
    # Burial termsfor CH2O, FeS, and reduced terms
    burial.CH2O.f <- tran.1D(C=CH2O.f,flux.up=F.CH2O.f,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)$flux.down
    burial.CH2O.s <- tran.1D(C=CH2O.s,flux.up=F.CH2O.s,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)$flux.down
    burial.FeS    <- tran.1D(C=FeS,flux.up=F.FeS,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid)$flux.down
    burial.NH4    <- tran.1D(C=NH4, C.up=NH4.ow, v=u.grid,D=D.NH4.grid, VF=por.grid,dx=grid)$flux.down
    burial.HS     <- tran.1D(C=HS,  C.up=HS.ow,  v=u.grid,D=D.HS.grid,  VF=por.grid,dx=grid)$flux.down 
    
    # O2 irrigation flux
    int.irr.O2 <- sum(irr.O2*grid$dx)
    O2.diff.flux.up <- tran.1D(C=O2, C.up=O2.ow, v=u.grid, D=D.O2.grid, VF=por.grid, dx=grid)$flux.up
    
    if (!full.output) return(list(c(dCH2O.f,dCH2O.s,dO2,dSO4,dFeS,dHCO3,dNH4,dNO3,dHS),
                                    'int.irr.O2'=int.irr.O2, 'O2.diff.flux.up'=O2.diff.flux.up, 
                                    'CH2O.f.burial'=burial.CH2O.f, 'CH2O.s.burial'=burial.CH2O.s,
                                    'FeS.burial'=burial.FeS, 'NH4.burial'=burial.NH4, 'HS.burial'=burial.HS)) else
    {
      rate.matrix <- matrix(nrow=N.var,ncol=N.rate)
      dimnames(rate.matrix) <- list(var.names,rate.names)
      
      rate.matrix["CH2O.f" ,] <- check.balance(tran.CH2O.f,reac.CH2O.f,irr.CH2O.f,dCH2O.f,svf.grid)
      rate.matrix["CH2O.s" ,] <- check.balance(tran.CH2O.s,reac.CH2O.s,irr.CH2O.s,dCH2O.s,svf.grid)
      rate.matrix["O2"   ,]   <- check.balance(tran.O2,reac.O2,irr.O2,dO2,por.grid)
      rate.matrix["SO4"  ,]   <- check.balance(tran.SO4,reac.SO4,irr.SO4,dSO4,por.grid)
      rate.matrix["FeS" ,]    <- check.balance(tran.FeS,reac.FeS,irr.FeS,dFeS,svf.grid)
      rate.matrix["HCO3" ,]   <- check.balance(tran.HCO3,reac.HCO3,irr.HCO3,dHCO3,por.grid)
      rate.matrix["NH4"  ,]   <- check.balance(tran.NH4,reac.NH4,irr.NH4,dNH4,por.grid)
      rate.matrix["NO3",]     <- check.balance(tran.NO3,reac.NO3,irr.NO3,dNO3,por.grid)
      rate.matrix["HS",]      <- check.balance(tran.HS,reac.HS,irr.HS,dHS,por.grid)
      
      reac.matrix <- matrix(nrow=N.reac,ncol=1)
      dimnames(reac.matrix) <- list(reac.names,names(reaction.summary))
      
      reac.matrix["R1","rate"] <- sum(R1*grid$dx)
      reac.matrix["R2","rate"] <- sum(R2*grid$dx)
      reac.matrix["R3","rate"] <- sum(R3*grid$dx)
      reac.matrix["R4","rate"] <- sum(R4*grid$dx)
      reac.matrix["R5","rate"] <- sum(R5*grid$dx)
      reac.matrix["R6","rate"] <- sum(R6*grid$dx)
      reac.matrix["R7","rate"] <- sum(R7*grid$dx)
      
      return(list(c(dCH2O.f,dCH2O.s,dO2,dSO4,dFeS,dHCO3,dNH4,dNO3,dHS),
                  R1 = R1,
                  R2 = R2,
                  R3 = R3,
                  R4 = R4,
                  R5 = R5,
                  R6 = R6,
                  R7 = R7,
                  rate.matrix = rate.matrix,
                  reac.matrix = reac.matrix
      ))
    }
  })
}
