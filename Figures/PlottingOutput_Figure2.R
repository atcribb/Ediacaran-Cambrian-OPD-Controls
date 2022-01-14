#Author: Alison Cribb
#Plotting output to re-create Figure 2

#==== REQUIRED PACKAGES ====#
library(readr)
library(ggplot2)
library(egg)
library(tayloRswift)

#===== PROCESS AND PLOT OUTPUT =====#
#Load data: 
load('BiodiffusionSensitivityTestOutput.RData')
View(results_df)

bioturbation_data <- read_csv("~/Desktop/Manucripts/OPD Controls/Submission_ScienceAdvances/Dataset/Db_MLD_foranlysis.csv")

#get biodiffusion coefficients of interest from bioturbation_data
Dbs_all <- as.data.frame(matrix(NA, nrow=4, ncol=3))
colnames(Dbs_all) <- c('Age', 'Source', 'Db')
Dbs_all$Age    <- bioturbation_data[1:4,]$Age
Dbs_all$Source <- bioturbation_data[1:4,]$Source
Dbs_all$Db     <- bioturbation_data[1:4,]$Db_final

Dbs_globalavs <- as.data.frame(matrix(NA, nrow=2, ncol=2))
colnames(Dbs_globalavs) <- c('Age', 'Db')
Dbs_globalavs$Age <- c('Ediacaran', 'Terreneuvian')
Dbs_globalavs$Db  <- bioturbation_data[7:8,]$Db_final

Dbpalette <- swift_palettes$taylor1989

Db_SA_out <- ggplot() + 
  #geom_vline(data=Dbs_globalavs, aes(xintercept=Db), color=Dbpalette[3:4], linetype='longdash', size=1.2) +
  geom_vline(xintercept=c(0.1,0.98), color=Dbpalette[3:4], linetype='longdash', size=1.2) +
  geom_path(data=results_df, aes(x=Db, y=OPD, color=as.factor(O2_bw))) + #plot points 
  geom_point(data=results_df, aes(x=Db, y=OPD, color=as.factor(O2_bw))) + #plot path 
  scale_color_taylor(palette='taylorRed') +
  xlab(expression(paste(D[B]~(cm^{2}~yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  theme_classic()     
Db_SA_out  
