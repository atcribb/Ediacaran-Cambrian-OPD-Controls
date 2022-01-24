#Author: Alison Cribb
#Plotting output to re-create Figure 2

#==== REQUIRED PACKAGES ====#
library(readr)
library(ggplot2)
library(egg)
library(tayloRswift)

#===== PROCESS AND PLOT OUTPUT =====#
#Load data: 
load('BiodiffusionSensitivityTestOutput_low.RData')
low_results <- results_df
load('BiodiffusionSensitivityTestOutput_med.RData')
med_results <- results_df
load('BiodiffusionSensitivityTestOutput_high.RData')
high_results <- results_df


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

Db_SA_out_low <- ggplot() + 
  geom_vline(xintercept=c(0.1,0.98), color=Dbpalette[3:4], linetype='longdash', size=1.2) +
  geom_path(data=low_results, aes(x=Db, y=OPD, color=as.factor(O2_bw))) + #plot points 
  geom_point(data=low_results, aes(x=Db, y=OPD, color=as.factor(O2_bw))) + #plot path 
  scale_color_taylor(palette='taylorRed') +
  xlab(expression(paste(D[B]~(cm^{2}~yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  ggtitle(expression(paste('Organic matter flux = 150', ~mu, mol~cm^{-2}, yr^{-1}))) +
  theme_classic()     


Db_SA_out_med <- ggplot() + 
  geom_vline(xintercept=c(0.1,0.98), color=Dbpalette[3:4], linetype='longdash', size=1.2) +
  geom_path(data=med_results, aes(x=Db, y=OPD, color=as.factor(O2_bw))) + #plot points 
  geom_point(data=med_results, aes(x=Db, y=OPD, color=as.factor(O2_bw))) + #plot path 
  scale_color_taylor(palette='taylorRed') +
  xlab(expression(paste(D[B]~(cm^{2}~yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  ggtitle(expression(paste('Organic matter flux = 400', ~mu, mol~cm^{-2}, yr^{-1}))) +
  theme_classic()     

Db_SA_out_high <- ggplot() + 
  geom_vline(xintercept=c(0.1,0.98), color=Dbpalette[3:4], linetype='longdash', size=1.2) +
  geom_path(data=high_results, aes(x=Db, y=OPD, color=as.factor(O2_bw))) + #plot points 
  geom_point(data=high_results, aes(x=Db, y=OPD, color=as.factor(O2_bw))) + #plot path 
  scale_color_taylor(palette='taylorRed') +
  xlab(expression(paste(D[B]~(cm^{2}~yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  ggtitle(expression(paste('Organic matter flux = 700', ~mu, mol~cm^{-2}, yr^{-1}))) +
  theme_classic()     


ggarrange(Db_SA_out_low, Db_SA_out_med, Db_SA_out_high, ncol=1)
