#Author: Alison Cribb
#Plotting output to re-create Figure 2

#==== REQUIRED PACKAGES ====#
library(readr)
library(ggplot2)
library(egg)
library(tayloRswift)

#==== LOAD DATA ====#
load('Ediacaran_BioirrigationSensitivityTest_Out.RData')
load('Terreneuvian_BioirrigationSensitivityTest_Out.RData')

#==== PLOTTING OUTPUT ====#
irr_palette <- swift_palettes$taylor1989

Ediacaran_irrtest <- ggplot() +
  geom_vline(xintercept=c(3.65,35.77), linetype='dashed', color=c(irr_palette[3:4])) +
  geom_path(data=Ediacaran_out, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  geom_point(data=Ediacaran_out, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  scale_color_taylor(palette='taylorRed') +
  ggtitle('Ediacaran Db = 0.1') +
  xlab(expression(paste(irr[0]~(yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  theme_classic()

Terreneuvian_irrtest <- ggplot() +
  geom_vline(xintercept=c(3.65,35.77), linetype='dashed', color=c(irr_palette[3:4])) +
  geom_path(data=Terreneuvian_out, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  geom_point(data=Terreneuvian_out, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  scale_color_taylor(palette='taylorRed') +
  ggtitle('Terreneuvian Db = 0.98') +
  xlab(expression(paste(irr[0]~(yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  theme_classic()

ggarrange(Ediacaran_irrtest, Terreneuvian_irrtest, ncol=1)
