#Author: Alison Cribb
#Plotting output to re-create Figure 2

#==== REQUIRED PACKAGES ====#
library(readr)
library(ggplot2)
library(egg)
library(tayloRswift)

#==== LOAD DATA ====#
#Ediacaran Db
load('Ediacaran_BioirrigationSensitivityTest_low.RData')
Ediacaran_low <- Ediacaran_out
load('Ediacaran_BioirrigationSensitivityTest_med.RData')
Ediacaran_med <- Ediacaran_out
load('Ediacaran_BioirrigationSensitivityTest_high.RData')
Ediacaran_high <- Ediacaran_out

#Terreneuvian Db
load('Terreneuvian_BioirrigationSensitivityTest_low.RData')
Terreneuvian_low <- Terreneuvian_out
load('Terreneuvian_BioirrigationSensitivityTest_med.RData')
Terreneuvian_med <- Terreneuvian_out
load('Terreneuvian_BioirrigationSensitivityTest_high.RData')
Terreneuvian_high <- Terreneuvian_out

#==== PLOTTING OUTPUT ====#
irr_palette <- swift_palettes$taylor1989

Ediacaran_irrtest_low <- ggplot() +
  geom_vline(xintercept=c(3.65,35.77), linetype='dashed', color=c(irr_palette[3:4])) +
  geom_path(data=Ediacaran_low, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  geom_point(data=Ediacaran_low, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  scale_color_taylor(palette='taylorRed') +
  ggtitle(expression(paste('Organic matter flux = 150', ~mu, mol~cm^{-2}, yr^{-1}))) +
  xlab(expression(paste(irr[0]~(yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  theme_classic() +
  theme(
    legend.position='none'
  )


Ediacaran_irrtest_med <- ggplot() +
  geom_vline(xintercept=c(3.65,35.77), linetype='dashed', color=c(irr_palette[3:4])) +
  geom_path(data=Ediacaran_med, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  geom_point(data=Ediacaran_med, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  scale_color_taylor(palette='taylorRed') +
  ggtitle(expression(paste('Organic matter flux = 400', ~mu, mol~cm^{-2}, yr^{-1}))) +
  xlab(expression(paste(irr[0]~(yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  theme_classic() +
  theme(
    legend.position='none'
  )


Ediacaran_irrtest_high <- ggplot() +
  geom_vline(xintercept=c(3.65,35.77), linetype='dashed', color=c(irr_palette[3:4])) +
  geom_path(data=Ediacaran_high, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  geom_point(data=Ediacaran_high, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  scale_color_taylor(palette='taylorRed') +
  ggtitle(expression(paste('Organic matter flux = 700', ~mu, mol~cm^{-2}, yr^{-1}))) +
  xlab(expression(paste(irr[0]~(yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  theme_classic() +
  theme(
    legend.position='none'
  )


Terreneuvian_irrtest_low <- ggplot() +
  geom_vline(xintercept=c(3.65,35.77), linetype='dashed', color=c(irr_palette[3:4])) +
  geom_path(data=Terreneuvian_low, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  geom_point(data=Terreneuvian_low, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  scale_color_taylor(palette='taylorRed') +
  ggtitle(expression(paste('Organic matter flux = 150', ~mu, mol~cm^{-2}, yr^{-1}))) +
  xlab(expression(paste(irr[0]~(yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  theme_classic() +
  theme(
    legend.position='none'
  )


Terreneuvian_irrtest_med <- ggplot() +
  geom_vline(xintercept=c(3.65,35.77), linetype='dashed', color=c(irr_palette[3:4])) +
  geom_path(data=Terreneuvian_med, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  geom_point(data=Terreneuvian_med, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  scale_color_taylor(palette='taylorRed') +
  ggtitle(expression(paste('Organic matter flux = 400', ~mu, mol~cm^{-2}, yr^{-1}))) +
  xlab(expression(paste(irr[0]~(yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  theme_classic() +
  theme(
    legend.position='none'
  )


Terreneuvian_irrtest_high <- ggplot() +
  geom_vline(xintercept=c(3.65,35.77), linetype='dashed', color=c(irr_palette[3:4])) +
  geom_path(data=Terreneuvian_high, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  geom_point(data=Terreneuvian_high, aes(x=irr, y=OPD, color=as.factor(O2_bw))) +
  scale_color_taylor(palette='taylorRed') +
  ggtitle(expression(paste('Organic matter flux = 700', ~mu, mol~cm^{-2}, yr^{-1}))) +
  xlab(expression(paste(irr[0]~(yr^{-1})))) +
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(15,0)) +
  theme_classic() +
  theme(
    legend.position='none'
  )

ggarrange(Ediacaran_irrtest_low, Ediacaran_irrtest_med, Ediacaran_irrtest_high,
          Terreneuvian_irrtest_low, Terreneuvian_irrtest_med, Terreneuvian_irrtest_high,
          ncol=3, nrow=2)
