#Author: Alison Cribb
#Plotting output to re-create figure to O2 fluxes in supplemental 

library(ggplot2)
library(egg)
library(tidyr)
library(RColorBrewer)

# Intput data
load('O2irrflux_lowO2_081022.RData')
load('O2irrflux_highO2_081022.RData')

# Plotting output -- O2flux vs bioturbation cases, 6 plots (3 CH2O, 2 O2)
# Extract data
lowO2_lowOC <- subset(lowO2_irr_results, lowO2_irr_results$F.CH2O==150)
lowO2_lowOC <- pivot_longer(lowO2_lowOC, cols=6:7, values_to='flux', names_to='flux_source')
lowO2_medOC <- subset(lowO2_irr_results, lowO2_irr_results$F.CH2O==300)
lowO2_medOC <- pivot_longer(lowO2_medOC, cols=6:7, values_to='flux', names_to='flux_source')
lowO2_highOC <- subset(lowO2_irr_results, lowO2_irr_results$F.CH2O==450)
lowO2_highOC <- pivot_longer(lowO2_highOC, cols=6:7, values_to='flux', names_to='flux_source')
#lowO2_highOC <- as.data.frame(lowO2_highOC)

highO2_lowOC <- subset(highO2_irr_results, highO2_irr_results$F.CH2O==150)
highO2_lowOC <- pivot_longer(highO2_lowOC, cols=6:7, values_to='flux', names_to='flux_source')
highO2_medOC <- subset(highO2_irr_results, highO2_irr_results$F.CH2O==300)
highO2_medOC <- pivot_longer(highO2_medOC, cols=6:7, values_to='flux', names_to='flux_source')
highO2_highOC <- subset(highO2_irr_results, highO2_irr_results$F.CH2O==450)
highO2_highOC <- pivot_longer(highO2_highOC, cols=6:7, values_to='flux', names_to='flux_source')

# Plots
#  low O2 
#PLOT 1 - low OC
lowO2_lowOC$Bioturbation <- factor(lowO2_lowOC$Bioturbation, levels=c('None', 'Ediacaran', 'Terreneuvian', 'Modern'))
lowO2_lowOC$flux_source <- factor(lowO2_lowOC$flux_source, levels=c('irr.O2', 'diff.O2'))
plot1 <- ggplot(data=lowO2_lowOC) +
  geom_bar(aes(fill=flux_source, y=flux, x=Bioturbation), position='stack', stat='identity') +
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(limits=c(0, 450), expand=c(0,0)) +
  ggtitle('O2=0.07 mM, CH2O=150 umol cm-2 yr-1') +
  ylab('O2 flux') +
  xlab('Bioturbation Intensity') +
  theme_bw() +
  theme(legend.position='none')

#PLOT 2 - med OC
lowO2_medOC$Bioturbation <- factor(lowO2_medOC$Bioturbation, levels=c('None', 'Ediacaran', 'Terreneuvian', 'Modern'))
lowO2_medOC$flux_source <- factor(lowO2_medOC$flux_source, levels=c('irr.O2', 'diff.O2'))
plot2 <- ggplot(data=lowO2_medOC) +
  geom_bar(aes(fill=flux_source, y=flux, x=Bioturbation), position='stack', stat='identity') +
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(limits=c(0, 450), expand=c(0,0)) +
  ggtitle('O2=0.07 mM, CH2O=300 umol cm-2 yr-1') +
  ylab('O2 flux') +
  xlab('Bioturbation Intensity') +
  theme_bw() +
  theme(legend.position='none')

#PLOT 3 - high OC
lowO2_highOC$Bioturbation <- factor(lowO2_highOC$Bioturbation, levels=c('None', 'Ediacaran', 'Terreneuvian', 'Modern'))
lowO2_highOC$flux_source <- factor(lowO2_highOC$flux_source, levels=c('irr.O2', 'diff.O2'))
plot3 <- ggplot(data=lowO2_highOC) +
  geom_bar(aes(fill=flux_source, y=flux, x=Bioturbation), position='stack', stat='identity') +
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(limits=c(0, 450), expand=c(0,0)) +
  ggtitle('O2=0.07 mM, CH2O=450 umol cm-2 yr-1') +
  ylab('O2 flux') +
  xlab('Bioturbation Intensity') +
  theme_bw() +
  theme(legend.position='none')

#High O2s
#PLOT 4 - low OC
highO2_lowOC$Bioturbation <- factor(highO2_lowOC$Bioturbation, levels=c('None', 'Ediacaran', 'Terreneuvian', 'Modern'))
highO2_lowOC$flux_source <- factor(highO2_lowOC$flux_source, levels=c('irr.O2', 'diff.O2'))
plot4 <- ggplot(data=highO2_lowOC) +
  geom_bar(aes(fill=flux_source, y=flux, x=Bioturbation), position='stack', stat='identity') +
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(limits=c(0, 450), expand=c(0,0)) +
  ggtitle('O2=0.14 mM, CH2O=150 umol cm-2 yr-1') +
  ylab('O2 flux') +
  xlab('Bioturbation Intensity') +
  theme_bw()

#PLOT 5 - med OC
highO2_medOC$Bioturbation <- factor(highO2_medOC$Bioturbation, levels=c('None', 'Ediacaran', 'Terreneuvian', 'Modern'))
highO2_medOC$flux_source <- factor(highO2_medOC$flux_source, levels=c('irr.O2', 'diff.O2'))
plot5 <- ggplot(data=highO2_medOC) +
  geom_bar(aes(fill=flux_source, y=flux, x=Bioturbation), position='stack', stat='identity') +
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(limits=c(0, 450), expand=c(0,0)) +
  ggtitle('O2=0.14 mM, CH2O=300 umol cm-2 yr-1') +
  ylab('O2 flux') +
  xlab('Bioturbation Intensity') +
  theme_bw() +
  theme(legend.position='none')


#PLOT 6 - high OC
highO2_highOC$Bioturbation <- factor(highO2_highOC$Bioturbation, levels=c('None', 'Ediacaran', 'Terreneuvian', 'Modern'))
highO2_highOC$flux_source <- factor(highO2_highOC$flux_source, levels=c('irr.O2', 'diff.O2'))
plot6 <- ggplot(data=highO2_highOC) +
  geom_bar(aes(fill=flux_source, y=flux, x=Bioturbation), position='stack', stat='identity') +
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(limits=c(0, 450), expand=c(0,0)) +
  ggtitle('O2=0.14 mM, CH2O=450 umol cm-2 yr-1') +
  ylab('O2 flux') +
  xlab('Bioturbation Intensity') +
  theme_bw() +
  theme(legend.position='none')

O2fluxplots <- ggarrange(plot1, plot4,
          plot2, plot5,
          plot3, plot6,
          ncol=2, nrow=3)

ggsave(plot=O2fluxplots, 'O2fluxplots.pdf', width=10, height=12)
