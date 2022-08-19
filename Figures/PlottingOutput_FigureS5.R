#Author: Alison Cribb
#Plotting output to re-create burial rates figures in supplemental
#****This script will plot burial of CH2O, FeS, which are in the paper, as well as HS and NH4, to include all reduced compounds

#load data
load('terreneuvian_burial_df.RData')
terr_results <- results_df
load('ediacaran_burial_df.RData')
ediacaran_results <- results_df
load('nobio_burial_df.RData')
nobio_results <- results_df

full_results <- rbind(nobio_results, ediacaran_results, terr_results)

CH2O_low <- subset(full_results, full_results$CH2O_tot==150)
CH2O_mid <- subset(full_results, full_results$CH2O_tot==300)
CH2O_high <- subset(full_results, full_results$CH2O_tot==450)
CH2O_veryhigh <- subset(full_results, full_results$CH2O_tot==700)

#=== plotting output =====#
#low CH2O plots
low_bCH2O <- ggplot(data=CH2O_low) +
  geom_line(aes(x=O2_bw, y=(burial_CH2O.s + burial_CH2O.f), color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=(burial_CH2O.s + burial_CH2O.f), color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.6)) +
  scale_color_brewer(palette='Blues') +
  ylab('CH2O(tot) burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

low_bFeS <- ggplot(data=CH2O_low) +
  geom_line(aes(x=O2_bw, y=burial_FeS, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_FeS, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,16)) +
  scale_color_brewer(palette='Blues') +
  ylab('FeS burial (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')


low_bHS <- ggplot(data=CH2O_low) +
  geom_line(aes(x=O2_bw, y=burial_HS, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_HS, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.20)) +
  scale_color_brewer(palette='Blues') +
  ylab('HS burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')


low_bNH4 <- ggplot(data=CH2O_low) +
  geom_line(aes(x=O2_bw, y=burial_NH4, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_NH4, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.06)) +
  scale_color_brewer(palette='Blues') +
  ylab('NH4 burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw()

#mid CH2O plots
mid_bCH2O <- ggplot(data=CH2O_mid) +
  geom_line(aes(x=O2_bw, y=(burial_CH2O.s + burial_CH2O.f), color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=(burial_CH2O.s + burial_CH2O.f), color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.6)) +
  scale_color_brewer(palette='Blues') +
  ylab('CH2O(tot) burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

mid_bFeS <- ggplot(data=CH2O_mid) +
  geom_line(aes(x=O2_bw, y=burial_FeS, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_FeS, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,16)) +
  scale_color_brewer(palette='Blues') +
  ylab('FeS burial (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

mid_bHS <- ggplot(data=CH2O_mid) +
  geom_line(aes(x=O2_bw, y=burial_HS, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_HS, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.20)) +
  scale_color_brewer(palette='Blues') +
  ylab('HS burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

mid_bNH4 <- ggplot(data=CH2O_mid) +
  geom_line(aes(x=O2_bw, y=burial_NH4, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_NH4, color=as.factor(Db))) +
  scale_color_brewer(palette='Blues') +
  scale_y_continuous(limits=c(0,0.06)) +
  ylab('NH4 burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

#high CH2O plots
high_bCH2O <- ggplot(data=CH2O_high) +
  geom_line(aes(x=O2_bw, y=(burial_CH2O.s + burial_CH2O.f), color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=(burial_CH2O.s + burial_CH2O.f), color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.6)) +
  scale_color_brewer(palette='Blues') +
  ylab('CH2O(tot) burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

high_bFeS <- ggplot(data=CH2O_high) +
  geom_line(aes(x=O2_bw, y=burial_FeS, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_FeS, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,16)) +
  scale_color_brewer(palette='Blues') +
  ylab('FeS burial (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

high_bHS <- ggplot(data=CH2O_high) +
  geom_line(aes(x=O2_bw, y=burial_HS, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_HS, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.20)) +
  scale_color_brewer(palette='Blues') +
  ylab('HS burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

high_bNH4 <- ggplot(data=CH2O_high) +
  geom_line(aes(x=O2_bw, y=burial_NH4, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_NH4, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.06)) +
  scale_color_brewer(palette='Blues') +
  ylab('NH4 burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

#very high CH2O plots
vhigh_bCH2O <- ggplot(data=CH2O_veryhigh) +
  geom_line(aes(x=O2_bw, y=(burial_CH2O.s + burial_CH2O.f), color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=(burial_CH2O.s + burial_CH2O.f), color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.6)) +
  scale_color_brewer(palette='Blues') +
  ylab('CH2O(tot) burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

vhigh_bFeS <- ggplot(data=CH2O_veryhigh) +
  geom_line(aes(x=O2_bw, y=burial_FeS, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_FeS, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,16)) +
  scale_color_brewer(palette='Blues') +
  ylab('FeS burial (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

vhigh_bHS <- ggplot(data=CH2O_veryhigh) +
  geom_line(aes(x=O2_bw, y=burial_HS, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_HS, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.20)) +
  scale_color_brewer(palette='Blues') +
  ylab('HS burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

vhigh_bNH4 <- ggplot(data=CH2O_veryhigh) +
  geom_line(aes(x=O2_bw, y=burial_NH4, color=as.factor(Db))) +
  geom_point(aes(x=O2_bw, y=burial_NH4, color=as.factor(Db))) +
  scale_y_continuous(limits=c(0,0.06)) +
  scale_color_brewer(palette='Blues') +
  ylab('NH4 burial  (umol cm-2 yr-1)') +
  xlab('O2_bw (mM)') +
  theme_bw() +
  theme(legend.position='none')

#full plot
ggarrange(low_bCH2O, low_bFeS, low_bHS, low_bNH4,
          mid_bCH2O, mid_bFeS, mid_bHS, mid_bNH4,
          high_bCH2O, high_bFeS, high_bHS, high_bNH4,
          vhigh_bCH2O, vhigh_bFeS, vhigh_bHS, vhigh_bNH4,
          ncol=4, nrow=4
          )


