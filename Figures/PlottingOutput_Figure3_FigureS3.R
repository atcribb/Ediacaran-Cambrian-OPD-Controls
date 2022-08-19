#Author: Alison Cribb
#plotting output script for Corg, O2bw, and oxygen penetration depth sensitivty anlaysis
#Figure 3 for different bioturbation intensities and Figure S2 for Biomixing vs Bioirrigation effects

#==== PACAKGES ====#
library(tidyverse)
library(ggmap)
library(ggplot2)
library(egg)
library(MASS)
library(metR)
library(paletteer)

#===== LOAD DATA =====#
#Uncomment whatever the data is that needs loading in 

# load('Sim1_NoBioturbation.RData')
# load('Sim2_EdiacaranBiomixing.RData')
# load('Sim3_EdiacaranBioirrigation.RData')
# load('Sim4_EdiacaranBioturbation.RData')
# load('Sim5_TerreneuvianBiomixing.RData')
# load('Sim6_TerreneuvianBioirrigation.RData')
# load('Sim7_TerreneuvianBioturbation.RData')
# load('Sim8_ModernBiomixing.RData')
# load('Sim9_ModernBioirrigation.RData')
# load('Sim10_ModernBioturbation.RData')

#take subset of data to plot Corg=100-200 [umol cm-2 yr-1] and O2.ow=0.014-0.14 mM
results.df$Corg_total <- results.df$Corg_fast + results.df$Corg_slow
corgsubset <- results.df[which(results.df$Corg_total>100 & results.df$Corg_total<450),]
finalsubset <- corgsubset[which(corgsubset$O2_bw>0.014 & corgsubset$O2_bw<0.28),]

#use geom_contour_fill -- change break sequence in geom_contour_fill if necessary 
label_names <- c('0 - 0.5 cm',
                 '0.5 - 1 cm',
                 '1 - 1.5 cm',
                 '1.5 - 2 cm',
                 '2 - 2.5 cm', 
                 '2.5 - 3 cm',
                 '3 - 3.5 cm', 
                 '3.5 - 4 cm', 
                 '4 - 4.5 cm', 
                 '4.5 - 5 cm', 
                 '5 - 5.5 cm', 
                 '5.5 - 6 cm',
                 '6 - 6.5 cm',
                 '6.5 - 7 cm',
                 '7 - 7.5 cm',
                 '7.5 - 8 cm',
                 '8 - 8.5 cm',
                 '8.5 - 9 cm',
                 '9 - 9.5 cm',
                 '9.5 - 10 cm',
                 '>10 cm')

contourplot <- ggplot(data=finalsubset) +
  geom_contour_fill(aes(O2_bw, Corg_total, z=O2_depth, fill=stat(level)), finalsubset, na.fill=TRUE, breaks=c(seq(0, 10, 0.5), 15)) +
  xlab(expression(paste(Bottom~Water~O[2]~(mM)))) +
  ylab(expression(paste(Organic~Matter~Flux~(mu~mol~cm^{2}~yr^{'-1'})))) +
  labs(fill='Oxygen Penetration Depth') +
  scale_fill_viridis_d(option='mako', labels=label_names) +
  scale_y_continuous(breaks=seq(from=100, to=450, by=50)) +
  scale_x_continuous(breaks=seq(from=0.014, to=0.28, length.out=5), labels=scales::number_format(accuracy=0.01)) +
  coord_cartesian(expand=c(0,1)) +
  # ggtitle('No bioturbation') +                  #change to appropriate title
  # ggtitle('Ediacaran biomixing') +
  # ggtitle('Ediacaran bioirrigation') +
  # ggtitle('Ediacaran bioturbation') +
  # ggtitle('Terreneuvian biomixing') +
  # ggtitle('Terreneuvian bioirrigation') +
  # ggtitle('Terreneuvian bioturbation') +
  # ggtitle('Modern biomixing') +
  # ggtitle('Modern bioirrigation') +
  # ggtitle('Modern bioturbation') +
  theme_classic() +
  theme(
    #text=element_text(family='Avenir', color='black'),
    plot.title=element_text(size=18, hjust=0.5)
  )

#contourplot

#Uncomment appropriate ones to arrange into a final plot
#Figure 3 contour plots
#nobioturbation <- contourplot 
#ediacaranboth <- contourplot 
#terreneuvianboth <- contourplot 
#modernboth <- contourplot

#Other contour plots in Figure S2
#ediacaranbiomix <- contourplot
#terreneuvianbiomix <- contourplot
#modernbiomix <- contourplot 
#ediacaranbioirr <- contourplot
#terreneuvianbioirr <- contourplot
#modernbioirr <- contourplot



Figure3 <- ggarrange(nobioturbation, ediacaranboth, terreneuvianboth, modernboth,
          ncol=2, nrow=2)
ggsave('biomix_bioirr_fandiagrams_raw_120621.pdf', all, width=16, height=18, units='in')

FigureS2 <- ggarrange(ediacaranbiomix, ediacaranbioirr, ediacaranboth,
                 terreneuvianbiomix, terreneuvianbioirr, terreneuvianboth,
                 modernbiomix, modernbioirr, modernboth,
                 ncol=3, nrow=3)
ggsave('allbio_fandiagrams_raw_120621.pdf', all.all, width=24, height=18, units='in')


#Save as appropriate file name: examples --
#ggsave('NoBioturbation.pdf', plot=contourplot, width=7.33, height=5.17, units='in')
#ggsave('EdiacaranBioturbation.pdf', plot=contourplot, width=7.33, height=5.17, units='in')
#ggsave('TerreneuvianBioturbation.pdf', plot=contourplot, width=7.33, height=5.17, units='in')
#ggsave('ModernBioturbation.pdf', plot=contourplot, width=7.33, height=5.17, units='in')

