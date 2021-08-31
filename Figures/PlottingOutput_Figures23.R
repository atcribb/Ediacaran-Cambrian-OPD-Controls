#plotting output script for Corg, O2bw, and oxygen penetration depth sensitivty anlaysis


#==== PACAKGES ====#
library(tidyverse)
library(ggmap)
library(ggplot2)
library(egg)
library(MASS)
library(metR)


#===== LOAD DATA =====#
#Uncomment whatever the data is that needs loading in 

 load('Sim1_NoBioturbation.RData')
# load('Sim2_WeakBiomixing.RData')
# load('Sim3_WeakBioirrigation.RData')
# load('Sim4_WeakBioturbation.RData')
# load('Sim5_StrongBiomixing.RData')
# load('Sim6_StrongBioirrigation.RData')
# load('Sim7_StrongBioturbation.RData')

#====SET UP PLOTTING OUTPUT====#
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
  ylab(expression(paste(Organic~Matter~Flux~(mu~mol~g^{'-1'})))) +
  labs(fill='Oxygen Penetration Depth') +
  scale_fill_viridis_d(labels=label_names) +
  scale_y_continuous(breaks=seq(from=100, to=450, by=50)) +
  scale_x_continuous(breaks=seq(from=0.014, to=0.28, length.out=5), labels=scales::number_format(accuracy=0.01)) +
  coord_cartesian(expand=c(0,1)) +
  ggtitle('No bioturbation') +                  #change to appropriate title
  #ggtitle('Weak Biomixing (Db = 1.0)') +
  #ggtitle('Weak bioirrigation (a0=36.5)') +
  #ggtitle('Weak biomixing and bioirrigation (Db=1, a0=36.5)') +
  #ggtitle('Strong Biomixing (Db = 10)') +
  #ggtitle('Strong bioirrigation (a0=365)') +
  #ggtitle('Strong biomixing and bioirrigation (Db=10, a0=365)') +
  theme_classic() +
  theme(
    text=element_text(family='Avenir', color='black'),
    plot.title=element_text(family='Avenir Heavy', size=18, hjust=0.5)
  )

contourplot

#save with appropriate simuation name!
