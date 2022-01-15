#Author: Alison Cribb

#Code to create figure 1 - relative occurrences of biomixing vs bioirrigation trace fossils
#Dataset primarily from Buatois et al. (2020) Science Advances, with additions noted in dataset .csv file
#Some trace fossils have been removed based on functional group -- see manuscript text for details 

#==== REQUIRED PACKAGES ====#
library(readr)
library(ggplot2)
library(egg)
library(tayloRswift)

#==== DATA INPUT ====#
ichnodata <- read_csv('Ichnotaxa_dataset.csv') #load the ichnotaxa_dataset.csv dataet 

#Break data up into Ediacaran and Terrneuvian time bins
#ages according to Buatois et al. (2020) or otherwise noted sources
Ediacaran_ichnodata <- ichnodata[which(ichnodata$Period=='Ediacaran'),]
Cambrian_ichnodata <- ichnodata[which(ichnodata$Period=='Terreneuvian'),]

#==== ANALYSIS -> occurrences in stratigraphic units  ====#
#calculate relative abundance of biomixing vs. bioirrigation in the Ediacaran vs. Cambrian
#dataframe set up
relabund_df <- as.data.frame(matrix(NA, nrow=4, ncol=5))
colnames(relabund_df) <- c('Age', 'Behavior', 'Strat_occurrences', 'Total_ichnotaxa', 'Percent')
rownames(relabund_df) <- c('Ediacaran_biomix', 'Ediacaran_bioirr', 'Terreneuvian_biomix', 'Terreneuvian_bioirr')
relabund_df$Age <- c('Ediacaran', 'Ediacaran', 'Terreneuvian', 'Terreneuvian')
relabund_df$Behavior <- rep(c('Biomixing', 'Bioirrigation'),2)
relabund_df

#Fill in dataframe
#total ichnotaxa stratigraphic occurrences in each time bin
tot_ediacaran <- length(Ediacaran_ichnodata$Ichnogenera)
tot_cambrian  <- length(Cambrian_ichnodata$Ichnogenera)

#total biomixing stratigraphic occurrences in each time bin 
biomix_ediacaran <- nrow(Ediacaran_ichnodata[which(Ediacaran_ichnodata$Bioturbation_behavior=='Biomixing'),])
biomix_cambrian  <- nrow(Cambrian_ichnodata[which(Cambrian_ichnodata$Bioturbation_behavior=='Biomixing'),])

#total bioirrigation stratigraphic occurrences in each time bin 
bioirr_ediacaran <- nrow(Ediacaran_ichnodata[which(Ediacaran_ichnodata$Bioturbation_behavior=='Bioirrigation'),])
bioirr_cambrian  <- nrow(Cambrian_ichnodata[which(Cambrian_ichnodata$Bioturbation_behavior=='Bioirrigation'),])

#add to dataframe
relabund_df['Ediacaran_biomix', 'Strat_occurrences'] <- biomix_ediacaran
relabund_df['Ediacaran_bioirr', 'Strat_occurrences'] <- bioirr_ediacaran
relabund_df['Ediacaran_biomix', 'Total_ichnotaxa'] <- tot_ediacaran
relabund_df['Ediacaran_bioirr', 'Total_ichnotaxa'] <- tot_ediacaran
relabund_df['Ediacaran_biomix', 'Percent'] <- (biomix_ediacaran/tot_ediacaran)*100
relabund_df['Ediacaran_bioirr', 'Percent'] <- (bioirr_ediacaran/tot_ediacaran)*100

relabund_df['Terreneuvian_biomix', 'Strat_occurrences'] <- biomix_cambrian
relabund_df['Terreneuvian_bioirr', 'Strat_occurrences'] <- bioirr_cambrian
relabund_df['Terreneuvian_biomix', 'Total_ichnotaxa'] <- tot_cambrian
relabund_df['Terreneuvian_bioirr', 'Total_ichnotaxa'] <- tot_cambrian
relabund_df['Terreneuvian_biomix', 'Percent'] <- (biomix_cambrian/tot_cambrian)*100
relabund_df['Terreneuvian_bioirr', 'Percent'] <- (bioirr_cambrian/tot_cambrian)*100

relabund_df

#==== PLOTTING OUTPUT ====#
strat_relabund <- ggplot(data=relabund_df, aes(fill=Behavior, y=Percent, x=Age)) +
  geom_bar(position='stack', stat='identity') +
  scale_fill_taylor(palette='taylor1989', reverse=TRUE) +
  ylab('Relative Abudance (%)') +
  ggtitle('From unique stratigraphic occurrences') +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=12, color='black'),
    axis.title.y=element_text(size=16, color='black'),
    axis.text.y=element_text(size=12, color='black'),
    legend.key.size=unit(1, 'cm'),
    legend.title=element_text(size=14),
    legend.text=element_text(size=12))
  
strat_relabund

#==== trace fossil diversity for biomixing vs. bioirrgation  ====#
#dataframe set up
behaviordiv_df <- as.data.frame(matrix(NA, nrow=4, ncol=5))
colnames(behaviordiv_df) <- c('Age', 'Behavior', 'No_ichnogenera', 'Total_ichnogenera', 'Percent')
rownames(behaviordiv_df) <- c('Ediacaran_biomix', 'Ediacaran_bioirr', 'Terreneuvian_biomix', 'Terreneuvian_bioirr')
behaviordiv_df$Age <- c('Ediacaran', 'Ediacaran', 'Terreneuvian', 'Terreneuvian')
behaviordiv_df$Behavior <- rep(c('Biomixing', 'Bioirrigation'),2)
behaviordiv_df

#analysis
#total number of unique ichnogenera in each time bin
tot_ediacaran <- length(unique(Ediacaran_ichnodata$Ichnogenera))
tot_cambrian <- length(unique(Cambrian_ichnodata$Ichnogenera))

#number of unique biomixing trace fossils
all_ediacaranbiomixers <- Ediacaran_ichnodata[which(Ediacaran_ichnodata$Bioturbation_behavior=='Biomixing'),]
ediacaran_biomixers <- length(unique(all_ediacaranbiomixers$Ichnogenera))
all_ediacaranbioirrigators <- Ediacaran_ichnodata[which(Ediacaran_ichnodata$Bioturbation_behavior=='Bioirrigation'),]
ediacaran_bioirrigators <- length(unique(all_ediacaranbioirrigators$Ichnogenera))

all_cambrianbiomixers <- Cambrian_ichnodata[which(Cambrian_ichnodata$Bioturbation_behavior=='Biomixing'),]
cambrian_biomixers <- length(unique(all_cambrianbiomixers$Ichnogenera))
all_cambrianbioirrigators <- Cambrian_ichnodata[which(Cambrian_ichnodata$Bioturbation_behavior=='Bioirrigation'),]
cambrian_bioirrigators <- length(unique(all_cambrianbioirrigators$Ichnogenera))

behaviordiv_df['Ediacaran_biomix', 'Total_ichnogenera'] <- tot_ediacaran
behaviordiv_df['Ediacaran_bioirr', 'Total_ichnogenera'] <- tot_ediacaran
behaviordiv_df['Ediacaran_biomix', 'No_ichnogenera'] <- ediacaran_biomixers
behaviordiv_df['Ediacaran_bioirr', 'No_ichnogenera'] <- ediacaran_bioirrigators
behaviordiv_df['Ediacaran_biomix', 'Percent'] <- (ediacaran_biomixers/tot_ediacaran)*100
behaviordiv_df['Ediacaran_bioirr', 'Percent'] <- (ediacaran_bioirrigators/tot_ediacaran)*100

behaviordiv_df['Terreneuvian_biomix', 'Total_ichnogenera'] <- tot_cambrian
behaviordiv_df['Terreneuvian_bioirr', 'Total_ichnogenera'] <- tot_cambrian
behaviordiv_df['Terreneuvian_biomix', 'No_ichnogenera'] <- cambrian_biomixers
behaviordiv_df['Terreneuvian_bioirr', 'No_ichnogenera'] <- cambrian_bioirrigators
behaviordiv_df['Terreneuvian_biomix', 'Percent'] <- (cambrian_biomixers/tot_cambrian)*100
behaviordiv_df['Terreneuvian_bioirr', 'Percent'] <- (cambrian_bioirrigators/tot_cambrian)*100

#==== PLOTTING OUTPUT ====#
behavior_relabund <- ggplot(data=behaviordiv_df, aes(fill=Behavior, y=Percent, x=Age)) +
  geom_bar(position='stack', stat='identity') +
  scale_fill_taylor(palette='taylor1989', reverse=TRUE) +
  ylab('Relative Abundance (%)') +
  ggtitle('From total unique trace fossils') +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=12, color='black'),
    axis.title.y=element_text(size=16, color='black'),
    axis.text.y=element_text(size=12, color='black'),
    legend.position='none'
    )
behavior_relabund

ggarrange(behavior_relabund, strat_relabund, ncol=2)

