#Plotting output to recreate oxygen consumption profiles in SM6

#==== REQUIRED PACKAGES ====#
library(ggplot2)
library(egg)

#====LOAD DATA====#
load('SM6_nobioturbation.RData')
load('SM6_weakbiomixing.RData')
load('SM6_weakbioirrigation.RData')
load('SM6_weakbioturbation.RData')
load('SM6_strongbiomixing.RData')
load('SM6_strongbioirrigation.RData')
load('SM6_strongbioturbation.RData')

#======== MODEL DOMAIN AND GRID DEFINITION =========#
L <- 15   # depth of sediment domain [cm]
N <- 400  # number of grid layers
grid <- setup.grid.1D(x.up=0, x.down=L, N=N, dx.1=L/2000, p.dx.1=1.1)
Depth <- grid$x.mid 

#====== PLOTTING OUTPUTS ======#
#plot aerobic respiration, oxidation, and total oxygen consumption

nobio.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(nobio.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
nobio.cons$Depth <- Depth
nobio.cons$AR <- results.nobio$R1
nobio.cons$CSO <- results.nobio$R4
nobio.cons$Nox <- results.nobio$R5
nobio.cons$FeSox <- results.nobio$R6
nobio.cons$Tot_consumption <- results.nobio$R4 + results.nobio$R5 + results.nobio$R6 + results.nobio$R1
nobio.cons$Oxidations <- results.nobio$R4 + results.nobio$R5 + results.nobio$R6
nobio.cons$Bioturbation <- 'None'

weakbiomix.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(weakbiomix.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
weakbiomix.cons$Depth <- Depth
weakbiomix.cons$AR <- results.weakbiomix$R1
weakbiomix.cons$CSO <- results.weakbiomix$R4
weakbiomix.cons$Nox <- results.weakbiomix$R5
weakbiomix.cons$FeSox <- results.weakbiomix$R6
weakbiomix.cons$Tot_consumption <- results.weakbiomix$R4 + results.weakbiomix$R5 + results.weakbiomix$R6 + results.weakbiomix$R1
weakbiomix.cons$Oxidations <- results.weakbiomix$R4 + results.weakbiomix$R5 + results.weakbiomix$R6
weakbiomix.cons$Bioturbation <- 'Weak Biomixing'

weakbioirr.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(weakbioirr.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
weakbioirr.cons$Depth <- Depth
weakbioirr.cons$AR <- results.weakbioirr$R1
weakbioirr.cons$CSO <- results.weakbioirr$R4
weakbioirr.cons$Nox <- results.weakbioirr$R5
weakbioirr.cons$FeSox <- results.weakbioirr$R6
weakbioirr.cons$Tot_consumption <- results.weakbioirr$R4 + results.weakbioirr$R5 + results.weakbioirr$R6 + results.weakbioirr$R1
weakbioirr.cons$Oxidations <- results.weakbioirr$R4 + results.weakbioirr$R5 + results.weakbioirr$R6
weakbioirr.cons$Bioturbation <- 'Weak Bioirrigation'

weakbio.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(weakbio.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
weakbio.cons$Depth <- Depth
weakbio.cons$AR <- results.weak$R1
weakbio.cons$CSO <- results.weak$R4
weakbio.cons$Nox <- results.weak$R5
weakbio.cons$FeSox <- results.weak$R6
weakbio.cons$Tot_consumption <- results.weak$R4 + results.weak$R5 + results.weak$R6 + results.weak$R1
weakbio.cons$Oxidations <- results.weak$R4 + results.weak$R5 + results.weak$R6
weakbio.cons$Bioturbation <- 'Weak Bioturbation'

strongbiomix.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(strongbiomix.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
strongbiomix.cons$Depth <- Depth
strongbiomix.cons$AR <- results.strongbiomix$R1
strongbiomix.cons$CSO <- results.strongbiomix$R4
strongbiomix.cons$Nox <- results.strongbiomix$R5
strongbiomix.cons$FeSox <- results.strongbiomix$R6
strongbiomix.cons$Tot_consumption <- results.strongbiomix$R4 + results.strongbiomix$R5 + results.strongbiomix$R6 + results.strongbiomix$R1
strongbiomix.cons$Oxidations <- results.strongbiomix$R4 + results.strongbiomix$R5 + results.strongbiomix$R6
strongbiomix.cons$Bioturbation <- 'Strong Bioturbation'

strongbioirr.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(strongbioirr.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
strongbioirr.cons$Depth <- Depth
strongbioirr.cons$AR <- results.strongbioirr$R1
strongbioirr.cons$CSO <- results.strongbioirr$R4
strongbioirr.cons$Nox <- results.strongbioirr$R5
strongbioirr.cons$FeSox <- results.strongbioirr$R6
strongbioirr.cons$Tot_consumption <- results.strongbioirr$R4 + results.strongbioirr$R5 + results.strongbioirr$R6 + results.strongbioirr$R1
strongbioirr.cons$Oxidations <- results.strongbioirr$R4 + results.strongbioirr$R5 + results.strongbioirr$R6
strongbioirr.cons$Bioturbation <- 'Strong Bioturbation'

strongbio.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(strongbio.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
strongbio.cons$Depth <- Depth
strongbio.cons$AR <- results.strong$R1
strongbio.cons$CSO <- results.strong$R4
strongbio.cons$Nox <- results.strong$R5
strongbio.cons$FeSox <- results.strong$R6
strongbio.cons$Tot_consumption <- results.strong$R4 + results.strong$R5 + results.strong$R6 + results.strong$R1
strongbio.cons$Oxidations <- results.strong$R4 + results.strong$R5 + results.strong$R6
strongbio.cons$Bioturbation <- 'Strong Bioturbation'


#=====GGPLOT CONSTRUCTIO=====#

nobioturbation.consumption <- ggplot(data=nobio.cons, aes(x=AR, y=Depth)) + 
  ggtitle('No Bioturbation') +
  xlab('Consumption Rate (umol cm2 y-1)') +
  ylab('Depth (cm)') +
  geom_path(data=nobio.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=nobio.cons, aes(xmin=0, xmax=AR, y=Depth), fill='lightblue', color='darkblue', alpha=0.5) +
  geom_ribbon(data=nobio.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill='lightgreen', color='darkgreen', alpha=0.5) +
  #xlim(c(0, 1750)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() 

weakbiomixing.consumption <- ggplot(data=weakbiomix.cons, aes(x=AR, y=Depth)) + 
  ggtitle('Weak Biomixing') +
  xlab('Consumption Rate (umol cm2 y-1)') +
  ylab('Depth (cm)') +
  geom_path(data=weakbiomix.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=weakbiomix.cons, aes(xmin=0, xmax=AR, y=Depth), fill='lightblue', color='darkblue', alpha=0.5) +
  geom_ribbon(data=weakbiomix.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill='lightgreen', color='darkgreen', alpha=0.5) +
  #xlim(c(0, 1750)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )

weakbioirrigation.consumption <- ggplot(data=weakbioirr.cons, aes(x=AR, y=Depth)) + 
  ggtitle('Weak Bioirrigation') +
  xlab('Consumption Rate (umol cm2 y-1)') +
  ylab('Depth (cm)') +
  geom_path(data=weakbioirr.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=weakbioirr.cons, aes(xmin=0, xmax=AR, y=Depth), fill='lightblue', color='darkblue', alpha=0.5) +
  geom_ribbon(data=weakbioirr.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill='lightgreen', color='darkgreen', alpha=0.5) +
  #xlim(c(0, 1750)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )

weakbioturbation.consumption <- ggplot(data=weakbio.cons, aes(x=AR, y=Depth)) + 
  ggtitle('Weak Bioturbation') +
  xlab('Consumption Rate (umol cm2 y-1)') +
  ylab('Depth (cm)') +
  geom_path(data=weakbio.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=weakbio.cons, aes(xmin=0, xmax=AR, y=Depth), fill='lightblue', color='darkblue', alpha=0.5) +
  geom_ribbon(data=weakbio.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill='lightgreen', color='darkgreen', alpha=0.5) +
  #xlim(c(0, 1750)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )


strongbiomixing.consumption <- ggplot(data=strongbiomix.cons, aes(x=AR)) + 
  ggtitle('Strong Biomixing') +
  xlab('Consumption Rate (umol cm2 y-1)') +
  ylab('Depth (cm)') +
  geom_path(data=strongbiomix.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=strongbiomix.cons, aes(xmin=0, xmax=AR, y=Depth), fill='lightblue', color='darkblue', alpha=0.5) +
  geom_ribbon(data=strongbiomix.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill='lightgreen', color='darkgreen', alpha=0.5) +
  #xlim(c(0, 1750)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )

strongbioirrigation.consumption <- ggplot(data=strongbioirr.cons, aes(x=AR)) + 
  ggtitle('Strong Bioirrigation') +
  xlab('Consumption Rate (umol cm2 y-1)') +
  ylab('Depth (cm)') +
  geom_path(data=strongbioirr.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=strongbioirr.cons, aes(xmin=0, xmax=AR, y=Depth), fill='lightblue', color='darkblue', alpha=0.5) +
  geom_ribbon(data=strongbioirr.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill='lightgreen', color='darkgreen', alpha=0.5) +
  #xlim(c(0, 1750)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )


strongbioturbation.consumption <- ggplot(data=strongbio.cons, aes(x=AR)) + 
  ggtitle('Strong Bioturbation') +
  xlab('Consumption Rate (umol cm2 y-1)') +
  ylab('Depth (cm)') +
  geom_path(data=strongbio.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=strongbio.cons, aes(xmin=0, xmax=AR, y=Depth), fill='lightblue', color='darkblue', alpha=0.5) +
  geom_ribbon(data=strongbio.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill='lightgreen', color='darkgreen', alpha=0.5) +
  #xlim(c(0, 1750)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )

blank <- ggplot() + geom_point() + theme_classic() +
  theme(
    axis.line.x=element_blank(),
    axis.line.y=element_blank()
  )


ggarrange(nobioturbation.consumption, weakbiomixing.consumption, weakbioirrigation.consumption, weakbioturbation.consumption,
          blank, strongbiomixing.consumption, strongbioirrigation.consumption, strongbioturbation.consumption, ncol=4)

