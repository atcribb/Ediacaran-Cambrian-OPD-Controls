#Author: Alison Cribb
#Plotting output to recreate oxygen consumption profiles in Figure 5 and Figure S3

#==== REQUIRED PACKAGES ====#
library(ggplot2)
library(egg)
library(ReacTran)

#====LOAD DATA====#
load('O2cons_nobioturbation.RData')
load('O2cons_ediacaranbiomixing.RData')
load('O2cons_ediacaranbioirrigation.RData')
load('O2cons_ediacaranbioturbation.RData')
load('O2cons_terreneuvianbiomixing.RData')
load('O2cons_terreneuvianbioirrigation.RData')
load('O2cons_terreneuvianbioturbation.RData')
load('O2cons_modernbiomixing.RData')
load('O2cons_modernbioirrigation.RData')
load('O2cons_modernbioturbation.RData')

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

edibiomix.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(edibiomix.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
edibiomix.cons$Depth <- Depth
edibiomix.cons$AR <- results.edibiomix$R1
edibiomix.cons$CSO <- results.edibiomix$R4
edibiomix.cons$Nox <- results.edibiomix$R5
edibiomix.cons$FeSox <- results.edibiomix$R6
edibiomix.cons$Tot_consumption <- results.edibiomix$R4 + results.edibiomix$R5 + results.edibiomix$R6 + results.edibiomix$R1
edibiomix.cons$Oxidations <- results.edibiomix$R4 + results.edibiomix$R5 + results.edibiomix$R6
edibiomix.cons$Bioturbation <- 'Ediacaran Biomixing'

edibioirr.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(edibioirr.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
edibioirr.cons$Depth <- Depth
edibioirr.cons$AR <- results.edibioirr$R1
edibioirr.cons$CSO <- results.edibioirr$R4
edibioirr.cons$Nox <- results.edibioirr$R5
edibioirr.cons$FeSox <- results.edibioirr$R6
edibioirr.cons$Tot_consumption <- results.edibioirr$R4 + results.edibioirr$R5 + results.edibioirr$R6 + results.edibioirr$R1
edibioirr.cons$Oxidations <- results.edibioirr$R4 + results.edibioirr$R5 + results.edibioirr$R6
edibioirr.cons$Bioturbation <- 'Ediacaran Bioirrigation'

ediacaran.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(ediacaran.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
ediacaran.cons$Depth <- Depth
ediacaran.cons$AR <- results.ediacaran$R1
ediacaran.cons$CSO <- results.ediacaran$R4
ediacaran.cons$Nox <- results.ediacaran$R5
ediacaran.cons$FeSox <- results.ediacaran$R6
ediacaran.cons$Tot_consumption <- results.ediacaran$R4 + results.ediacaran$R5 + results.ediacaran$R6 + results.ediacaran$R1
ediacaran.cons$Oxidations <- results.ediacaran$R4 + results.ediacaran$R5 + results.ediacaran$R6
ediacaran.cons$Bioturbation <- 'Ediacaran Bioturbation'

terrebiomix.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(terrebiomix.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
terrebiomix.cons$Depth <- Depth
terrebiomix.cons$AR <- results.terrebiomix$R1
terrebiomix.cons$CSO <- results.terrebiomix$R4
terrebiomix.cons$Nox <- results.terrebiomix$R5
terrebiomix.cons$FeSox <- results.terrebiomix$R6
terrebiomix.cons$Tot_consumption <- results.terrebiomix$R4 + results.terrebiomix$R5 + results.terrebiomix$R6 + results.terrebiomix$R1
terrebiomix.cons$Oxidations <- results.terrebiomix$R4 + results.terrebiomix$R5 + results.terrebiomix$R6
terrebiomix.cons$Bioturbation <- 'Terreneuvian Biomixing'

terrebioirr.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(terrebioirr.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
terrebioirr.cons$Depth <- Depth
terrebioirr.cons$AR <- results.terrebioirr$R1
terrebioirr.cons$CSO <- results.terrebioirr$R4
terrebioirr.cons$Nox <- results.terrebioirr$R5
terrebioirr.cons$FeSox <- results.terrebioirr$R6
terrebioirr.cons$Tot_consumption <- results.terrebioirr$R4 + results.terrebioirr$R5 + results.terrebioirr$R6 + results.terrebioirr$R1
terrebioirr.cons$Oxidations <- results.terrebioirr$R4 + results.terrebioirr$R5 + results.terrebioirr$R6
terrebioirr.cons$Bioturbation <- 'Terreneuvian Bioirrigation'

terreneuvian.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(terreneuvian.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
terreneuvian.cons$Depth <- Depth
terreneuvian.cons$AR <- results.terreneuvian$R1
terreneuvian.cons$CSO <- results.terreneuvian$R4
terreneuvian.cons$Nox <- results.terreneuvian$R5
terreneuvian.cons$FeSox <- results.terreneuvian$R6
terreneuvian.cons$Tot_consumption <- results.terreneuvian$R4 + results.terreneuvian$R5 + results.terreneuvian$R6 + results.terreneuvian$R1
terreneuvian.cons$Oxidations <- results.terreneuvian$R4 + results.terreneuvian$R5 + results.terreneuvian$R6
terreneuvian.cons$Bioturbation <- 'Terreneuvian Bioturbation'

modernbiomix.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(modernbiomix.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
modernbiomix.cons$Depth <- Depth
modernbiomix.cons$AR <- results.modernbiomix$R1
modernbiomix.cons$CSO <- results.modernbiomix$R4
modernbiomix.cons$Nox <- results.modernbiomix$R5
modernbiomix.cons$FeSox <- results.modernbiomix$R6
modernbiomix.cons$Tot_consumption <- results.modernbiomix$R4 + results.modernbiomix$R5 + results.modernbiomix$R6 + results.modernbiomix$R1
modernbiomix.cons$Oxidations <- results.modernbiomix$R4 + results.modernbiomix$R5 + results.modernbiomix$R6
modernbiomix.cons$Bioturbation <- 'Modern Biomixing'

modernbioirr.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(modernbioirr.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
modernbioirr.cons$Depth <- Depth
modernbioirr.cons$AR <- results.modernbioirr$R1
modernbioirr.cons$CSO <- results.modernbioirr$R4
modernbioirr.cons$Nox <- results.modernbioirr$R5
modernbioirr.cons$FeSox <- results.modernbioirr$R6
modernbioirr.cons$Tot_consumption <- results.modernbioirr$R4 + results.modernbioirr$R5 + results.modernbioirr$R6 + results.modernbioirr$R1
modernbioirr.cons$Oxidations <- results.modernbioirr$R4 + results.modernbioirr$R5 + results.modernbioirr$R6
modernbioirr.cons$Bioturbation <- 'Modern Bioirrigation'

modern.cons <- as.data.frame(matrix(data=NA, nrow=length(Depth), ncol=8))
colnames(modern.cons) <- c('Depth', 'AR', 'CSO', 'Nox', 'FeSox', 'Oxidations', 'Tot_consumption', 'Bioturbation')
modern.cons$Depth <- Depth
modern.cons$AR <- results.modern$R1
modern.cons$CSO <- results.modern$R4
modern.cons$Nox <- results.modern$R5
modern.cons$FeSox <- results.modern$R6
modern.cons$Tot_consumption <- results.modern$R4 + results.modern$R5 + results.modern$R6 + results.modern$R1
modern.cons$Oxidations <- results.modern$R4 + results.modern$R5 + results.modern$R6
modern.cons$Bioturbation <- 'Modern Bioturbation'


#=====GGPLOT CONSTRUCTIO=====#
fillcols <- c('#A6836F', '#445188')
colorcols<- c('#675050', '#0F1F7A')

nobioturbation.consumption <- ggplot(data=nobio.cons, aes(x=AR, y=Depth)) + 
  ggtitle('No Bioturbation') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{2}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=nobio.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=nobio.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=nobio.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() 

ediacaranbiomixing.consumption <- ggplot(data=edibiomix.cons, aes(x=AR, y=Depth)) + 
  ggtitle('Ediacaran Biomixing') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{2}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=edibiomix.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=edibiomix.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=edibiomix.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )

ediacaranbioirrigation.consumption <- ggplot(data=edibioirr.cons, aes(x=AR, y=Depth)) + 
  ggtitle('Ediacaran Bioirrigation') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{2}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=edibioirr.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=edibioirr.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=edibioirr.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )

ediacaranbioturbation.consumption <- ggplot(data=ediacaran.cons, aes(x=AR, y=Depth)) + 
  ggtitle('Ediacaran Bioturbation') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{3}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=ediacaran.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=ediacaran.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=ediacaran.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )

terreneuvianbiomixing.consumption <- ggplot(data=terrebiomix.cons, aes(x=AR, y=Depth)) + 
  ggtitle('Terreneuvian Biomixing') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{2}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=terrebiomix.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=terrebiomix.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=terrebiomix.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    #axis.title.y = element_blank()
  )

terreneuvianbioirrigation.consumption <- ggplot(data=terrebioirr.cons, aes(x=AR, y=Depth)) + 
  ggtitle('Terreneuvian Bioirrigation') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{2}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=terrebioirr.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=terrebioirr.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=terrebioirr.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )

terreneuvianbioturbation.consumption <- ggplot(data=terreneuvian.cons, aes(x=AR, y=Depth)) + 
  ggtitle('Terreneuvian Bioturbation') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{3}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=terreneuvian.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=terreneuvian.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=terreneuvian.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )

modernbiomixing.consumption <- ggplot(data=modernbiomix.cons, aes(x=AR)) + 
  ggtitle('Modern Biomixing') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{3}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=modernbiomix.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=modernbiomix.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=modernbiomix.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    #axis.title.y = element_blank()
  )

modernbioirrigation.consumption <- ggplot(data=modernbioirr.cons, aes(x=AR)) + 
  ggtitle('Modern Bioirrigation') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{3}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=modernbioirr.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=modernbioirr.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=modernbioirr.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
  scale_y_reverse(limits=c(1, -0.01)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )


modernbioturbation.consumption <- ggplot(data=modern.cons, aes(x=AR)) + 
  ggtitle('Modern Bioturbation') +
  xlab(expression(paste(Consumption~Rate~(mu~mol~cm^{3}~yr^{-1})))) +
  ylab('Depth (cm)') +
  geom_path(data=modern.cons, aes(x=Tot_consumption, y=Depth), color='black', linetype='dotted') +
  geom_ribbon(data=modern.cons, aes(xmin=0, xmax=AR, y=Depth), fill=fillcols[1], color=colorcols[1], alpha=0.5) +
  geom_ribbon(data=modern.cons, aes(xmin=0, xmax=Oxidations, y=Depth), fill=fillcols[2], color=colorcols[2], alpha=0.5) +
  xlim(c(0, 8000)) +
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





Figure5_bioturbationconsprofiles <- ggarrange(nobioturbation.consumption, ediacaranbioturbation.consumption, terreneuvianbioturbation.consumption, modernbioturbation.consumption, ncol=4)
ggsave('O2consumptionprofiles_main.pdf', Figure5_bioturbationconsprofiles, width=11, height=5)


FigureS3_allconsprofiles <- ggarrange(nobioturbation.consumption, ediacaranbiomixing.consumption, ediacaranbioirrigation.consumption, ediacaranbioturbation.consumption,
                             blank,                      terreneuvianbiomixing.consumption, terreneuvianbioirrigation.consumption, terreneuvianbioturbation.consumption,
                             blank,                      modernbiomixing.consumption, modernbioirrigation.consumption, modernbioturbation.consumption,
                             ncol=4)
ggsave('O2consumptionprofiles_all.pdf', FigureS3_allconsprofiles, width=11, height=15)
