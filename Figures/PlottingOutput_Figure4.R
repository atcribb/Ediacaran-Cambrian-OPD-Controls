#Author: Alison Cribb
#Plotting output to re-create Figure 4

#==== REQUIRED PACKAGES ====#
library(ggplot2)
library(egg)
library(tayloRswift)

#==== DATA INPUT ====#
#load low endmember (0.014 mM) results
load('Figure4_lowend.RData')
results.low$Bioturbation <- factor(results.low$Bioturbation, levels=c('None', 'Ediacaran', 'Terreneuvian', 'Modern'))
#load high endmember (0.15 mM) results
load('Figure4_highend.RData')
results.high$Bioturbation <- factor(results.low$Bioturbation, levels=c('None', 'Ediacaran', 'Terreneuvian', 'Modern'))

#===== PLOT SET UP ====#
#set color scales: low end is blue, high end is red
#high end (FOAM):
nobio_col     <- '#BFBCAA'
ediacaran_col   <- '#A6836F'
terre_col <- '#73564C'
modern_col <- '#400303'
high_biolevel_cols <- c(nobio_col, ediacaran_col, terre_col, modern_col)
plot(1:4, 1:4, pch=19, cex=4, col=high_biolevel_cols)

#low end:
nobio_col     <- '#AAB7BF'
ediacaran_col   <- '#6F93A6'
terre_col <- '#4C6573'
modern_col <- '#032440'
low_biolevel_cols <- c(nobio_col, ediacaran_col, terre_col, modern_col)
plot(1:4, 1:4, pch=19, cex=4, col=low_biolevel_cols)


#===== GGPLOT =====#
highend_plot <- ggplot(data=results.high, aes(x=F.CH2O, y=OPD, fill=Bioturbation)) +
  geom_point(size=6, pch=21, colour='black', alpha=0.8, aes(shape=Bioturbation)) +
  scale_fill_manual(values=high_biolevel_cols) +
  scale_x_continuous(expression(paste(CH[2],'O'~Flux~'(',mu,mol~cm^{'-2'}~yr^{'-1'},')')), breaks=c(100, 250, 450)) +
  #xlab('Organic Matter Flux (cm^2 yr^-1)') + 
  ylab('Oxygen Penetration Depth (cm)') +
  scale_y_reverse(limits=c(10, 0)) +
  ggtitle(expression(paste(O['2,BW']~'='~'0.15'~'mM'))) +
  theme_bw() +
  theme(
    axis.text.x=element_text(size=12, color='black'),
    axis.text.y=element_blank(),
    axis.title.x=element_text(size=15, color='black'),
    axis.title.y=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13),
    legend.position=c(0.82,0.15)
  )

highend_plot


lowend_plot <- ggplot(data=results.low, aes(x=F.CH2O, y=OPD, fill=Bioturbation)) +
  geom_point(size=6, pch=21, colour='black', alpha=0.8) +
  scale_fill_manual(values=low_biolevel_cols) +
  scale_x_continuous(expression(paste(CH[2],'O'~Flux~'(',mu,mol~cm^{'-2'}~yr^{'-1'},')')), breaks=c(100, 250, 450)) +
  #xlab('Organic Matter Flux (cm^2 yr^-1)') + 
  ylab('Oxygen Penetration Depth (cm)') +
  #ylim(c(0, 10)) +
  scale_y_reverse(limits=c(10, 0)) +
  ggtitle(expression(paste(O['2,BW']~'='~'0.014'~'mM'))) +
  theme_bw() +
  theme(
    axis.text.x=element_text(size=12, color='black'),
    axis.text.y=element_text(size=12, color='black'),
    axis.title.x=element_text(size=15, color='black'),
    axis.title.y=element_text(size=15, color='black'),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13),
    legend.position=c(0.82,0.15)
  )

lowend_plot

combined <- ggarrange(lowend_plot, highend_plot, ncol=2)
combined

ggsave('O2EndMembers_raw.PDF', combined, width=12, height=6.5)
