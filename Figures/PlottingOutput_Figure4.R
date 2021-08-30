#Plotting output to re-create Figure 4

library(ggplot2)
library(egg)

#==== DATA INPUT ====#
#load low endmember (0.014 mM) results
load('Figure4_lowend.RData')
results.low$Bioturbation <- factor(results.low$Bioturbation, levels=c('None', 'Weak', 'Strong'))
#load high endmember (0.15 mM) results
load('Figure4_highend.RData')
results.high$Bioturbation <- factor(results.high$Bioturbation, levels=c('None', 'Weak', 'Strong'))

#===== PLOT SET UP ====#
#set color scales: low end is blue, high end is red
#high end (FOAM):
nobio_col <- '#F9DEDC'
weakbio_col <- '#F7B538'
strongbio_col <- '#B50223'
high_biolevel_cols <- c(nobio_col, weakbio_col, strongbio_col)

#low end:
nobio_col     <- '#D8D8D8' 
weakbio_col   <- '#AEC3B0' 
strongbio_col <- '#124559' 
low_biolevel_cols <- c(nobio_col, weakbio_col, strongbio_col)


#===== GGPLOT =====#
highend_plot <- ggplot(data=results.high, aes(x=F.CH2O, y=OPD, fill=Bioturbation)) +
  geom_point(size=6, pch=21, colour='black', alpha=0.8) +
  scale_fill_manual(values=high_biolevel_cols) +
  scale_x_continuous(expression(paste(CH[2],'O'~Flux~'(',mu,mol~cm^{'-2'}~yr^{'-1'},')')), breaks=c(100, 250, 450)) +
  #xlab('Organic Matter Flux (cm^2 yr^-1)') + 
  ylab('Oxygen Penetration Depth (cm)') +
  #ylim(c(0, 10)) +
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

ggarrange(lowend_plot, highend_plot, ncol=2)
