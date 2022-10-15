library(ggplot2)

setwd("/media/barbitoff/DATA/Working issues/FizGen/Hsp40 mechanisms in vitro/PhD/Disser_v0/Data analysis/Ssa1-NM-validation/")

dens = read.table('dens.csv', sep=',', header=T)
head(dens)

ggplot(dens, aes(x=j, y=pct_pellet, fill=j)) + 
  geom_boxplot(outlier.shape=NA) +
  scale_fill_brewer(palette = "Accent") + theme_bw() +
  facet_wrap(~condition, nrow=1) + guides(fill=F) +
  geom_jitter(pch=21)

p1 <- wilcox.test(pct_pellet~j, dens[dens$j %in% c('Sis1', 'none') & 
                                       dens$condition == 'NM+AMP-PNP',  ])$p.value
p2 <- wilcox.test(pct_pellet~j, dens[dens$j %in% c('Sis1', 'none') & 
                                       dens$condition == 'NM+ATP',  ])$p.value
p3 <- wilcox.test(pct_pellet~j, dens[dens$j %in% c('Sis1', 'Sis1dDD') & 
                                 dens$condition == 'NM+ATP',  ])$p.value
#wilcox.test(pct_pellet~j, dens[(dens$j == 'Sis1dDD' & dens$condition == 'NM+AMP-PNP') |
#                    (dens$j == 'none' & dens$condition == 'NM+ATP'),  ])
p4 <- wilcox.test(pct_pellet~j, dens[dens$j %in% c('Sis1', 'none') & 
                                       dens$condition == 'NM+Ssa1-21',  ])$p.value
p5 <- wilcox.test(pct_pellet~j, dens[dens$j %in% c('Sis1', 'Sis1dDD') & 
                                       dens$condition == 'NM+Ssa1-21',  ])$p.value
p6 <- wilcox.test(pct_pellet~condition, dens[dens$j == 'none' & 
                                               dens$condition %in% c('NM+ATP', 'NM+Ssa1-21'),  ])$p.value
p7 <- wilcox.test(pct_pellet~condition, dens[dens$j == 'Sis1' & 
                                               dens$condition %in% c('NM+ATP', 'NM+Ssa1-21'),  ])$p.value


pvals <- p.adjust(c(p1, p2, p3, p4, p5, p6, p7), method = "fdr")
pvals

aggregate(pct_pellet~j+condition, dens, median)


