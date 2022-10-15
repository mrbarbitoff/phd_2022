setwd("/media/barbitoff/DATA/Working issues/FizGen/Hsp40 mechanisms in vitro/Refolding/2022-04-19-Refolding-rep1")
library(ggplot2)
library(cowplot)

lucdata <- read.table('./refolding_rep1.tsv', sep='\t', header=T)
lucdata$Chaperones = factor(lucdata$Chaperones, 
        levels=c('Sis1+Ssa1', 'Ssa1', 'none'))

a <- ggplot(lucdata, aes(x=Chaperones, y=MeanActivity, 
                    ymax=MeanActivity+SD, ymin=MeanActivity-SD)) +
  geom_bar(stat='identity', col='black', fill='gray') + 
  geom_errorbar(width=0.4, lwd=0.5) + theme_bw() + coord_flip()
a

plot_grid(a, a, nrow=2)
