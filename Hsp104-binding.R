setwd("/media/barbitoff/DATA/Working issues/FizGen/Hsp40 mechanisms in vitro/PhD/Disser_v0/Data analysis/Hsp104-binding/")
library(ggplot2)
mycol1 = rgb(85, 138, 221, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 121, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

# Change to 'data_nm.csv'
data_q1 = read.table('data_hsp104.csv', sep=',', header=T)
head(data_q1)

data_q1$bound = data_q1$Pellet / (data_q1$Super + data_q1$Pellet)

ggplot(data_q1, aes(x=Ligand, y=bound, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]')

baseline1 = apply(data_q1, 1, function(x) mean(data_q1[data_q1$Chaperone == x[5] & data_q1$Ligand == 0, 'bound']))
baseline1

data_q1$corr = ifelse(data_q1$Ligand == 0, 0, (data_q1$bound - baseline1)/(1-baseline1))
data_q1$corr[data_q1$corr < 0] = 0

myformula = formula(bound~(1-a)*(Ligand/(K+Ligand))+a)

bestfit_wt_atp <- nls(myformula, data_q1[data_q1$Chaperone == 'wt+ATP', ], start=list(K=2, a = 0))
summary(bestfit_wt_atp)
bestfit_dn_atp <- nls(myformula, data_q1[data_q1$Chaperone == 'deltaN+ATP', ], start=list(K=2, a=0))
summary(bestfit_dn_atp)
bestfit_wt_amp <- nls(myformula, data_q1[data_q1$Chaperone == 'wt+AMP-PNP', ], start=list(K=2, a=0))
summary(bestfit_wt_amp)
bestfit_wt_atp_sis <- nls(myformula, data_q1[data_q1$Chaperone == 'wt+Sis1+ATP', ], start=list(K=2, a=0))
summary(bestfit_wt_atp_sis)

coeffs1 = data.frame(K = c(summary(bestfit_wt_atp)$coeff[1], 
                           summary(bestfit_wt_amp)$coeff[1],
                           summary(bestfit_dn_atp)$coeff[1],
                           summary(bestfit_wt_atp_sis)$coeff[1]), 
                     se = c(summary(bestfit_wt_atp)$coeff[3], 
                            summary(bestfit_wt_amp)$coeff[3],
                            summary(bestfit_dn_atp)$coeff[3],
                            summary(bestfit_wt_atp_sis)$coeff[3]),
                     Chaperone = c('wt+ATP', 'wt+AMP-PNP', 'deltaN+ATP', 'wt+Sis1+ATP'))
coeffs1


reg.t.test <- function(k.a, k.b, se.a, se.b, df) {
  z = (k.a - k.b)/sqrt(se.a^2 + se.b^2)
  return(pt(abs(z), lower.tail=F, df = df))
}

data_q1$j = c(rep('none', 32), rep('Sis1', 12))

ggplot(data_q1, aes(x=Ligand, y=corr, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', 
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2)), se=FALSE) + 
  scale_color_manual(values=c(mycol1, mycol4, mycol3, mycol6))

kd_104 <- ggplot(coeffs1, aes(x=Chaperone, y=K, fill=Chaperone)) + geom_bar(col='black', stat='identity') +
  geom_errorbar(ymin=coeffs1$K-coeffs1$se, ymax=coeffs1$K+coeffs1$se, width=0.5) +
  scale_y_continuous(limits=c(0, 1000)) +
  scale_fill_manual(values=c(mycol2, mycol1, mycol3, mycol5)) + theme_bw() +
  ylab('Estimated Kd (uM)') + coord_flip()
kd_104
