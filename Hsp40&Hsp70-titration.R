setwd("/media/barbitoff/DATA/Working issues/FizGen/Hsp40 mechanisms in vitro/PhD/Disser_v0/Data analysis/Hsp40&Hsp70-titration/")
library(ggplot2)
mycol1 = rgb(190, 174, 212, maxColorValue = 255)
mycol2 = rgb(253, 192, 134, maxColorValue = 255)
mycol3 = rgb(255, 121, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

# Change to 'data_nm.csv'
data_q1 = read.table('data_nm_iter2.csv', sep=',', header=T)
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

bestfit_sis1_nm <- nls(myformula, data_q1[data_q1$Chaperone == 'Sis1', ], start=list(K=2, a=0))
summary(bestfit_sis1_nm)
bestfit_dd_nm <- nls(myformula, data_q1[data_q1$Chaperone == 'deltaDD', ], start=list(K=2, a=0))
summary(bestfit_dd_nm)

reg.t.test <- function(k.a, k.b, se.a, se.b, df) {
  z = (k.a - k.b)/sqrt(se.a^2 + se.b^2)
  return(pt(abs(z), lower.tail=F, df = df))
}

reg.t.test(16.81, 76.13, 2.28, 23.90, df=20)

#### Another replicate

data_q1_r1 = read.table('data_nm.csv', sep=',', header=T)
head(data_q1_r1)

data_q1_r1$bound = data_q1_r1$Pellet / (data_q1_r1$Super + data_q1_r1$Pellet)

bestfit_sis1_nm_r1 <- nls(myformula, data_q1_r1[data_q1_r1$Chaperone == 'Sis1', ], start=list(K=2, a=0))
summary(bestfit_sis1_nm_r1)
bestfit_dd_nm_r1 <- nls(myformula, data_q1_r1[data_q1_r1$Chaperone == 'deltaDD', ], start=list(K=2, a=0))
summary(bestfit_dd_nm_r1)
bestfit_ydj_nm_r1 <- nls(myformula, data_q1_r1[data_q1_r1$Chaperone == 'Ydj1', ], start=list(K=2, a=0))
summary(bestfit_ydj_nm_r1)


reg.t.test(2.53, 22.77, 0.73, 5.34, df=32)



coeffs1 = data.frame(K = c(summary(bestfit_sis1_nm)$coeff[1], summary(bestfit_dd_nm)$coeff[1]), 
                    se = c(summary(bestfit_sis1_nm)$coeff[3], summary(bestfit_dd_nm)$coeff[3]),
                    Chaperone = c('Sis1', 'deltaDD'))
coeffs1$Fibril = 'NM'
coeffs1

ggplot(data_q1, aes(x=Ligand, y=corr, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', 
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2)), se=FALSE) + 
  scale_color_manual(values=c(mycol2, mycol1, mycol3))

ggplot(coeffs1, aes(x=Chaperone, y=K, fill=Chaperone)) + geom_bar(col='black', stat='identity') +
  geom_errorbar(ymin=coeffs1$K-coeffs1$se, ymax=coeffs1$K+coeffs1$se, width=0.5) +
  scale_y_continuous(limits=c(0, 120)) +
  scale_fill_manual(values=c(mycol2, mycol1, mycol3)) + theme_bw() +
  ylab('Estimated Kd (uM)')


# Rnq1 goes here

# Change tp 'data_rnq.csv'
data_q2 = read.table('data_rnq1_iter2.csv', sep=',', header=T)
head(data_q2)

data_q2$bound = data_q2$Pellet / (data_q2$Super + data_q2$Pellet)

ggplot(data_q2, aes(x=Ligand, y=bound, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]')

baseline2 = apply(data_q2, 1, function(x) mean(data_q2[data_q2$Chaperone == x[5] & data_q2$Ligand == 0, 'bound']))
baseline2

data_q2$corr = ifelse(data_q2$Ligand == 0, 0, (data_q2$bound - baseline2)/(1-baseline2))
data_q2$corr[data_q2$corr < 0] = 0

data_q2[data_q2$corr == 0 & data_q2$Ligand > 0, 'corr'] = 0.001


bestfit_sis1_rnq <- nls(myformula, data_q2[data_q2$Chaperone == 'Sis1', ], start=list(K=2, a=0))
summary(bestfit_sis1_rnq)
bestfit_dd_rnq <- nls(myformula, data_q2[data_q2$Chaperone == 'deltaDD', ], start=list(K=2, a=0))
summary(bestfit_dd_rnq)

ggplot(data_q2, aes(x=Ligand, y=corr, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', 
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2)), se=FALSE) + 
  scale_color_manual(values=c(mycol2, mycol1, mycol3))


# Hsp40 binding vs amyloidpgenic proteins
data_all= rbind(data_q1[, c(1:5, 7, 8)], data_q2)
data_all[data_all$corr == 0 & data_all$Ligand > 0, 'corr'] = 0.001

ggplot(data_all, aes(x=Ligand, y=corr, col=Chaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', 
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2),
                                                  control=nls.control(maxiter=100)),
                                 se=FALSE) + 
  scale_color_manual(values=c(mycol2, mycol1, mycol6)) + facet_wrap(~Fibril, nrow=2)


# SSA1 FOR NM

data_q3 = read.table('Ssa1_data.csv', sep=',', header=T)
head(data_q3)

data_q3$bound = data_q3$Pellet / (data_q3$Super + data_q3$Pellet)

ggplot(data_q3, aes(x=Ligand, y=bound, col=Cochaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]')

baseline3 = apply(data_q3, 1, function(x) mean(data_q3[data_q3$Cochaperone == x[5] & data_q3$Ligand == 0, 'bound']))
baseline3

data_q3$corr = ifelse(data_q3$Ligand == 0, 0, (data_q3$bound - baseline3)/(1-baseline3))
data_q3$corr[data_q3$corr < 0] = 0

bestfit_sis1_nm_ssa <- nls(myformula, data_q3[data_q3$Cochaperone == 'Sis1', ], start=list(K=2, a=0))
summary(bestfit_sis1_nm_ssa)
bestfit_dd_nm_ssa <- nls(myformula, data_q3[data_q3$Cochaperone == 'deltaDD', ], start=list(K=2, a=0))
summary(bestfit_dd_nm_ssa)
bestfit_none_nm_ssa <- nls(myformula, data_q3[data_q3$Cochaperone == 'none', ], start=list(K=2, a=0))
summary(bestfit_none_nm_ssa)

ggplot(data_q3, aes(x=Ligand, y=corr, col=Cochaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', 
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2)), se=FALSE) + 
  scale_color_manual(values=c(mycol2, mycol4, mycol1))


# Ssa1 with Rnq1

data_q4 = read.table('Ssa1_data_Rnq1.csv', sep=',', header=T)
head(data_q4)

data_q4$bound = data_q4$Pellet / (data_q4$Super + data_q4$Pellet)

ggplot(data_q4, aes(x=Ligand, y=bound, col=Cochaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]')

baseline4 = apply(data_q4, 1, 
                  function(x) mean(data_q4[data_q4$Cochaperone == x[5] & data_q4$Ligand == 0, 'bound']))
baseline4

data_q4$corr = ifelse(data_q4$Ligand == 0, 0, (data_q4$bound - baseline4)/(1-baseline4))
data_q4$corr[data_q4$corr < 0] = 0

bestfit_sis1_rnq_ssa <- nls(myformula, data_q4[data_q4$Cochaperone == 'Sis1', ], start=list(K=2, a=0))
summary(bestfit_sis1_rnq_ssa)
bestfit_dd_rnq_ssa <- nls(myformula, data_q4[data_q4$Cochaperone == 'deltaDD', ], start=list(K=2, a=0))
summary(bestfit_dd_rnq_ssa)
bestfit_none_rnq_ssa <- nls(myformula, data_q4[data_q4$Cochaperone == 'none', ], start=list(K=2, a=0))
summary(bestfit_none_rnq_ssa)

ggplot(data_q4, aes(x=Ligand, y=corr, col=Cochaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', fullrange = TRUE,
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2)), se=FALSE) + 
  scale_color_manual(values=c(mycol2, mycol4, mycol1))


# Merging Ssa1
data_ssa_all = rbind(data_q3, data_q4)

ggplot(data_ssa_all, aes(x=Ligand, y=corr, col=Cochaperone)) + geom_point(size=2.3) + 
  scale_y_continuous(limits=c(0, 1)) + theme_bw() + xlab('Fibril concentration (uM)') +
  ylab('[CF]/[C]') + geom_smooth(method='nls', fullrange = TRUE,
                                 formula = y~1*x/(K+x),
                                 method.args=list(start=c(K=2)), se=FALSE) + 
  scale_color_manual(values=c(mycol2, mycol4, mycol1)) + facet_wrap(~Fibril, nrow=2)

