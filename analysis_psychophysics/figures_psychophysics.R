
# paths
usrPath <- '/Users/bear/Downloads' # change to your local
wdPath <- file.path(usrPath, 'dynamic_contrast', 'analysis_psychophysics')
fitPath <- file.path(usrPath, 'dynamic_contrast', 'fitdata_psychophysics',
                     'psychophysics_fits.csv')
funcPath <- file.path(usrPath, 'dynamic_contrast', 'r_helper_functions.R')
# sourcing r_helper_functions has all the functions needed for analysis,
# also sets plot colors/pchs

setwd(wdPath)
source(funcPath)


load('descriptives_psychophysics.RData')
# loads data frame called, imaginatively, data

fitdata <- read.csv(fitPath)
data <- merge(data,fitdata, by = 'sID')
data$group <- factor(data$group)

# map between left/right and amb/fellow
data[,c('acuityAE', 'acuityFE')] <- map2amb(data$acuityLE, data$acuityRE, data$AE)
data[,c('AULCSF_AE', 'AULCSF_FE')] <- map2amb(data$l.AULCSF, data$r.AULCSF, data$AE)
data[,c('sens2cpd_AE', 'sens2cpd_FE')] <- map2amb(data$l.sens2cpd, data$r.sens2cpd, data$AE)
data[,c('LS', 'LS_fellow')] <- map2amb(data$LS_LE, data$LS_RE, data$AE)

data$IODacuity <- abs(data$acuityAE - data$acuityFE)

# bino summation/inhibition:  best minus bino (negative = inhibition, positive = summation)
data$bino_acuity_summation = apply(
  data[,c('acuityLE', 'acuityRE')], 1, FUN=min) - data$acuityBE

data$AULCSF_summation = data$b.AULCSF / apply(data[,c('AULCSF_AE', 'AULCSF_FE')], 1, FUN=max)


# clean descriptives not used in analysis - remove from this list to keep
data <- subset(data, select = -c(
  animals, shapes, diff_proportion,
  LE_1cpd, RE_1cpd,
  LE_2cpd, RE_2cpd,
  LE_4cpd, RE_4cpd,
  LE_8cpd, RE_8cpd,
  LS_fellow,
  l.AULCSF, r.AULCSF,
  b.peakgain, l.peakgain, r.peakgain,
  b.peakfreq, l.peakfreq, r.peakfreq,
  b.bandwidth, l.bandwidth, r.bandwidth,
  b.truncation, l.truncation, r.truncation,
  l.sens2cpd, r.sens2cpd,
  b.highSFcutoff, l.highSFcutoff, r.highSFcutoff,
  LS_RE, LS_LE))

# clean workspace
rm(fitdata)

# reduced data
redPath <- file.path(usrPath, 'dynamic_contrast', 'fitdata_psychophysics',
                     'psychophysics_reduced_fits.csv')
reddata <- read.csv(redPath)

if (!all(data$sID == reddata$sID)) {
  print('WARNING, SIDs DO NOT MATCH')
}

# Plots 

## Model parameters 

# kAE:
stripchart_w_means(round(data$kAE*3,1)/3, data$group, 'kAE', c(0,1),3)
dev.copy(pdf,width=4, height = 4.5, 'plot_kae.pdf')
dev.off()

# normalization parameters:
# remove the outlying value of 34
data$uAEcln <- data$uAE
data$uAEcln[data$uAEcln>33] <- NA

stripchart_w_means_byeye(round(data$uAEcln,2),
                         round(data$uFE,2),
                          data$group, data$group, 'u', c(0,10))
dev.copy(pdf,width=4.5, height = 4.5, 'plot_u.pdf')
dev.off()

# sigma:

stripchart_w_means(round(data$sigma*2,1)/2, data$group, 'sigma', c(0,2),2)
dev.copy(pdf,width=4, height = 4.5, 'plot_sigma.pdf')
dev.off()

# MSE:

stripchart_w_means(data$softmaxErr, data$group, 'MSE', c(0,0.3),2)
dev.copy(pdf,width=4, height = 4.5, 'plot_mse.pdf')
dev.off()


## kAE relationships

RegressPlot(data$kAE, log10(data$circles), data$group, legend=TRUE,
             'kAE', 'log stereo (arcsec)', c(0, 1), c(1, 3))
abline(h = 2.75, lty=2, col = 'black')
text(1, 2.81, 'nil', col='black', cex = 0.75)
dev.copy(pdf,width=4, height = 4.5, 'plot_kae_stereo_regress.pdf')
dev.off()

RegressPlot(data$kAE, data$LS, data$group, legend=FALSE,
             'kAE', 'balance point', c(0, 1), c(0.4, 1))
abline(h = 0.5, lty=2, col = 'black')
dev.copy(pdf,width=4, height = 4.5, 'plot_kae_balancept_regress.pdf')
dev.off()

### Supplementary Figure 2

# kAE
correlPlot(data$kAE, reddata$kAE, data$group, TRUE, 'original dataset', 'reduced dataset', c(0,1), c(0,1))
dev.copy(pdf,width=4, height = 4.5, 'plot_reduced_kae.pdf')
dev.off()

# sigma
correlPlot(data$sigma, reddata$sigma, data$group, FALSE, 'original dataset', 'reduced dataset', c(0.4,1.6), c(0.4,1.6))
dev.copy(pdf,width=4, height = 4.5, 'plot_reduced_sigma.pdf')
dev.off()

# uAE
correlPlot(data$uAEcln, reddata$uAE, data$group, FALSE, 'original dataset', 'reduced dataset', c(0,10), c(0,10))
dev.copy(pdf,width=4, height = 4.5, 'plot_reduced_uae.pdf')
dev.off()

# uFE
correlPlot(data$uFE, reddata$uFE, data$group, FALSE, 'original dataset', 'reduced dataset', c(0,3), c(0,3))
dev.copy(pdf,width=4, height = 4.5, 'plot_reduced_ufe.pdf')
dev.off()