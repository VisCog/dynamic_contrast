library('psych')
library('tidyr')
library('rstatix')
library('lm.beta')
library('kableExtra')
library('effectsize')
library('readxl')

setwd("/Users/bear/Downloads/dynamic_contrast")

# read in data
allMotorData <- read.csv("/Users/bear/Downloads/dynamic_contrast/fitdata_vep_psychophysics/vep_psychophysics_fits.csv", header = TRUE)
allVepData <- read.csv("/Users/bear/Downloads/dynamic_contrast/fitdata_vep/vep_fits.csv", header = TRUE)

# drop BD
allMotorData <- subset(allMotorData, group != 'BD')
allVepData <- subset(allVepData, group != 'BD')
# 
# # these two have AE acuity better than 0.2 but no strab, talk to KTH on whether to include
# allMotorData <- subset(allMotorData, sID != 'AM_RE_XV_19')
# allMotorData <- subset(allMotorData, sID != 'AM_LE_QQ_37')
# 
# allVepData <- subset(allVepData, sID != 'AM_RE_XV_19')
# allVepData <- subset(allVepData, sID != 'AM_LE_QQ_37')

# factorize
allMotorData$group <- factor(allMotorData$group)
allVepData$group <- factor(allVepData$group)

# merge
allData <- merge(allMotorData, allVepData, 'sID', suffixes = c('.mtr', '.vep'))

# test: attenuation for controls vs amb, motor
t_result_motor <- t.test(allData$kAE.mtr[allData$group.mtr=='AM'], allData$kAE.mtr[allData$group.mtr=='NS'])
paste(
  'amblyopia M =', round(mean(allData$kAE.mtr[allData$group.mtr=='AM']),2), 
  'control M =', round(mean(allData$kAE.mtr[allData$group.mtr=='NS']),2))

paste('t(', abs(round(t_result_motor$parameter,1)),
      ') =', abs(round(t_result_motor$statistic,2)),
      'p =', abs(round(t_result_motor$p.value,3)))


# test: attenuation for controls vs amb, vep
t_result_vep <- t.test(allData$kAE.vep[allData$group.vep=='AM'], allData$kAE.vep[allData$group.mtr=='NS'])
paste(
  'amblyopia M =', round(mean(allData$kAE.vep[allData$group.mtr=='AM']),2), 
  'control M =', round(mean(allData$kAE.vep[allData$group.mtr=='NS']),2))

paste('t(', abs(round(t_result_vep$parameter,1)),
      ') =', abs(round(t_result_vep$statistic,2)),
      'p =', abs(round(t_result_vep$p.value,3)))




### COLORS AND PCHS

grey_transp <- rgb(.1, .1, .1, .1)

colAM <- rgb(1,0,0, alpha = 0.6)
colNS <- rgb(0,0,1, alpha = 0.6)
pchAM <- 16
pchNS <- 18

colAM_solid <- rgb(1,0,0)
colNS_solid <- rgb(0,0,1)


# Plot motor on x-axis, eeg on y=axis
x <- allData$kAE.mtr
y <- allData$kAE.vep

group <- allData$group.mtr
  # plot the x,y data
  plot(x,y,                 
       xlab = 'behavior',
       ylab = 'VEP signal',
       xlim = c(0,1),
       ylim = c(0,1),
       type = 'n',
       cex = 1,
       cex.lab = 1.2,
       cex.axis = 1.2
  )
  
  grid(nx = NULL, ny = NULL)
  
  # add points
  
  # add amblyopia
  points(x[group == 'AM'], y[group == 'AM'],
         pch = pchAM, col=colAM)
  
  # add control
  points(x[group == 'NS'], y[group == 'NS'],
         pch = pchNS, col=colNS)
  
  # add means
  t.am.x <- t.test(x[group == 'AM'])
  t.am.y <- t.test(y[group == 'AM'])
  t.ns.x <- t.test(x[group == 'NS'])
  t.ns.y <- t.test(y[group == 'NS'])
  
  m.am.x <- t.am.x$estimate
  m.am.y <- t.am.y$estimate
  m.ns.x <- t.ns.x$estimate
  m.ns.y <- t.ns.y$estimate
  
  se.am.x <- t.am.x$stderr
  se.am.y <- t.am.y$stderr
  se.ns.x <- t.ns.x$stderr
  se.ns.y <- t.ns.y$stderr
  # 
  # drawlinesX <- c(m.am.x - se.am.x, m.am.x + se.am.x, NA,
  #                 m.ns.x - se.ns.x, m.ns.x + se.ns.x, NA,
  #                 m.am.x, m.am.x, NA,
  #                 m.ns.x, m.ns.x)
  # 
  # drawlinesY <- c(m.am.y, m.am.y, NA,
  #                 m.ns.y, m.ns.y, NA,
  #                 m.am.y - se.am.y, m.am.y + se.am.y, NA,
  #                 m.ns.y - se.ns.y, m.ns.y + se.ns.y)
  # 
  # 
  # drawlinesX <- c(m.am.x - se.am.x, m.am.x + se.am.x, NA,
  #                 m.ns.x - se.ns.x, m.ns.x + se.ns.x, NA,
  #                 0, 0, NA,
  #                 0, 0)
  # 
  # drawlinesY <- c(0, 0, NA,
  #                 0, 0, NA,
  #                 m.am.y - se.am.y, m.am.y + se.am.y, NA,
  #                 m.ns.y - se.ns.y, m.ns.y + se.ns.y)
  # 
  # lines(drawlinesX, drawlinesY, col='black', lwd=1.5)
  # 
  # points(m.am.x, 0,
  #        pch = c(pchAM),
  #        col = c(colAM_solid),
  #        cex = 1.5)
  # 
  # points(m.ns.x, 0,
  #        pch = c(pchNS),
  #        col = c(colNS_solid),
  #        cex = 1.5)

 # lines(c(m.am.x,(m.am.x+m.ns.x)/2), c(0.025, 0.05),col='black', lwd=1.5)

  # points(0, m.ns.y,
  #        pch = c(pchNS),
  #        col = c(colNS_solid),
  #        cex = 1.5)
  # 
  # 
  # points(0, m.am.y,
  #        pch = c(pchAM),
  #        col = c(colAM_solid),
  #        cex = 1.5)

  #text(0,(m.am.y+m.ns.y)/2, paste0('p =', round(t_result_vep$p.value,3)), pos = 4, col='black', cex=0.7)
  #text((m.am.x+m.ns.x)/2, 0.06, paste0('p =', round(t_result_vep$p.value,3)), adj = 0.5, col='black', cex=0.7)
  #text((m.am.x+m.ns.x)/2, 0.02, '*', adj = 0.5, col='black', cex=1.5)
  
  lines(c(-1, 2), c(-1,2), lty=2)
  
  legend(x=c(0.1-.07,0.5-0.07), y=c(0.85, 0.65),
         legend=c('amblyopia', 'control'),
         col=c(colAM_solid, colNS_solid),
         #bty='n',
         y.intersp = 1,
         pch=c(pchAM, pchNS), cex=.85)
  
  
  dev.copy(pdf,width=4, height = 4.5, 'motor_vs_eeg2.pdf')
  dev.off()
