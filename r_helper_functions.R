## Contains commonly-used functions for R analyses on model fits


map2amb <- function(leCol, reCol, aeIndex) {
  # input: left, right eye column of data
  #        aeIndex should be 'left' or 'right'
  # list out: amblyopic eye column, fellow eye column
  le_index <- aeIndex == 'left'
  re_index <- aeIndex == 'right'
  
  aeCol <- rep(NA, length(aeIndex))
  feCol <- rep(NA, length(aeIndex))
  
  aeCol[le_index] <- leCol[le_index]
  aeCol[re_index] <- reCol[re_index] 
  
  feCol[le_index] <- reCol[le_index] 
  feCol[re_index] <- leCol[re_index]
  
  return(as.data.frame(cbind(aeCol, feCol)))
}


## -- Printing summary tables/stats -- ##
printMeansAndSDs <- function(dvColumn, groupColumn, roundto) {
  # groupColumn must be factor with levels 'NS', 'AM', 'BD'
  # list out:
  # 1: control mean
  # 2: amblyopia mean
  # 3: non-amblyopic binocular disorder mean
  m_AM <- mean(dvColumn[groupColumn=='AM'], na.rm=T)
  sd_AM <- sd(dvColumn[groupColumn=='AM'], na.rm=T)
  m_BD <- mean(dvColumn[groupColumn=='BD'], na.rm=T)
  sd_BD <- sd(dvColumn[groupColumn=='BD'], na.rm=T)
  m_NS <- mean(dvColumn[groupColumn=='NS'], na.rm=T)
  sd_NS <- sd(dvColumn[groupColumn=='NS'], na.rm=T)
  
  output1 <- paste0('Control: ', round(m_NS,roundto), ' (', round(sd_NS,roundto), ')')
  output2 <- paste0('Amblyopia: ', round(m_AM,roundto), ' (', round(sd_AM,roundto), ')')
  output3 <- paste0('Strabismus with equal acuity: ', round(m_BD,2), ' (', round(sd_BD,roundto), ')')
  return(list(output1,output2,output3))
}
printTTestResult <- function(dvColumn, groupColumn) {
  # groupColumn must be factor with levels 'NS', 'AM' (ignores any 'BD')
  # list out:
  # 1: formatted t-test result with glass delta effect size
  # 2: p-value in string
  t_result <- t.test(dvColumn[groupColumn=='AM'], dvColumn[groupColumn=='NS'])
  gd <- glass_delta(dvColumn[groupColumn=='AM'], dvColumn[groupColumn=='NS'])
  
  if (t_result$p.value < 0.0001) { pstring <- 'p < 0.0001'
  } else{pstring <- paste('p =', signif(t_result$p.value,2))
  }
  
  out1 <- paste0("t(", round(t_result$parameter,1),
                 ") = ", abs(round(t_result$statistic,2)),
                 ", ", pstring,  
                 "; Glass's ∆ = ", abs(round(gd$Glass_delta,2)))
  return(list(out1, pstring))
}
printAnovaByEyeTable <- function(data, selectCols) {
  #  example selectCols = c('uAE', 'uFE')
  # assumes your condition/group var is called group
  
  data_subset <- subset(data, select = c('sID', 'group', selectCols))
  
  #  subset out the BD data
  data_subset <- subset(data_subset, group != 'BD')
  data_subset$group <- droplevels(data_subset$group)
  
  data_long <- gather(data_subset, eye, value, 3:4)
  data_long$sID <- factor(data_long$sID)
  data_long$eye <- factor(data_long$eye)
  
  res.aov <- anova_test(data=data_long, dv=value, wid = sID, between = group, within = eye)
  tmp <- get_anova_table(res.aov)
  
  posthoc1 <- data_long %>% group_by(eye) %>%
    anova_test(dv = value, wid = sID, between = group) %>%
    get_anova_table() %>%
    adjust_pvalue(method = 'bonferroni')
  
  indGrp <- which(tmp[,'Effect'] == 'group')
  indEye <- which(tmp[,'Effect'] == 'eye')
  indIxn <- which(tmp[,'Effect'] == 'group:eye')
  
  
  # build text output
  if (tmp[indGrp,'p'] < 0.05) {
    # sig
    strGrp <- 'Main effect of group'
  } else {
    # non-sig
    strGrp <- 'No main effect of group'
  }
  
  if (tmp[indGrp,'p'] < 0.0001) {
    pGrp <- 'p < 0.0001'
  } else {
    pGrp <- paste('p =', signif(tmp[indGrp,'p'],2))
  }
  
  if (tmp[indEye,'p'] < 0.05) {
    # sig
    strEye <- 'main effect of eye'
  } else {
    # non-sig
    strEye <- 'no main effect of eye'
  }
  
  if (tmp[indEye,'p'] < 0.0001) {
    pEye <- 'p < 0.0001'
  } else {
    pEye <- paste('p =', signif(tmp[indEye,'p'],2))
  }
  
  
  if (tmp[indIxn,'p'] < 0.05) {
    # sig
    strIxn <- 'a significant group by eye interaction'
  } else {
    # non-sig
    strIxn <- 'no group by eye interaction'
  }
  
  if (tmp[indIxn,'p'] < 0.0001) {
    pIxn <- 'p < 0.0001'
  } else {
    pIxn <- paste('p =', signif(tmp[indIxn,'p'],2))
  }
  
  
  # posthoc: eye 1
  if (posthoc1$p[1] < 0.05) {
    # sig
    strPh1 <- paste('difference in', selectCols[1])
  } else {
    # non-sig
    strPh1 <- paste('no difference in', selectCols[1])
  }
  
  if (posthoc1$p[1] < 0.0001) {
    pPh1 <- 'p < 0.0001'
  } else {
    pPh1 <- paste('p =', signif(posthoc1$p[1],2))
  }
  
  
  # posthoc: eye 2
  if (posthoc1$p[2] < 0.05) {
    # sig
    strPh2 <- paste('difference in', selectCols[2])
  } else {
    # non-sig
    strPh2 <- paste('no difference in', selectCols[2])
  }
  
  if (posthoc1$p[2] < 0.0001) {
    pPh2 <- 'p < 0.0001'
  } else {
    pPh2 <- paste('p =', signif(posthoc1$p[2],2))
  }
  
  
  out1 <- paste0(strGrp, ", F(", tmp[indGrp,'DFn'], ",", tmp[indGrp,'DFd'],
                 ") = ", round(tmp[indGrp,'F'],2), ", ", pGrp, ", ges = ", signif(tmp[indGrp, 'ges'],2))
  
  out2 <- paste0(strEye, ", F(", tmp[indEye,'DFn'], ",", tmp[indEye,'DFd'],
                 ") = ", round(tmp[indEye,'F'],2), ", ", pEye,
                 ", ges = ", signif(tmp[indEye, 'ges'],2))
  
  out3 <- paste0(strIxn, ", F(", tmp[indIxn,'DFn'], ",", tmp[indIxn,'DFd'],
                 ") = ", round(tmp[indIxn,'F'],2), ", ", pIxn, ", ges = ", signif(tmp[indIxn, 'ges'],2))
  
  out4 <- paste0("Simple main effects: ", strPh1, 
                 ", F(", posthoc1[1,'DFn'], ",", posthoc1[1,'DFd'],
                 ") = ", round(posthoc1[1,'F'],2), ", ", pPh1,
                 ", ges = ", signif(posthoc1[1, 'ges'],2))
  
  
  out5 <- paste0("Simple main effects: ", strPh2, ", F(", posthoc1[2,'DFn'], ",", posthoc1[2,'DFd'],
                 ") = ", round(posthoc1[2,'F'],2), ", ", pPh2,
                 ", ges = ", signif(posthoc1[2, 'ges'],2))
  
  out6 <- paste0("F(", tmp[indGrp,'DFn'], ",", tmp[indGrp,'DFd'],
                 ") = ", round(tmp[indGrp,'F'],2), ", ", pGrp)
  
  out7 <- paste0("F(", tmp[indEye,'DFn'], ",", tmp[indEye,'DFd'],
                 ") = ", round(tmp[indEye,'F'],2), ", ", pEye)
  
  out8 <- paste0("F(", tmp[indIxn,'DFn'], ",", tmp[indIxn,'DFd'],
                 ") = ", round(tmp[indIxn,'F'],2), ", ", pIxn)
  
  
  returnStr <- paste(out1,out2,out3,out4,out5)
  return(list(returnStr, pGrp, pEye, pIxn, pPh1, pPh2, out6, out7, out8))
}
printRegressionResults <-function(dvCol,groupCol,modelParamCol){
  # output order: group alone, parameter alone, parameter after controlling for grp
  mod_group <- lm(dvCol ~ groupCol)
  mod_param <- lm(dvCol ~ modelParamCol)
  mod_group_param <- lm(dvCol ~ groupCol + modelParamCol)
  
  smod_group <- summary(mod_group)
  smod_param <- summary(mod_param)
  smod_group_param <- Anova(mod_group_param)
  
  # format pstring
  
  p_group <-   pf(smod_group$fstatistic[1],
                  smod_group$fstatistic[2],
                  smod_group$fstatistic[3], lower.tail = F)
  if (p_group < 0.0001) {
    p1 <- 'p < 0.0001'
  } else {
    p1 <- paste('p =', signif(p_group,2))
  }
  if (p_group < 0.05) {
    p1 <- paste(p1, '*')
  }
  
  
  p_param <-   pf(smod_param$fstatistic[1],
                  smod_param$fstatistic[2],
                  smod_param$fstatistic[3], lower.tail = F)
  if (p_param < 0.0001) {
    p2 <- 'p < 0.0001'
  } else {
    p2 <- paste('p =', signif(p_param,2))
  }
  if (p_param < 0.05) {
    p2 <- paste(p2, '*')
  }
  
  
  
  p_group_param <- smod_group_param$P[2]
  
  if (p_group_param < 0.0001) {
    p3 <- 'p < 0.0001'
  } else {
    p3 <- paste('p =', signif(p_group_param,2))
  }
  if (p_group_param < 0.05) {
    p3 <- paste(p3, '*')
  }
  
  
  out1 <- paste0("F(", smod_group$fstatistic["numdf"], ",",
                 smod_group$fstatistic["dendf"], 
                 ") = ", round(smod_group$fstatistic['value'],2),
                 ", ", p1)
  out2 <- paste0("F(", smod_param$fstatistic["numdf"],",",
                 smod_param$fstatistic["dendf"], 
                 ") = ", round(smod_param$fstatistic['value'],2),
                 ", ", p2)
  
  out3 <-  paste0("F(", smod_group_param$Df[2],",",
                  smod_group_param$Df[3], 
                  ") = ", round(smod_group_param$F[2],2),
                  ", ", p3)
  return(list(out1,out2,out3))
}
printModelAnovaTable <- function(data, selectCols) {
  #  compares models
  #  label out is always second model listed
  #  example selectCols = c('softmax', 'meanJoystick')
  #  assumes your condition/group var is called group
  data_subset <- subset(data, select = c('sID', 'group', selectCols))
  
  data_long <- gather(data_subset, modeltype, mse, 3:4)
  data_long$sID <- factor(data_long$sID)
  data_long$modeltype <- factor(data_long$modeltype)
  
  # lmer/lme4/lmerTest notes:
  # call would be 
  # mod <- lmer(mse ~ modeltype + (1|sID) + (1|group) + (1|modeltype:group), data=data_long)
  # fixed <- anova(mod)
  # random <- ranova(mod)
  # however not sure how to get F-tests for random (i.e. ranova() does not provide them)
  
  
  res.aov <- anova_test(data=data_long, dv=mse, wid = sID, between = group, within = modeltype)
  tmp <- get_anova_table(res.aov)
  
  indGrp <- which(tmp[,'Effect'] == 'group')
  indMod <- which(tmp[,'Effect'] == 'modeltype')
  indIxn <- which(tmp[,'Effect'] == 'group:modeltype')
  
  if (tmp[indGrp,'p'] < 0.0001) {
    pGrp <- 'p < 0.0001'
  } else {
    pGrp <- paste('p =', signif(tmp[indGrp,'p'],2))
  }
  
  if (tmp[indMod,'p'] < 0.0001) {
    pMod <- 'p < 0.0001'
  } else {
    pMod <- paste('p =', signif(tmp[indMod,'p'],2))
  }
  
  if (tmp[indIxn,'p'] < 0.0001) {
    pIxn <- 'p < 0.0001'
  } else {
    pIxn <- paste('p =', signif(tmp[indIxn,'p'],2))
  }
  
  out1 <- paste0("F(", tmp[indGrp,'DFn'], ",", tmp[indGrp,'DFd'],
                 ") = ", round(tmp[indGrp,'F'],2), ", ", pGrp)
  
  out2 <- paste0("F(", tmp[indMod,'DFn'], ",", tmp[indMod,'DFd'],
                 ") = ", round(tmp[indMod,'F'],2), ", ", pMod)
  
  out3 <- paste0("F(", tmp[indIxn,'DFn'], ",", tmp[indIxn,'DFd'],
                 ") = ", round(tmp[indIxn,'F'],2), ", ", pIxn)
  
  dfout <- data.frame(model=out2, group=out1, interaction=out3)
  rownames(dfout) <- selectCols[2]
  return(dfout)
}


## -- Plots -- ##

stripchart_w_means <- function(DV, group, strDV, ylimits, roundto) {
  # DV = column for DV (eg data$delay)
  # group = column for group (eg data$group)
  # strDV = string, DV name (for plot)
  # ylimits = c(min, max)
  # roundto = for plotting, nearest (helps with stacking)
  
  # deal with missing values
  if (any(is.na(DV))) {
    ind <- !is.na(DV)
    DV <- DV[ind]
    group <- group[ind]
  }
  
  desc_matrix <- groupMeans(DV, group) # helper function below
  
  # rounding for stripcharts - allows prettier stacking
  DV <- round(DV,roundto) 
  
  # plot
  stripchart(DV ~ group,
             group.names = c(' ', ' ', ' '),
             xlab = ' ',
             ylab = paste(strDV, '[SE]'),
             xlim = c(0.5,3.5),
             ylim = ylimits,
             method = 'stack', offset = 0.3,
             pch = c(pchAM, pchBD, pchNS),
             col = c(colAM, colBD, colNS),
             bg = grey_transp, vertical=TRUE,
             cex.axis=1.2, cex.lab = 1.2) 
  
  grid(nx = NA, ny = NULL)
  
  offset <- 0.2
  drawlinesX <- c(1, 1, NA, 2, 2, NA, 3, 3) - offset
  drawlinesY <- c(desc_matrix[1,1] + desc_matrix[1,2],
                  desc_matrix[1,1] - desc_matrix[1,2], NA,
                  desc_matrix[2,1] + desc_matrix[2,2],
                  desc_matrix[2,1] - desc_matrix[2,2], NA,
                  desc_matrix[3,1] + desc_matrix[3,2],
                  desc_matrix[3,1] - desc_matrix[3,2])
  
  lines(drawlinesX, drawlinesY, col="black", lwd=3)
  
  points(c(1,2,3)-offset , desc_matrix[,1], 
         pch = c(pchAM, pchBD, pchNS), 
         col = c(colAM_solid, colBD_solid, colNS_solid),
         cex = 1.5)
  
  mtext('amblyopia', side = 1, at=1, font=1, line = 1, cex=.9, col=colAM_solid)
  mtext('strabismus', side = 1, at=2, font=1, line=1,cex=.9, col=colBD_solid)
  mtext('control', side = 1, at=3, font=1, line=1,cex=.9, col=colNS_solid)
}
stripchart_w_means_byeye <- function(DV_AE, DV_FE, group_AE, group_FE, strDV, ylimits) {
  # DV_AE = column for DV, amblyopic eye (eg data$uAE)
  # DV_FE = column for DV, fellow eye (eg data$uFE)
  # group = column for group (eg data$group)
  # strDV = string, DV name (for plot)
  # ylimits = c(min, max)
  
  # missing data checks
  
  if (any(is.na(DV_AE))) {
    ind <- !is.na(DV_AE)
    DV_AE <- DV_AE[ind]
    group_AE <- group_AE[ind]
  }
  
  if (any(is.na(DV_FE))) {
    ind <- !is.na(DV_FE)
    DV_FE <- DV_FE[ind]
    group_FE <- group_FE[ind]
  }
  
  group_AE <- paste0(group_AE, ' (AE)')
  group_FE <- paste0(group_FE, ' (FE)')
  
  desc_matrix <- groupMeans_byEye(c(DV_AE, DV_FE), c(group_AE, group_FE))
  xvar <- round(c(DV_AE, DV_FE),4)

  # plot
  stripchart(xvar ~ c(group_AE, group_FE),
             group.names = c('AE', 'FE', 'AE', 'FE', 'AE', 'FE'),
             xlab = ' ',
             ylab = paste(strDV, '[SE]'),
             xlim = c(0.5,6.5),
             ylim = ylimits,
             method = 'stack', offset = 0.3,
             pch = c(pchAM, pchAM, pchBD, pchBD, pchNS, pchNS),
             col = c(colAM, colAM, colBD, colBD, colNS, colNS),
             bg=grey_transp, vertical=TRUE,
             cex.axis=1.2, cex.lab=1.2) 
  
  grid(nx = NA, ny = NULL)
  
  offset <- 0.2
  drawlinesX <- c(1, 1, NA, 2, 2, NA, 
                  3, 3, NA, 4, 4, NA,
                  5, 5, NA, 6, 6) - offset
  drawlinesY <- c(desc_matrix[1,1] + desc_matrix[1,2],
                  desc_matrix[1,1] - desc_matrix[1,2], NA,
                  desc_matrix[2,1] + desc_matrix[2,2],
                  desc_matrix[2,1] - desc_matrix[2,2], NA,
                  desc_matrix[3,1] + desc_matrix[3,2],
                  desc_matrix[3,1] - desc_matrix[3,2], NA,
                  desc_matrix[4,1] + desc_matrix[4,2],
                  desc_matrix[4,1] - desc_matrix[4,2], NA,
                  desc_matrix[5,1] + desc_matrix[5,2],
                  desc_matrix[5,1] - desc_matrix[5,2], NA,
                  desc_matrix[6,1] + desc_matrix[6,2],
                  desc_matrix[6,1] - desc_matrix[6,2])
  
  lines(drawlinesX, drawlinesY, col="black", lwd=3)
  
  points(c(1,2,3,4,5,6)-offset , desc_matrix[,1],
         pch = c(pchAM, pchAM, pchBD, pchBD, pchNS, pchNS),
         col = c(colAM_solid, colAM_solid,
                 colBD_solid, colBD_solid,
                 colNS_solid, colNS_solid),
         cex = 1.2)
  
  
  mtext('amblyopia', side = 1, at=1.5, font=1, line = 2, cex=1, col=colAM_solid)
  mtext('strabismus', side = 1, at=3.5, font=1, line=2,cex=1, col=colBD_solid)
  mtext('control', side = 1, at=5.5, font=1, line=2,cex=1, col=colNS_solid)
}
correlPlot <- function(x, y, group, legendTF, xlab, ylab, xlim, ylim) {
  
  # Get coefficient for the plot
  mod_full <- lm(y ~ x)           
  xycor <- round(cor(x, y),2)
  
  # plot the x,y data
  plot(x,y,                 
       xlab = xlab,
       ylab = ylab,
       xlim = xlim,
       ylim = ylim,
       type = 'n',
       cex = 0.8,
       cex.lab = 0.8,
       cex.axis = 0.8
  )
  
  grid(nx = NULL, ny = NULL)
  
  # Add the unity line
  abline(a=0, b=1, col='black')
  
  # add points
  
  # add amblyopia
  points(x[group == 'AM'], y[group == 'AM'],
         pch = pchAM, col=colAM)
  
  # add strab
  points(x[group == 'BD'], y[group == 'BD'],
         pch = pchBD, col=colBD)
  
  # add control
  points(x[group == 'NS'], y[group == 'NS'],
         pch = pchNS, col=colNS)
  
  # add means
  t.am.x <- t.test(x[group == 'AM'])
  t.am.y <- t.test(y[group == 'AM'])
  t.bd.x <- t.test(x[group == 'BD'])
  t.bd.y <- t.test(y[group == 'BD'])
  t.ns.x <- t.test(x[group == 'NS'])
  t.ns.y <- t.test(y[group == 'NS'])
  
  m.am.x <- t.am.x$estimate
  m.am.y <- t.am.y$estimate
  m.bd.x <- t.bd.x$estimate
  m.bd.y <- t.bd.y$estimate
  m.ns.x <- t.ns.x$estimate
  m.ns.y <- t.ns.y$estimate
  
  se.am.x <- t.am.x$stderr
  se.am.y <- t.am.y$stderr
  se.bd.x <- t.bd.x$stderr
  se.bd.y <- t.bd.y$stderr
  se.ns.x <- t.ns.x$stderr
  se.ns.y <- t.ns.y$stderr
  
  drawlinesX <- c(m.am.x - se.am.x, m.am.x + se.am.x, NA,
                  m.bd.x - se.bd.x, m.bd.x + se.bd.x, NA,
                  m.ns.x - se.ns.x, m.ns.x + se.ns.x, NA,
                  m.am.x, m.am.x, NA,
                  m.bd.x, m.bd.x, NA,
                  m.ns.x, m.ns.x)
  
  drawlinesY <- c(m.am.y, m.am.y, NA,
                  m.bd.y, m.bd.y, NA,
                  m.ns.y, m.ns.y, NA,
                  m.am.y - se.am.y, m.am.y + se.am.y, NA,
                  m.bd.y - se.bd.y, m.bd.y + se.bd.y, NA,
                  m.ns.y - se.ns.y, m.ns.y + se.ns.y)
  
  lines(drawlinesX, drawlinesY, col='black', lwd=1.5)
  
  points(m.am.x, m.am.y, 
         pch = c(pchAM), 
         col = c(colAM_solid),
         cex = 1.5)
  
  points(m.bd.x, m.bd.y, 
         pch = c(pchBD), 
         col = c(colBD_solid),
         cex = 1.5)
  
  points(m.ns.x, m.ns.y, 
         pch = c(pchNS), 
         col = c(colNS_solid),
         cex = 1.5)
  
  if (legendTF == TRUE) {
    legendloc <- 'bottomright'
    legend(legendloc, legend=c('amblyopia', 'strabismus', 'control'),
           col=c(colAM_solid, colBD_solid, colNS_solid),
           pch=c(pchAM, pchBD, pchNS),
           bty='n',
           inset = c(0.05, 0),
           y.intersp = 0.7,
           lwd=c(1,1), lty=c(2,3,4), cex=0.7)
  }
}
RegressPlot <- function(x, y, group, legendTF, xlab, ylab, xlim, ylim) {
  
  group <- relevel(group, 'NS') # normally-sighted folks are the ref
  mod_all <-lm(y ~ group + x)
  coefs_all <- coef(mod_all)
  
  b_all <- coefs_all[4]
  a_NS <- coefs_all[1]
  a_AM <- coefs_all[1] + coefs_all[2]
  a_BD <- coefs_all[1] + coefs_all[3]
  
  # plot the x,y data
  plot(x,y,                 
       xlab = xlab,
       ylab = ylab,
       xlim = xlim,
       ylim = ylim,
       type = 'n',
       cex = 1,
       cex.lab = 1.2,
       cex.axis = 1
  )
  
  grid(nx = NULL, ny = NULL)
  
  # Add the correlation lines
  
  # amblyopia
  abline(a=a_AM, b = b_all, 
         lwd=1.5, 
         lty=2, 
         col=colAM_solid)   
  
  # strab
  abline(a=a_BD, b=b_all, 
         lwd=1.5, 
         lty=3, 
         col=colBD_solid)   
  
  # control
  abline(a = a_NS, b = b_all, 
         lwd=1.5, 
         lty=4, 
         col=colNS_solid)   
  
  # add points
  
  # add amblyopia
  points(x[group == 'AM'], y[group == 'AM'],
         pch = pchAM, col=colAM)
  
  # add strab
  points(x[group == 'BD'], y[group == 'BD'],
         pch = pchBD, col=colBD)
  
  # add control
  points(x[group == 'NS'], y[group == 'NS'],
         pch = pchNS, col=colNS)
  
  if (legendTF == TRUE) {
    legendloc <- 'bottomleft'
    legend(legendloc, legend=c('amblyopia', 'strabismus', 'control'),
           col=c(colAM_solid, colBD_solid, colNS_solid),
           pch=c(pchAM, pchBD, pchNS),
           bty='n',
           inset = c(0.05, 0),
           y.intersp = 0.7,
           lwd=c(1,1), lty=c(2,3,4), cex=0.7)
  }
}


## Helper functions
groupMeans <- function(DV, group) {
  # DV = column for DV (eg data$delay)
  # group = column for group (eg data$group)
  
  # get means and SEs
  t_AM <- t.test(DV[group == 'AM'])
  desc_AM <- c(t_AM$estimate, t_AM$stderr, t_AM$conf.int)
  
  t_BD <- t.test(DV[group == 'BD'])
  desc_BD <- c(t_BD$estimate, t_BD$stderr, t_BD$conf.int)
  
  # workaround (error if controls have all same [eg stereo] values)
  t_NS <- try(t.test(DV[group == 'NS']))
  desc_NS <- c(t_NS$estimate, t_NS$stderr, t_NS$conf.int)
  
  descriptives <- rbind(desc_AM, desc_BD, desc_NS)
  rownames(descriptives) <- c('AM', 'BD', 'NS')
  colnames(descriptives) <- c('mean', 'stderr', '95ci-low', '95ci-upp')
  
  return(descriptives)
}
groupMeans_byEye <- function(DV, group) {
  
  # get means and SEs
  t_AM_AE <- t.test(DV[group == 'AM (AE)'])
  desc_AM_AE <- c(t_AM_AE$estimate, t_AM_AE$stderr, t_AM_AE$conf.int)
  t_AM_FE <- t.test(DV[group == 'AM (FE)'])
  desc_AM_FE <- c(t_AM_FE$estimate, t_AM_FE$stderr, t_AM_FE$conf.int)
  
  
  # workaround for controls having stereo values that are all the same (good stereo) is needed
  
  t_NS_AE <- try(t.test(DV[group == 'NS (AE)']))
  desc_NS_AE <- c(t_NS_AE$estimate, t_NS_AE$stderr, t_NS_AE$conf.int)
  
  t_NS_FE <- try(t.test(DV[group == 'NS (FE)']))
  desc_NS_FE <- c(t_NS_FE$estimate, t_NS_FE$stderr, t_NS_FE$conf.int)
  
  
  # 
  t_BD_AE <- t.test(DV[group == 'BD (AE)'])
  desc_BD_AE <- c(t_BD_AE$estimate, t_BD_AE$stderr, t_BD_AE$conf.int)
  t_BD_FE <- t.test(DV[group == 'BD (FE)'])
  desc_BD_FE <- c(t_BD_FE$estimate, t_BD_FE$stderr, t_BD_FE$conf.int)
  
  
  descriptives <- rbind(desc_AM_AE, desc_AM_FE,
                        desc_BD_AE, desc_BD_FE,
                        desc_NS_AE, desc_NS_FE)
  rownames(descriptives) <- c('AM (AE)','AM (FE)', 
                              'BD (AE)', 'BD (FE)',
                              'NS (AE)', 'NS (FE)')
  colnames(descriptives) <- c('mean', 'stderr', '95ci-low', '95ci-upp')
  
  return(descriptives)
}
ceil_round <- function(x, nDigits) {
  # ceilings to the nearest n digits
  n <- 10^nDigits
  return(ceiling(x*n)/n)
}
floor_round <- function(x, nDigits) {
  # floors to the nearest n digits
  n <- 10^nDigits
  return(floor(x*n)/n)
}

## colours & pchs
grey_transp <- rgb(.1, .1, .1, .1)

colAM <- rgb(1,0,0, alpha = 0.6)
colNS <- rgb(0,0,1, alpha = 0.6)
colBD <- rgb(0,.5, 0, alpha = 0.6)
pchAM <- 16
pchNS <- 18
pchBD <- 15

colAM_solid <- rgb(1,0,0)
colNS_solid <- rgb(0,0,1)
colBD_solid <- rgb(0,.5,0)
