 #week4 ass
library(tidyverse)
library(foreign)
library(plm)
library(lmtest)

setwd("C:/Users/Samuel/OneDrive/UvA/S01P03/Advanced Econometrics II/Assignments")

fig <- read.dta('FIG.dta')

##### Panel regression
fig_panel <- pdata.frame(fig, index=c('country', 'period'))
panel_regression <- plm(growth ~ (lly + btot + privo) + dby + log(initial) + log(initial) * lly + log(initial) * btot + 
                          log(initial) * privo + log(initial) * dby, data=fig_panel)
summary(panel_regression)

panel_regression <- plm(growth ~ (lly + btot + privo) + dby + log(initial) + log(initial) * lly + log(initial) * btot + 
                          log(initial) * privo + log(initial) * dby + school + gov + log1p(pi) + log1p(bmp), data=fig_panel)

betas <- coef(panel_regression)

lly_fitted <- fig_panel[,'lly'] * betas[1]
btot_fitted <- fig_panel[,'btot'] * betas[2]
privo_fitted <- fig_panel[,'privo'] * betas[3]

growth_fitted <- lly_fitted + btot_fitted + privo_fitted

resettest(growth ~ log(lly) + log(btot) + log(privo) + log(dby), power = 2:5, data = fig_panel, type = c('fitted'))

