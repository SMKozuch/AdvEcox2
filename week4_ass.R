 #week4 ass
library(tidyverse)
library(foreign)
library(plm)
library(lmtest)

setwd("C:/Users/Samuel/OneDrive/UvA/S01P03/Advanced Econometrics II/Assignments")

fig <- read.dta('FIG.dta')

##### Panel regression
fig_panel <- pdata.frame(fig, index=c('country', 'period'))

normal_formula = formula(growth ~ (lly + btot + privo + dby) + log(initial) + log(initial) * (lly + btot + privo + dby))
extended_formula = formula(growth ~ (lly + btot + privo + dby) + log(initial) + log(initial) * (lly + btot + privo + dby) +
                             school + gov + log1p(pi) + log1p(bmp))

panel_regression <- plm(normal_formula, data=fig_panel)
summary(panel_regression)

panel_regression_control <- plm(extended_formula, data=fig_panel)
summary(panel_regression_control)
