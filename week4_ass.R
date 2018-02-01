 #week4 ass
library(tidyverse)
library(foreign)
library(plm)
library(lmtest)

setwd("C:/Users/Samuel/OneDrive/UvA/S01P03/Advanced Econometrics II/Assignments")

fig <- read.dta('FIG.dta')

##### Panel regression
fig_panel <- pdata.frame(fig, index=c('country', 'period'))

normal_formula = formula(growth ~ privo + log(initial) + log(initial) * (lly + btot + privo + dby))
extended_formula = formula(growth ~ privo + log(initial) + log(initial) * (lly + btot + privo + dby) +
                             school + gov + log1p(pi) + log1p(bmp))

panel_regression <- plm(normal_formula, data=fig_panel, effect = 'twoways')
summary(panel_regression)

panel_regression_control <- plm(extended_formula, data=fig_panel, effect = 'twoways')
summary(panel_regression_control)

residuals <- panel_regression_control$residuals
plot(residuals, type='l')

formula_gpmm_core <- as.formula(growth ~ lag(growth, 1:2) + privo
                                + log(initial) 
                                + lly * log(initial) + btot * log(initial) 
                                + privo * log(initial)
                                + dby * log(initial) | lag(growth, 1:5))

pgmm_regression <- pgmm(formula_gpmm_core, data = fig_panel, transformation = 'd', model = 'twosteps',
                        robust = T)

summary(pgmm_regression)

formula_gpmm_core_2 <- as.formula(growth ~ lag(growth, 1:2) + privo
                                  + log(initial) 
                                  + lly * log(initial) + btot * log(initial) 
                                  + privo * log(initial)
                                  + dby * log(initial) + school
                                  + log(1+pi) + gov + log(1 + bmp)| lag(growth, 1:5))

pgmm_regression_2 <- pgmm(formula_gpmm_core_2, data = fig_panel, transformation = 'd', model = 'twosteps',
                          robust = T)

summary(pgmm_regression)
summary(pgmm_regression_2)
