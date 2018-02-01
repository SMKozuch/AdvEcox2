#week4 ass
#load packages
library(tidyverse)
library(foreign)
library(plm)
library(lmtest)

#setwd("C:/Users/Samuel/OneDrive/UvA/S01P03/Advanced Econometrics II/Assignments")

#load the data, 
fig <- read.dta('FIG.dta')

##### Panel regression
#load the data in the data frame, set the main axes (id which values correspond to n and which correspond to t)
fig_panel <- pdata.frame(fig, index=c('country', 'period'))

#specify the formula (w/ & w/0 control variables)
normal_formula = formula(growth ~ privo + log(initial) + log(initial) * (lly + btot + privo + dby))
extended_formula = formula(growth ~ privo + log(initial) + log(initial) * (lly + btot + privo + dby) +
                             school + gov + log1p(pi) + log1p(bmp))

#panel OLS w/o control vars
panel_regression <- plm(normal_formula, data=fig_panel, effect = 'twoways')
summary(panel_regression)

#panel OLS w/ control vars
panel_regression_control <- plm(extended_formula, data=fig_panel, effect = 'twoways')
summary(panel_regression_control)

##GMM
##gmm formula w/0 control vars
formula_gpmm_core <- as.formula(growth ~ lag(growth, 1:2) + privo
                                + log(initial) 
                                + lly * log(initial) + btot * log(initial) 
                                + privo * log(initial)
                                + dby * log(initial) | lag(growth, 1:5))

#panel GMM w/o control, (trans. = 'd' gives First Differences, model = 'twosteps' specifies two step variance estimation, 
#robust = True specifies robust standard errors)
pgmm_regression <- pgmm(formula_gpmm_core, data = fig_panel, transformation = 'd', model = 'twosteps',
                        robust = T)

summary(pgmm_regression)

#panel GMM w/ control
formula_gpmm_core_2 <- as.formula(growth ~ lag(growth, 1:2) + privo
                                  + log(initial) 
                                  + lly * log(initial) + btot * log(initial) 
                                  + privo * log(initial)
                                  + dby * log(initial) + school
                                  + log(1+pi) + gov + log(1 + bmp)| lag(growth, 1:5))

#panel GMM w/ control variables, (trans. = 'd' gives First Differences, model = 'twosteps' specifies two step variance estimation, 
#robust = True specifies robust standard errors)
pgmm_regression_2 <- pgmm(formula_gpmm_core_2, data = fig_panel, transformation = 'd', model = 'twosteps',
                          robust = T)


summary(pgmm_regression)
summary(pgmm_regression_2)
