#week4 ass
library(tidyverse)
library(foreign)
library(plm)
library(lmtest)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

fig <- read.dta('FIG.dta')

##### Panel regression
fig_panel <- pdata.frame(fig, index=c('country', 'period'))
panel_regression <- plm(growth ~ (lly + btot + privo) + dby + log(initial) + log(initial) * lly + log(initial) * btot + 
                          log(initial) * privo + log(initial) * dby, data = fig_panel)
summary(panel_regression)

panel_regression <- plm(growth ~ (lly + btot + privo) + dby + log(initial) + log(initial) * lly + log(initial) * btot + 
                          log(initial) * privo + log(initial) * dby + school + gov + log1p(pi) + log1p(bmp), data = fig_panel)

summary(panel_regression)

attributes(panel_regression)


residuals <- panel_regression$residuals
plot(residuals, type='l')

formula_gpmm_core <- as.formula(growth ~ lag(growth, 1:2) + lly + btot + privo
                                + dby + log(initial) 
                                + lly * log(initial) + btot * log(initial) 
                                + privo * log(initial)
                                + dby * log(initial) | lag(growth, 1:5))

pgmm_regression <- pgmm(formula_gpmm_core, data = fig_panel, transformation = 'd', model = 'twosteps',
     robust = T)

summary(pgmm_regression)

formula_gpmm_core_2 <- as.formula(growth ~ lag(growth, 1:2) + lly + btot + privo
                                + dby + log(initial) 
                                + lly * log(initial) + btot * log(initial) 
                                + privo * log(initial)
                                + dby * log(initial) + school
                                + log(1+pi) + gov + log(1 + bmp)| lag(growth, 1:5))

pgmm_regression_2 <- pgmm(formula_gpmm_core_2, data = fig_panel, transformation = 'd', model = 'twosteps',
                        robust = T)

summary(pgmm_regression_2)
