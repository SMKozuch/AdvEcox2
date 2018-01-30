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
hist(residuals, breaks = 20)


summary(panel_regression, vcov = vcovHC)
