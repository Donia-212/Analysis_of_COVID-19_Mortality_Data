rm(list=ls())
library(MASS)
library(stats4)
library(fitdistrplus)
library(tidyverse)
library(lubridate)
library(patchwork)
library(VGAM)
library(VGAMextra)


################# India Wave 3 Dataset #####################

india_wave_3 <-  unlist(read.csv('Ind_Wave_3.csv'))

# Visualizing the data:
ggplot(as_tibble(india_wave_3), aes(value)) +
  geom_histogram(bins = 30, color = 'black', fill = 'grey') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Frequency') +
  labs(title = "India - Third Wave")


# 5 data summary 
summary(india_wave_3)
var(india_wave_3)




################# Investigating the data #####################

# The Third wave: from 11/2021 to 3/2022  
date <- dmy('1/11/2021')

as_tibble(india_wave_3) %>%
  mutate(days = seq(1,length(india_wave_3)), date = date + days(days)) %>%
  ggplot(aes(date,value)) +
  geom_line() +
  labs(title = "India - Third Wave", y = 'Number of Cases', x = 'Date')

as_tibble(india_wave_3) %>%
  mutate(days = seq(1,length(india_wave_3)), date = date + days(days)) %>%
  arrange(desc(value)) %>%
  View()


# We have 2 outliers in our data.
# Around the beginning of dec (2021-12-06), it can be related to Diwali (2021-12-09). 
# Most likely people gathered which lead to the spread of the virus.
# Another reason is that the christmas holiday and new years are also approcahing so people can be preparing for the holidays.
# Another outlier occured towards the end of March (2022-03-26) as there was the Rama Navami holiday.


in3 <- as_tibble(india_wave_3) %>%
  filter(value != 0, value <= 1500)  # 96.69% of the data remained

summary(in3)
var(in3)

in3 %>%
  mutate(days = seq(1,length(unlist(in3))), date = date + days(days)) %>%
  ggplot(aes(date,value)) +
  geom_line()+
  labs(title = "India - Third Wave - Modified", y = 'Number of Cases', x = 'Date')




################# Exponential Distribution #####################

# Maximum-likelihood fitting of univariate distribution
Exp <- fitdistr(unlist(in3), "exponential")
est_lambda <- Exp$estimate
l_exp <- Exp$loglik


# Exponential Distribution
ggplot(in3, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in3),dexp(unlist(in3),est_lambda)),se = F, color = 'red') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Exponential Distribution")




################# Exponentiated Exponential Distribution #####################

fit_ee <- vglm(value ~ 1, fam = expexpff, data = in3, trace = TRUE, maxit = 99)
parameters <-  as.vector(Coef(fit_ee))
alph_hat <- parameters[2]
lamb_hat <- parameters[1]

# Extracting Log-Likelihood
l_ee <-  logLik(fit_ee)


# Exponentiated exponential distribution
eed <-  lamb_hat*alph_hat * (1 - exp(-1 * lamb_hat * unlist(in3)))^(alph_hat - 1) * exp(-1 * lamb_hat * unlist(in3))
ggplot(in3, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in3),eed),se = F, color = 'orange') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Exponentiated Exponential Distribution")




################# Gamma Distribution #####################

fit <- vglm(value ~ 1, fam = gammaR, data = in3, trace = TRUE)

a <- Coef(fit)[2]
b <- Coef(fit)[1]
l_ga <-  logLik(fit)

# Gamma Distribution
ggplot(in3, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in3),dgamma(unlist(in3),a,b)), color = 'green', se = F) +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Gamma Distribution")




################# Weibull Distribution #####################

fit_w <- vglm(value ~ 1, fam = weibullRff, data = in3, trace = TRUE, maxit = 99)
k <- Coef(fit_w)[2]
l <- Coef(fit_w)[1]
l_w <-  logLik(fit_w)


# Weibull distribution
ggplot(in3, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in3),dweibull(unlist(in3),k,l)),se = F, color = 'blue') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Weibull Distribution")





################# Lognormal Distribution #####################

fit_lg <- vglm(value ~ 1, fam = lognormal, data = in3, trace = TRUE, maxit = 99)
mu <- Coef(fit_lg)[1]
sd <- Coef(fit_lg)[2]
l_lg <-  logLik(fit_lg)


ggplot(in3, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in3),dlnorm(unlist(in3),mu,sd)),se = F, color = 'yellow') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Lognormal Distribution")




################# Generalized Poisson Distribution #####################

fit_genp <- vglm(value ~ 1, fam = genpoisson0, data = in3, trace = TRUE, maxit = 99)
th <- Coef(fit_genp)[1]
la <- Coef(fit_genp)[2]
l_genp <-  logLik(fit_genp)

ggplot(in3, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in3),dgenpois0(unlist(in3),th,la)),se = F, color = 'maroon') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Generalized Poisson Distribution")





################# Summary #####################

# Combined plot: Modified
ggplot(in3, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in3),dexp(unlist(in3),est_lambda), color = 'Exponential'),se = F) +
  geom_smooth(aes(unlist(in3),eed, color = 'Exponentiated Exponential'),se = F) +
  geom_smooth(aes(unlist(in3),dgamma(unlist(in3),a,b), color = 'Gamma'), se = F) +
  geom_smooth(aes(unlist(in3),dweibull(unlist(in3),k,l), color = "Weibull"),se = F) +
  geom_smooth(aes(unlist(in3),dlnorm(unlist(in3),mu,sd), color = 'Lognormal'),se = F) +
  geom_smooth(aes(unlist(in3),dgenpois0(unlist(in3),th,la), color = 'Generalized Poisson'),se = F) +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Combined Plot: Modified") +
  scale_color_discrete(name='Distribution:') +
  theme(legend.position = "bottom")


# Testing the function: Kolmogorov-Smirnov Tests
p1 <- ks.test(in3,"pexp",est_lambda)

FEE <- function(x,alp,lam){
  (1 - exp(-lam * x) )^alp
}

p2 <- ks.test(in3,"FEE",alph_hat,lamb_hat)
p3 <- ks.test(in3,"pgamma",a,b)
p4 <- ks.test(in3,"pweibull",k,l)
p5 <- ks.test(in3,"plnorm",mu,sd)
p6 <- ks.test(in3,"pgenpois0",th,la)


# Table of Comparison
Dist <-  c("Exponential","Exponentiated Exponential", "Gamma", "Weibull", "Lognormal", 'Generalized Poisson')
par1 <-  c(est_lambda,alph_hat,a, k,mu,th)
par2 <-  c(est_lambda[2],est_lambda,b,l,sd,la)
Likelihood <-  c(l_exp,l_ee, l_ga, l_w,l_lg,l_genp)
p_value <-  c(p1$p.value, p2$p.value, p3$p.value, p4$p.value,p5$p.value,p6$p.value)
cbind(Dist,par1,par2,Likelihood,p_value)

