rm(list=ls())
library(MASS)
library(stats4)
library(fitdistrplus)
library(tidyverse)
library(lubridate)
library(patchwork)
library(VGAM)
library(VGAMextra)


################# India Wave 2 Dataset #####################

india_wave_2 <-  unlist(read.csv('Ind_Wave_2.csv')) 
in2 <- as_tibble(india_wave_2) %>%
  filter(!is.na(value)) 

# Visualizing the data:
ggplot(as_tibble(in2), aes(value)) +
  geom_histogram(bins = 30, color = 'black', fill = 'grey') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Frequency') +
  labs(title = "India - Second Wave")

# 5 data summary 
summary(in2)
var(in2)




################# Investigating the data #####################

# The second wave: from 03/2021 to 10/2021  
date <- dmy('1/3/2021')

in2 %>%
  mutate(days = seq(1,length(unlist(in2))), date = date + days(days)) %>%
  ggplot(aes(date,value)) +
  geom_line()+
  labs(title = "India - Second Wave", y = 'Number of Cases', x = 'Date')


in2 %>%
  mutate(days = seq(1,length(unlist(in2))), date = date + days(days)) %>%
  arrange(desc(value)) %>%
  View()


# We have two outliers in our data.
# Towards the middle of June (2021-06-11), it can be related to Eid Al-Adha (2021-07-19). 
# Most likely people started to prepare for the Eid
# The other outlier is around the end of July (2021-07-20) which could be caused by the celebrations
# of Eid Al-Adha (2021-07-19). The people most likely gathered which lead to the spread of the virus.
# As part of the Muslim culture, younger members of the family visit the elderly and spend eid with
# them which explains the sudden spike of daily death cases after eid as elderly people as less immune
# to covid and without proper care, their chances of passing away due to its symptoms increase.


in22 <- in2 %>%
  filter(value != 2020, value < 3500)


in22 %>%
  mutate(days = seq(1,length(unlist(in22))), date = date + days(days)) %>%
  ggplot(aes(date,value)) +
  geom_line() +
  labs(title = "India - Second Wave - Modified", y = 'Number of Cases', x = 'Date')


summary(in22)
var(in22)




################# Exponential Distribution #####################

# Maximum-likelihood fitting of univariate distribution
Exp <- fitdistr(unlist(in22), "exponential")
est_lambda <- Exp$estimate
l_exp <- Exp$loglik


# Exponential Distribution
ggplot(in22, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in22),dexp(unlist(in22),est_lambda)),se = F, color = 'red') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Exponential Distribution")




################# Exponentiated Exponential Distribution #####################

fit_ee <- vglm(value ~ 1, fam = expexpff, data = in22, trace = TRUE, maxit = 99)
parameters <-  as.vector(Coef(fit_ee))
alph_hat <- parameters[2]
lamb_hat <- parameters[1]

# Extracting Log-Likelihood
l_ee <-  logLik(fit_ee)


# Exponentiated exponential distribution
eed <-  lamb_hat*alph_hat * (1 - exp(-1 * lamb_hat * unlist(in22)))^(alph_hat - 1) * exp(-1 * lamb_hat * unlist(in22))
ggplot(in22, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in22),eed),se = F, color = 'orange') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Exponentiated Exponential Distribution")




################# Gamma Distribution #####################

fit <- vglm(value ~ 1, fam = gammaR, data = in22, trace = TRUE)

a <- Coef(fit)[2]
b <- Coef(fit)[1]
l_ga <-  logLik(fit)

# Gamma Distribution
ggplot(in22, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in22),dgamma(unlist(in22),a,b)), color = 'green', se = F) +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Gamma Distribution")




################# Weibull Distribution #####################

fit_w <- vglm(value ~ 1, fam = weibullRff, data = in22, trace = TRUE, maxit = 99)
k <- Coef(fit_w)[2]
l <- Coef(fit_w)[1]
l_w <-  logLik(fit_w)


# Weibull distribution
ggplot(in22, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in22),dweibull(unlist(in22),k,l)),se = F, color = 'blue') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Weibull Distribution")





################# Lognormal Distribution #####################

fit_lg <- vglm(value ~ 1, fam = lognormal, data = in22, trace = TRUE, maxit = 99)
mu <- Coef(fit_lg)[1]
sd <- Coef(fit_lg)[2]
l_lg <-  logLik(fit_lg)


ggplot(in22, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in22),dlnorm(unlist(in22),mu,sd)),se = F, color = 'yellow') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Lognormal Distribution")




################# Generalized Poisson Distribution #####################

fit_genp <- vglm(value ~ 1, fam = genpoisson0, data = in22, trace = TRUE, maxit = 99)
th <- Coef(fit_genp)[1]
la <- Coef(fit_genp)[2]
l_genp <-  logLik(fit_genp)

ggplot(in22, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in22),dgenpois0(unlist(in22),th,la)),se = F, color = 'maroon') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Generalized Poisson Distribution")





################# Summary #####################

# Combined plot: Modified
ggplot(in22, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in22),dexp(unlist(in22),est_lambda), color = 'Exponential'),se = F) +
  geom_smooth(aes(unlist(in22),eed, color = 'Exponentiated Exponential'),se = F) +
  geom_smooth(aes(unlist(in22),dgamma(unlist(in22),a,b), color = 'Gamma'), se = F) +
  geom_smooth(aes(unlist(in22),dweibull(unlist(in22),k,l), color = "Weibull"),se = F) +
  geom_smooth(aes(unlist(in22),dlnorm(unlist(in22),mu,sd), color = 'Lognormal'),se = F) +
  geom_smooth(aes(unlist(in22),dgenpois0(unlist(in22),th,la), color = 'Generalized Poisson'),se = F) +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Combined Plot: Modified") +
  scale_color_discrete(name='Distribution:') +
  theme(legend.position = "bottom")


# Testing the function: Kolmogorov-Smirnov Tests
p1 <- ks.test(in22,"pexp",est_lambda)

FEE <- function(x,alp,lam){
  (1 - exp(-lam * x) )^alp
}

p2 <- ks.test(in22,"FEE",alph_hat,lamb_hat)
p3 <- ks.test(in22,"pgamma",a,b)
p4 <- ks.test(in22,"pweibull",k,l)
p5 <- ks.test(in22,"plnorm",mu,sd)
p6 <- ks.test(in22,"pgenpois0",th,la)


# Table of Comparison
Dist <-  c("Exponential","Exponentiated Exponential", "Gamma", "Weibull", "Lognormal", 'Generalized Poisson')
par1 <-  c(est_lambda,alph_hat,a, k,mu,th)
par2 <-  c(est_lambda[2],est_lambda,b,l,sd,la)
Likelihood <-  c(l_exp,l_ee, l_ga, l_w,l_lg,l_genp)
p_value <-  c(p1$p.value, p2$p.value, p3$p.value, p4$p.value,p5$p.value,p6$p.value)
cbind(Dist,par1,par2,Likelihood,p_value)
