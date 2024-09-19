rm(list=ls())
library(MASS)
library(stats4)
library(VGAM)
library(tidyverse)
library(fitdistrplus)
library(lubridate)
library(patchwork)
library(VGAM)
library(VGAMextra)
library(EstimationTools)
library(bbmle)
library(psych)


################# US Wave 1 Dataset #####################

us_wave_1 <-  unlist(read.csv('US_Wave_1.csv'))

# Visualizing the data
ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(bins = 30, color = 'black', fill = 'grey') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Frequency') +
  labs(title = "US - First Wave")

# 5 data summary 
summary(us_wave_1)
var(us_wave_1)
(d <- describe(us_wave_1))
(coef <- d$sd/d$mean)




################# Investigating the data #####################

# The first wave: from 4/2020 to 8/2020  
date <- dmy('1/4/2020')

as_tibble(us_wave_1) %>%
  mutate(days = seq(1,length(us_wave_1)), date = date + days(days)) %>%
  arrange(desc(value)) %>%
  View()

as_tibble(us_wave_1) %>%
  mutate(days = seq(1,length(us_wave_1)), date = date + days(days)) %>%
  ggplot(aes(date,value)) +
  geom_line()

# Why is the daily death cases gradually decreasing? Why did it increase again?
# 1st of Apr was the first Day of Stay at Home Order (COVID-19 Lockdown).
# 1st of Jun was the last Day of Stay at Home Order	(Michigan)
# 4th of Jun was the last Day of Stay at Home Order	(Pennsylvania)
# Jun 13th was the last Day of Stay at Home Order	(New York) 
# followed by other states so the cases increased again after lifting lockdown.
# Thus, we can split the data into two sets: during lockdown & after lifting lockdown 


us1_ApMy <- as_tibble(us_wave_1) %>%
  mutate(days = seq(1,length(us_wave_1)),
         date = date + days(days),
         month = month(date)) %>%
  filter(month <= 5)

ap <- us1_ApMy %>%
  ggplot(aes(date,value)) +
  geom_line() +
  labs(title = 'US - First Wave - Lockdown', y = 'Number of Cases', x = 'Date')

us1_JuAu <- as_tibble(us_wave_1) %>%
  mutate(days = seq(1,length(us_wave_1)),
         date = date + days(days),
         month = month(date)) %>%
  filter(month >= 6) 

ju <- us1_JuAu %>%
  ggplot(aes(date,value)) +
  geom_line() +
  labs(title = 'US - First Wave - Lockdown Lifted', y = 'Number of Cases', x = 'Date')

ap + ju




################# Exponential Distribution #####################

Exp <- fitdistr(us_wave_1, "exponential")
est_lambda <- Exp$estimate

# Maximum-likelihood of Exponential Distribution 
l_exp <- Exp$loglik

# Exponential Distribution
ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,dexp(us_wave_1,est_lambda)),se = F, color = 'red') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Exponential Distribution")



# Hazard function for Exponential Distribution

ha_exponential <- function(df,den,prob,par1){
  den(df,par1)/(1-prob(df,par1))
}

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_smooth(aes(us_wave_1,ha_exponential(us_wave_1,dexp,pexp,est_lambda)))+
  scale_y_continuous('h(t)')+
  scale_x_continuous('t')




################# Exponentiated Exponential Distribution #####################

fit_ee <- vglm(value ~ 1, fam = expexpff, data = as_tibble(us_wave_1), trace = TRUE, maxit = 99)

parameters <-  as.vector(Coef(fit_ee))
alph_hat <- parameters[2]
lamb_hat <- parameters[1]

# Extracting Log-Likelihood
l_ee <-  logLik(fit_ee)


# Exponentiated exponential distribution
eed <-  lamb_hat*alph_hat * (1 - exp(-1 * lamb_hat * us_wave_1))^(alph_hat - 1) * exp(-1 * lamb_hat * us_wave_1)

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,eed),se = F, color = 'orange') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Exponentiated Exponential Distribution")



# Hazard function for Exponentiated Exponential Distribution
FEE <- function(x,alp,lam){
  (1 - exp(-lam * x) )^alp
}

ha_par_2 <- function(df,den,prob,par1,par2){
  den(df,par1,par2)/(1-prob(df,par1,par2))
}

ha_ee <- function(df,den,prob,par1,par2){
  den/(1-prob)
}

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_smooth(aes(us_wave_1,ha_ee(us_wave_1,eed,FEE(us_wave_1,alph_hat,lamb_hat),lamb_hat,alph_hat)))+
  scale_y_continuous('h(t)')+
  scale_x_continuous('t')




################# Gamma Distribution #####################

fit <- vglm(value ~ 1, fam = gammaR, data = as_tibble(us_wave_1), trace = TRUE)

a <- Coef(fit)[2]
b <- Coef(fit)[1]
l_ga <-  logLik(fit)

# Gamma Distribution
ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,dgamma(us_wave_1,a,b)), color = 'green', se = F) +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Gamma Distribution")



# Hazard function for Gammma Distribution

ha_3 <- function(df,den,prob,par1,par2){
  den(df,par1,par2)/(1-prob(df,par1,par2))
}

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_smooth(aes(us_wave_1,ha_3(us_wave_1,dgamma,pgamma,a,b)))+
  scale_y_continuous('h(t)')+
  scale_x_continuous('t')




################# Weibull Distribution #####################

fit_w <- vglm(value ~ 1, fam = weibullRff, data = as_tibble(us_wave_1), trace = TRUE, maxit = 99)
k <- Coef(fit_w)[2]
l <- Coef(fit_w)[1]
l_w <-  logLik(fit_w)


# Weibull distribution
ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,dweibull(us_wave_1,k,l)),se = F, color = 'blue') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Weibull Distribution")



# Hazard function for Weibull Distribution

ha_weibull <- function(df,den,prob,par1,par2){
  den(df,par1,par2)/(1-prob(df,par1,par2))
}

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_smooth(aes(us_wave_1,ha_weibull(us_wave_1,dweibull,pweibull,k,l)))+
  scale_y_continuous('h(t)')+
  scale_x_continuous('t')




################# Lognormal Distribution #####################

fit_lg <- vglm(value ~ 1, fam = lognormal, data = as_tibble(us_wave_1), trace = TRUE, maxit = 99)
mu <- Coef(fit_lg)[1]
sd <- Coef(fit_lg)[2]
l_lg <-  logLik(fit_lg)


ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,dlnorm(us_wave_1,mu,sd)),se = F, color = 'yellow') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Lognormal Distribution")



# Hazard function for Lognormal Distribution

ha_lognormal <- function(df,den,prob,par1,par2){
  den(df,par1,par2)/(1-prob(df,par1,par2))
}

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_smooth(aes(us_wave_1,ha_lognormal(us_wave_1,dlnorm,plnorm,mu,sd)))+
  scale_y_continuous('h(t)')+
  scale_x_continuous('t')




################# Marshall-Olkin Exponential Distribution #####################

dMOE <-  function(x,l,t,log = FALSE){
  loglik <-  log(l*t) - l*(x) - 2*log(1 - (1-t)*exp(-l*x))
  
  if (log == FALSE)
    density <- exp(loglik)
  
  else density <- loglik
  
  return(density)
}


fit_moe <-  maxlogL(us_wave_1, dist = "dMOE", link = list(over = c("l","t"), fun = c("log_link","log_link")))
lae <- fit_moe$fit$par[1]
te <- fit_moe$fit$par[2]
l_moe <- logLik(fit_moe)


ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,dMOE(us_wave_1,lae,te)),se = F, color = 'cyan') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Marshall Olkin & Exponential Distribution")



# Hazard function for Marshall-Olkin Exponential Distribution

FMOE <- function(x,l,t) {
  (1-exp(-l*x)) / (1- (1-t)*exp(-l*x))
}


ha_moe <- function(df,den,prob,par1,par2){
  den(df,par1,par2)/(1-prob(df,par1,par2))
}

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_smooth(aes(us_wave_1,ha_moe(us_wave_1,dMOE,FMOE,lae,te)))+
  scale_y_continuous('h(t)')+
  scale_x_continuous('t')




################# Marshall-Olkin & Exponentiated Exponential Distribution #####################

dMOEE <-  function(x,a,l,t,log = FALSE){
  loglik <-  log(a*l*t) - l*(x) + (a-1)*log(1-exp(-l*x)) - 2*log(1 - (1-t)*(1 - (1-exp(-l*x))^a))
  
  if (log == FALSE)
    density <- exp(loglik)
  
  else density <- loglik
  
  return(density)
}


fit_moee <-  maxlogL(us_wave_1, dist = "dMOEE", link = list(over = c("a","l","t"), fun = c("log_link","log_link","log_link")))
al <- fit_moee$fit$par[1]
la <- fit_moee$fit$par[2]
t <- fit_moee$fit$par[3]
l_moee <- logLik(fit_moee)


ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,dMOEE(us_wave_1,al,la,t)),se = F, color = 'purple') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Marshall Olkin & EE Distribution")



# Hazard function for Marshall-Olkin Exponentiated Exponential Distribution

FMOEE <- function(x,a,l,t) {
  ((1-exp(-l*x))^a) / (1- (1-t)*(1 - (1 - exp(-l*x))^a))
}


ha_moee <- function(df,den,prob,par1,par2,par3){
  den(df,par1,par2,par3)/(1-prob(df,par1,par2,par3))
}

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_smooth(aes(us_wave_1,ha_moee(us_wave_1,dMOEE,FMOEE,al,la,t)))+
  scale_y_continuous('h(t)')+
  scale_x_continuous('t')




################# Marshall-Olkin Weibull Distribution #####################

dMOW <-  function(x,a,l,t,log = FALSE){
  loglik <-  log(a*l*t) - l*x^a + (a-1)*log(x) - 2*log(1 - (1-t)*exp(-l*x^a))
  
  if (log == FALSE)
    density <- exp(loglik)
  
  else density <- loglik
  
  return(density)
}


fit_mow <-  maxlogL(us_wave_1, dist = "dMOW", link = list(over = c("a","l","t"), fun = c("log_link","log_link","log_link")))
alw <- fit_mow$fit$par[1]
law <- fit_mow$fit$par[2]
tw <- fit_mow$fit$par[3]
l_mow <- logLik(fit_mow)


ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,dMOW(us_wave_1,alw,law,tw)),se = F, color = 'brown') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Marshall Olkin & Weibull Distribution")



# Hazard function for Marshall-Olkin Weibull Distribution

FMOW <- function(x,a,l,t) {
  (1-exp(-l*x^a)) / (1- (1-t)*exp(-l*x^a))
}


ha_mow <- function(df,den,prob,par1,par2,par3){
  den(df,par1,par2,par3)/(1-prob(df,par1,par2,par3))
}

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_smooth(aes(us_wave_1,ha_mow(us_wave_1,dMOW,FMOW,alw,law,tw)))+
  scale_y_continuous('h(t)')+
  scale_x_continuous('t')




################# Marshall-Olkin Poisson - EE Distribution #####################

dMOPEE <-  function(x,a,b,l,log = FALSE){
  loglik <-  log(a*b) - b*x + (a-1)*log(1-exp(-b*x)) + a*log(1 + l*(1- (1+exp(-b*x))^a)) - l*(1- exp(-b*x))^a
  
  if (log == FALSE)
    density <- exp(loglik)
  
  else density <- loglik
  
  return(density)
}


fit_mopee <-  maxlogL(us_wave_1, dist = "dMOPEE", link = list(over = c("a","b","l"), fun = c("log_link","log_link","log_link")))
alp <- fit_mopee$fit$par[1]
beta <- fit_mopee$fit$par[2]
lap <- fit_mopee$fit$par[3]
l_mopee <- logLik(fit_mopee)


ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,dMOPEE(us_wave_1,alp,beta,lap)),se = F, color = 'purple') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Marshall Olkin Poisson - EE Distribution")



# Hazard function for Marshall-Olkin Poisson - EE Distribution

FMOPEE <- function(x,a,b,l) {
  1- ((1- (1-exp(-b*x))^a) / exp(l*(1 - exp(-b*x))^a))
}


ha_mopee <- function(df,den,prob,par1,par2,par3){
  den(df,par1,par2,par3)/(1-prob(df,par1,par2,par3))
}

ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_smooth(aes(us_wave_1,ha_mow(us_wave_1,dMOPEE,FMOPEE,alp,beta,lap)))+
  scale_y_continuous('h(t)')+
  scale_x_continuous('t')




################# Summary #####################

# Combined plot
ggplot(as_tibble(us_wave_1), aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(us_wave_1,dexp(us_wave_1,est_lambda), color = 'Exponential'),se = F) +
  geom_smooth(aes(us_wave_1,eed, color = 'Exponentiated Exponential'),se = F) +
  geom_smooth(aes(us_wave_1,dgamma(us_wave_1,a,b), color = 'Gamma'), se = F) +
  geom_smooth(aes(us_wave_1,dweibull(us_wave_1,k,l), color = "Weibull"),se = F) +
  geom_smooth(aes(us_wave_1,dlnorm(us_wave_1,mu,sd), color = 'Lognormal'),se = F) +
  geom_smooth(aes(us_wave_1,dMOE(us_wave_1,lae,te), color = 'MOE'),se = F) +
  geom_smooth(aes(us_wave_1,dMOEE(us_wave_1,al,la,t), color = 'MOEE'),se = F) +
  geom_smooth(aes(us_wave_1,dMOW(us_wave_1,alw,law,tw), color = 'MOW'),se = F) +
  geom_smooth(aes(us_wave_1,dMOPEE(us_wave_1,alp,beta,lap), color = 'MOPEE'),se = F) +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Combined Plot") +
  scale_color_discrete(name='Distribution:') +
  theme(legend.position = "bottom")



# Testing the function: Kolmogorov-Smirnov Tests
p1 <- ks.test(us_wave_1,"pexp",est_lambda)

FEE <- function(x,alp,lam){
  (1 - exp(-lam * x) )^alp
}

p2 <- ks.test(us_wave_1,"FEE",alph_hat,lamb_hat)
p3 <- ks.test(us_wave_1,"pgamma",a,b)
p4 <- ks.test(us_wave_1,"pweibull",k,l)
p5 <- ks.test(us_wave_1,"plnorm",mu,sd)

FMOE <- function(x,l,t) {
  (1-exp(-l*x)) / (1- (1-t)*exp(-l*x))
}

p6 <-  ks.test(us_wave_1,"FMOE",lae,te)


FMOEE <- function(x,a,l,t) {
  ((1-exp(-l*x))^a) / (1- (1-t)*(1 - (1 - exp(-l*x))^a))
}

p7 <-  ks.test(us_wave_1,"FMOEE",al,la,t)


FMOW <- function(x,a,l,t) {
  (1-exp(-l*x^a)) / (1- (1-t)*exp(-l*x^a))
}

p8 <-  ks.test(us_wave_1,"FMOW",alw,law,tw)


FMOPEE <- function(x,a,b,l) {
  1 - ((1- (1-exp(-b*x))^a) / exp(l*(1 - exp(-b*x))^a))
}


p9 <-  ks.test(us_wave_1,"FMOPEE",alp,beta,lap)



# Table of Comparison
Dist <-  c("Exponential","Exponentiated Exponential", "Gamma", "Weibull", "Lognormal","MOE","MOEE","MOW","MOPEE")
par1 <-  c(est_lambda,alph_hat,a, k,mu,lae,al,alw,alp)
par2 <-  c(est_lambda[2],lamb_hat,b,l,sd,te,la,law,beta)
par3 <-  c(est_lambda[3],lamb_hat[2],b[2],l[2],sd[2],te[2],t,tw,lap)
Likelihood <-  c(l_exp,l_ee, l_ga, l_w,l_lg,l_moe,l_moee,l_mow,l_mopee)
p_value <-  c(p1$p.value, p2$p.value, p3$p.value, p4$p.value,p5$p.value,p6$p.value,p7$p.value,p8$p.value,p9$p.value)
cbind(Dist,par1,par2,par3,Likelihood,p_value)