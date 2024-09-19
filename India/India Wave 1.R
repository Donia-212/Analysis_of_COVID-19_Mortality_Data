rm(list=ls())
library(MASS)
library(stats4)
library(fitdistrplus)
library(tidyverse)
library(lubridate)
library(patchwork)
library(VGAM)
library(VGAMextra)


################# India Wave 1 Dataset #####################

india_wave_1 <-  unlist(read.csv('Ind_Wave_1.csv'))

# Visualizing the data:
ggplot(as_tibble(india_wave_1), aes(value)) +
  geom_histogram(bins = 30, color = 'black', fill = 'grey') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Frequency') +
  labs(title = "India - First Wave")

# 5 data summary 
summary(india_wave_1)
var(india_wave_1)




################# Investigating the data #####################

# The First wave: from 4/2020 to 10/2020  
date <- dmy('1/4/2020')

as_tibble(india_wave_1) %>%
  mutate(days = seq(1,length(india_wave_1)), date = date + days(days)) %>%
  ggplot(aes(date,value)) +
  geom_line() +
  labs(title = "India - First Wave", y = 'Number of Cases', x = 'Date')

as_tibble(india_wave_1) %>%
  mutate(days = seq(1,length(india_wave_1)), date = date + days(days)) %>%
  print(n = Inf)%>%
  arrange(desc(value)) %>%
  View()


# We have an outlier towards mid November when the death cases hit 2003. 
# This could be due to public holidays that occurred towards the end of October.

# List of holidays:
# 1- Oct	22nd: Maha Saptami
# 2- Oct 23rd: Maha Ashtami	
# 3- Oct 24th: Maha Navami
# 4- Oct 25th: Dussehra	
# 5- Oct 30th: Milad un-Nabi (15.5% of India's population are Muslims)
# 6- Oct 31st: Maharishi Valmiki Jayanti	Restricted Holiday

# It is possible that people began to leave their houses to visit each other thus causing the number
# of death cases to increase rapidly.

# Another outlier would be in the beginning of January where the death cases suddenly became zero.
# There are few reasons for such a thing to occur:
# 1- There are no deaths (which is very unlikely as the it is obvious from the graph that the number of deaths did not stop)
# 2- There was no data recorded on that day for the number of deaths.
# 3- The system was down


# Removing the values when they number of death cases are not zero (missing values)
in1 <- as_tibble(india_wave_1) %>%
  filter(value != 0, value != 2003, value!= 1129) 

in1 %>%
  mutate(days = seq(1,length(unlist(in1))), date = date + days(days)) %>%
  ggplot(aes(date,value)) +
  geom_line()+
  labs(title = "India - First Wave - Modified", y = 'Number of Cases', x = 'Date')

summary(in1)
var(in1)



################# Exponential Distribution #####################

# Maximum-likelihood fitting of univariate distribution
Exp <- fitdistr(unlist(in1), "exponential")
est_lambda <- Exp$estimate
l_exp <- Exp$loglik


# Exponential Distribution
ggplot(in1, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in1),dexp(unlist(in1),est_lambda)),se = F, color = 'red') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Exponential Distribution")




################# Exponentiated Exponential Distribution #####################

fit_ee <- vglm(value ~ 1, fam = expexpff, data = in1, trace = TRUE, maxit = 99)
parameters <-  as.vector(Coef(fit_ee))
alph_hat <- parameters[2]
lamb_hat <- parameters[1]

# Extracting Log-Likelihood
l_ee <-  logLik(fit_ee)


# Exponentiated exponential distribution
eed <-  lamb_hat*alph_hat * (1 - exp(-1 * lamb_hat * unlist(in1)))^(alph_hat - 1) * exp(-1 * lamb_hat * unlist(in1))
ggplot(in1, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in1),eed),se = F, color = 'orange') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Exponentiated Exponential Distribution")




################# Gamma Distribution #####################

fit <- vglm(value ~ 1, fam = gammaR, data = in1, trace = TRUE)

a <- Coef(fit)[2]
b <- Coef(fit)[1]
l_ga <-  logLik(fit)

# Gamma Distribution
ggplot(in1, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in1),dgamma(unlist(in1),a,b)), color = 'green', se = F) +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Gamma Distribution")




################# Weibull Distribution #####################

fit_w <- vglm(value ~ 1, fam = weibullRff, data = in1, trace = TRUE, maxit = 99)
k <- Coef(fit_w)[2]
l <- Coef(fit_w)[1]
l_w <-  logLik(fit_w)


# Weibull distribution
ggplot(in1, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in1),dweibull(unlist(in1),k,l)),se = F, color = 'blue') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Weibull Distribution")





################# Lognormal Distribution #####################

fit_lg <- vglm(value ~ 1, fam = lognormal, data = in1, trace = TRUE, maxit = 99)
mu <- Coef(fit_lg)[1]
sd <- Coef(fit_lg)[2]
l_lg <-  logLik(fit_lg)


ggplot(in1, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in1),dlnorm(unlist(in1),mu,sd)),se = F, color = 'yellow') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Lognormal Distribution")




################# Generalized Poisson Distribution #####################

fit_genp <- vglm(value ~ 1, fam = genpoisson0, data = in1, trace = TRUE, maxit = 99)
th <- Coef(fit_genp)[1]
la <- Coef(fit_genp)[2]
l_genp <-  logLik(fit_genp)

ggplot(in1, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in1),dgenpois0(unlist(in1),th,la)),se = F, color = 'maroon') +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Generalized Poisson Distribution")





################# Summary #####################

# Combined plot: Modified
ggplot(in1, aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'grey') +
  geom_smooth(aes(unlist(in1),dexp(unlist(in1),est_lambda), color = 'Exponential'),se = F) +
  geom_smooth(aes(unlist(in1),eed, color = 'Exponentiated Exponential'),se = F) +
  geom_smooth(aes(unlist(in1),dgamma(unlist(in1),a,b), color = 'Gamma'), se = F) +
  geom_smooth(aes(unlist(in1),dweibull(unlist(in1),k,l), color = "Weibull"),se = F) +
  geom_smooth(aes(unlist(in1),dlnorm(unlist(in1),mu,sd), color = 'Lognormal'),se = F) +
  geom_smooth(aes(unlist(in1),dgenpois0(unlist(in1),th,la), color = 'Generalized Poisson'),se = F) +
  scale_x_continuous('Cases') +
  scale_y_continuous('Density') +
  labs(title = "Combined Plot: Modified") +
  scale_color_discrete(name='Distribution:') +
  theme(legend.position = "bottom")


# Testing the function: Kolmogorov-Smirnov Tests
p1 <- ks.test(in1,"pexp",est_lambda)

FEE <- function(x,alp,lam){
  (1 - exp(-lam * x) )^alp
}

p2 <- ks.test(in1,"FEE",alph_hat,lamb_hat)
p3 <- ks.test(in1,"pgamma",a,b)
p4 <- ks.test(in1,"pweibull",k,l)
p5 <- ks.test(in1,"plnorm",mu,sd)
p6 <- ks.test(in1,"pgenpois0",th,la)


# Table of Comparison
Dist <-  c("Exponential","Exponentiated Exponential", "Gamma", "Weibull", "Lognormal", 'Generalized Poisson')
par1 <-  c(est_lambda,alph_hat,a, k,mu,th)
par2 <-  c(est_lambda[2],est_lambda,b,l,sd,la)
Likelihood <-  c(l_exp,l_ee, l_ga, l_w,l_lg,l_genp)
p_value <-  c(p1$p.value, p2$p.value, p3$p.value, p4$p.value,p5$p.value,p6$p.value)
cbind(Dist,par1,par2,Likelihood,p_value)
