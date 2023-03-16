############################################################
#Aim: to visualize the fitted models of Ward et al. (2022) BMJ
#Final edit: 16 March 2022
#Editor: Fumi Miura
############################################################
###Procedure
#0. Package 
#1. Data 
#2. Model parameters
#3. Visualization
############################################################

###0. Package -----
library(tidyverse)

###1. Data -----
## raw data of Fig-5
fig5_data <- data.frame( #Retrived from Fig-5 in Ward et al. 
  SI=c(3,3,1,4,5,2,10,6,7,6,8,11,10)
)
## posterior samples of the best fit model (gamma)
output_ukhsa <- read_csv("output_trace.csv") #source: authors' Github https://github.com/OvertonC/Transmission-Dynamics-of-Monkeypox-in-the-United-Kingdom

###2. Model parameters -----
##Gamma; shape=0.804, rate=0.101 (source: posterior means provided on Github)
gamma_output <- output_ukhsa %>%
  summarise(mean(alpha),median(alpha),mean(beta),median(beta)) #0.810 0.804 0.103 0.101 #all chains
shape_gam <- gamma_output$`median(alpha)` #shape = 0.804
rate_gam <- gamma_output$`median(beta)`# rate = 0.101

##Weibull; shape=0.912, scale=7.850 (source: Suup B table, Weibull (ICC))
sig <- 9.0
mu <- 8.2
f <- function(k){sig^2/mu^2 - gamma(1+2/k)/gamma(1+1/k)^2 + 1} #function to translate (sigma,mu) into shape para
shape_wei <- uniroot(f, c(0.1, 10))$root 
scale_wei <- mu/gamma(1+1/shape_wei) 

##Lognormal; para_sd=1.205, para_mu=1.401 (source: Suup B table, lognormal (ICC))
mean_ln <- 8.4
sd_ln <- 15.2
para_sd <- (log((sd_ln^2)/(mean_ln^2) + 1))^(1/2)
para_mu <- log(mean_ln) - (1/2)*para_sd^2

###3. Visualization -----
##data.frame for plot
x_plot <- seq(0,30,by=0.1)
no_model <- 3  #no. of tested models 
model_plot <- data.frame(
  x_plot = rep(x_plot, no_model),
  dens = c(dgamma(x_plot, shape = shape_gam, rate = rate_gam),
           dweibull(x_plot,shape= shape_wei,scale=scale_wei),
           dlnorm(x_plot,meanlog=para_mu,sdlog=para_sd)
  ),
  model = c(rep("Gamma",length(x_plot)),
            rep("Weibull",length(x_plot)),
            rep("LogNormal",length(x_plot)))
)
##ggplot
ggplot() + 
  geom_histogram(data=fig5_data, 
                 mapping = aes(x = SI, y = ..density..),
                 binwidth = 1,
                 colour = 1, fill = "grey") +
  geom_line(data = model_plot,
            mapping = aes(x=x_plot, y=dens, color=model),
            size=1.2) +
  theme_bw()+
  scale_colour_manual(values=c(RColorBrewer::brewer.pal(5, "RdBu")[1],RColorBrewer::brewer.pal(5, "RdBu")[2],RColorBrewer::brewer.pal(5, "RdBu")[5]))+
  labs(x = "Serial interval (days)", y = "Density")
