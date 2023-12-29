cat("\014")
rm(list = ls())

source('TD_Delta_Theta_function.R')

start <- c(0,0) #the intersection point
T <- 6 #the number of generations
N <- 2 #the number of children of each parent node
Tij <- 1 #the number of observations of each node
alpha <- 0.95 #the significance level

rho <- 0.5
sigma <- c(0.1,0.1)
mu <- c(1,1)

deltaanglefunc(rho,mu,sigma,T,N,Tij,alpha,start)