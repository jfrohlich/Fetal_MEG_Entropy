rm(list = ls()) # clear all existing vars
set.seed(111) # seed RNG
#install.packages("lme4")
#install.packages("mediation")
library(lme4)
library(MASS)
library("mediation")

# Read table generated from all patients, including burst suppression
T0 = read.csv("C:/Users/joelf/OneDrive/Documents/Research/Preissl/FetalData.csv",'header'=TRUE)
# Mean center data
T0$CTW  = scale(T0$CTW, scale=FALSE)
T0$LZC  = scale(T0$LZC, scale=FALSE)
T0$PE32 = scale(T0$PE32, scale=FALSE)
T0$PE64 = scale(T0$PE64, scale=FALSE)
T0$MSE  = scale(T0$MSE, scale=FALSE)
T0$SE   = scale(T0$SE, scale=FALSE)
T0$mval = scale(T0$mval, scale=FALSE)

nsm = 1000 # number of simulations to run
useme = T0 # data table to use

# LZC 
XpredY <- lmer(LZC ~ XChrom + (1|ID),data=useme) 
XMpredY <- lmer(LZC ~ XChrom + mval + (1|ID),data=useme)

med <- mediate(XpredY,XMpredY,treat="XChrom",mediator="mval",boot=F,sims=nsm) # run mediation model
summary(med) 

# CTW 
XpredY <- lmer(CTW ~ XChrom + (1|ID),data=useme) 
XMpredY <- lmer(CTW ~ XChrom + mval + (1|ID),data=useme)

med <- mediate(XpredY,XMpredY,treat="XChrom",mediator="mval",boot=F,sims=nsm) # run mediation model
summary(med) 


# PE32 
XpredY <- lmer(PE32 ~ XChrom + (1|ID),data=useme) 
XMpredY <- lmer(PE32 ~ XChrom + mval + (1|ID),data=useme)

med <- mediate(XpredY,XMpredY,treat="XChrom",mediator="mval",boot=F,sims=nsm) # run mediation model
summary(med) 


# PE64 
XpredY <- lmer(PE64 ~ XChrom + (1|ID),data=useme) 
XMpredY <- lmer(PE64 ~ XChrom + mval + (1|ID),data=useme)

med <- mediate(XpredY,XMpredY,treat="XChrom",mediator="mval",boot=F,sims=nsm) # run mediation model
summary(med) 


# MSE 
XpredY <- lmer(MSE ~ XChrom + (1|ID),data=useme) 
XMpredY <- lmer(MSE ~ XChrom + mval + (1|ID),data=useme)

med <- mediate(XpredY,XMpredY,treat="XChrom",mediator="mval",boot=F,sims=nsm) # run mediation model
summary(med) 


# SE 
XpredY <- lmer(SE ~ XChrom + (1|ID),data=useme) 
XMpredY <- lmer(SE ~ XChrom + mval + (1|ID),data=useme)

med <- mediate(XpredY,XMpredY,treat="XChrom",mediator="mval",boot=F,sims=nsm) # run mediation model
summary(med) 


