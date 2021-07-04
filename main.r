library(devtools)
load_all("G:/work/ben/ascr")
library(ascr)
library(mgcv)
library(stringr)
library(TMB)
library(MASS)

#compile("fit_ascr_version7.cpp")

setwd('C:/Users/haoli/Documents/GitHub/ascr_project')

source('create.capt.r')
#source('function in progress.r')

source('test_data_preparation.r')
source('fit_ascr_tmb.r')
#available test data:
#simple fitting -- half normal
#simple fitting -- hazard halfnormal
#simple fitting -- hazard rate
#bearing fitting
#dist fitting
#toa fitting
#joint bearing/dist fitting
#ss fitting
#joint ss/toa fitting
#Inhomogeneous density estimation
#Multi-session models
test_data("Multi-session models")

#the par.extend below is for "Multi-session models"
#par.extend = list(data = list(session = data.frame(session = 1:2, weather = c('rain', 'sunny')),
#                              trap = data.frame(session = rep(1:2, c(6,3)), trap = c(1:6, 1:3), 
#                                                brand = sample(c('a','b','c'), size = 9, replace = T))),
#                  model = list(g0 = "~weather + brand"))

capt_input = create.capt(captures, traps = traps)


out = fit.ascr.tmb(capt = capt_input, traps = traps, mask = mask, sv = sv_input, fix = fix_input,
                   local = local_input, bounds = bounds_input, cue.rates = cue.rates, 
                   sound.speed = sound.speed, ss.opts = ss.opts, survey.length = survey.length,
                   par.extend = par.extend)




#enter into fit.ascr.tmb()

capt = capt_input
sv = sv_input
fix = fix_input
local = local_input
bounds = bounds_input


