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
test_data("Inhomogeneous density estimation")
#short_name = "mul_s_with_extension"

#the par.extend below is for "Multi-session models"
#par.extend = list(data = list(session = data.frame(session = 1:2, weather = c('rain', 'sunny')),
#                             trap = data.frame(session = rep(1:2, c(6,3)), trap = c(1:6, 1:3), 
#                                                brand = sample(c('a','b','c'), size = 9, replace = T))),
#                  model = list(g0 = "~weather", sigma = "~ brand"))


#the par.extend below is for "Inhomogeneous density estimation"
#par.extend$data[['trap']] = data.frame(session = 1, trap = 1:6, brand = c(rep('a', 3), rep('b', 3)))
#par.extend$model[['sigma']] = "~brand"


capt_input = create.capt(captures, traps = traps)

source('fit_ascr_tmb.r')
fit = fit.ascr.tmb(capt = capt_input, traps = traps, mask = mask, sv = sv_input, fix = fix_input,
                   detfn = detfn, local = local_input, bounds = bounds_input, cue.rates = cue.rates, 
                   sound.speed = sound.speed, ss.opts = ss.opts, survey.length = survey.length,
                   par.extend = par.extend)




################################################################################################
#simulation
source('sim_fit_tmb.r')

#sim_capture = simcapt_from_fit_tmb(fit)


n.sim = 100

sim_coef_result = vector('list', n.sim)
sim_sd_result = vector('list', n.sim)

for(i in 1:n.sim){
  set.seed(i)
  sim_capture = simcapt_from_fit_tmb(fit)
  sim_capt_input = create.capt(sim_capture, traps = traps)
  #replace capt_input by simulated capture history, keep others to be the same
  args.input = fit$args
  args.input$capt = sim_capt_input
  fit_sim = do.call('fit.ascr.tmb', args.input)
  sim_coef_result[[i]] = fit_sim$output.tmb$coef_link
  sim_sd_result[[i]] = fit_sim$output.tmb$se_link
  print(paste0(i, "/", n.sim, " finished."))
}

sim_coef_agg = sim_result_agg(sim_coef_result)
#sim_sd_agg = sim_result_agg(sim_sd_result)

pre_fix = paste0("./simulation_result/improved_simulation_ver2/", short_name)

dir.create(pre_fix)

sim_result_plot(sim_coef_agg, fit$output.tmb$coef_link, pre_fix)


#apparently the estimation of std in simulation does not look as well as coef,
#however, it is possible that the std estimation from the model is not reliable
#in the first place
#sim_result_plot(sim_sd_agg, fit$output.tmb$se_link)


#########################################################################################

#plot detection function
source('detfn_tmb.r')
#sample of multi-session with g0 and sigma extension

#plot only intercept
show_detfn_tmb(fit)

#plot with some extension covariates
new_data = data.frame(weather = 'rain', brand = 'c')
show_detfn_tmb(fit, new_covariates = new_data)

#plot with multiple sets of extension covariates 
new_data = data.frame(weather = c('rain', 'rain', 'rain'), brand = c('b', 'c', 'a'))
show_detfn_tmb(fit, new_covariates = new_data)

#plot with extension of sigma, but ignore the extension on g0
new_data = data.frame(weather = c('rain', 'sunny', 'rain'), brand = c('b', 'c', 'a'))
show_detfn_tmb(fit, new_covariates = new_data, param_extend_skip = 'g0')


#sample of ihd
new_data = data.frame(x1 = 0.002, y1 = 0.05, brand = c('b', 'a'))
show_detfn_tmb(fit, new_covariates = new_data, col = c(2, 4))


##########################################################################################
#demo for method vcov and coef
#for ihd
source('methods.r')
vcov(fit)

new_data1 = data.frame(x1 = 0.002, y1 = 0.05, brand = 'b')
new_data2 = data.frame(brand = 'b')

vcov(fit, type = 'fitted')
vcov(fit, type = 'fitted', par = 'sigma')
vcov(fit, type = 'fitted', new_covariates = new_data1)
vcov(fit, type = 'fitted', par = 'sigma', new_covariates = new_data2)

coef(fit, types = 'fitted')
coef(fit, types = 'fitted', par = 'sigma')
coef(fit, types = 'fitted', new_covariates = new_data1)
coef(fit, types = 'fitted', par = 'sigma', new_covariates = new_data2)

stdEr(fit, types = 'fitted')
stdEr(fit, types = 'fitted', par = 'sigma')
stdEr(fit, types = 'fitted', new_covariates = new_data1)
stdEr(fit, types = 'fitted', par = 'sigma', new_covariates = new_data2)
