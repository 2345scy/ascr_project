library(devtools)
load_all("D:/work/ben/ascr")
library(ascr)
library(mgcv)
source('functions in fit.r')
source('create.capt.r')
#source('function in progress.r')
source('support functions.r')
set.seed(810)
captures.m <- data.frame(session = c(1, 1, 1, 1, 3, 3, 3, 3, 3),
                         id = c(1, 2, 1, 2, 1, 2, 3, 4, 4),
                         occasion = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
                         trap = c(2, 1, 3, 2, 4, 3, 2, 6, 5),
                         bearing = c(5.59814552863741, 1.55785557290341, 0.226989015597898,
                                     6.0614780856239, 4.45397150613391, 4.11997712380155,
                                     5.63614140398099, 4.54698079620703, 3.70723587716335),
                         toa = runif(9, 0, 3),
                         ss = rnorm(9, 5, 1))

traps.m <- list(data.frame(x = 1:3, y = rep(0, 3)),
                data.frame(x = 1:3, y = rep(0, 3)),
                data.frame(x = rep(1:4, 2), y = rep(c(0, 1), each = 4)),
                data.frame(x = 1:3, y = rep(0, 3)))

#mrds.locs.m = vector('list', 4)
#mrds.locs.m[[1]] = matrix(rnorm(4), ncol = 2)
#mrds.locs.m[[3]] = matrix(rnorm(8), ncol = 2)
#mrds.locs.m.name = lapply(mrds.locs.m, function(x) {if(!is.null(x)) {
#  x = as.data.frame(x)
#  x$ID = 1:nrow(x)
#  return(x)
#}})
#capt.m.mrds <- create.capt(captures.m, traps = traps.m, mrds.locs = mrds.locs.m)

capt.m = create.capt(captures.m, traps = traps.m)
mask.m <- create.mask(traps.m, 10)


#for the example of 1 session only
#captures.1 <- subset(captures.m, session == 1)
#traps.1 <- traps.m[[1]] 
#mrds.locs.1 = mrds.locs.m[[1]]
#mrds.locs.1.name = mrds.locs.m.name[[1]]
#capt.1 = create.capt(captures.1, traps = traps.1)
#capt.1.mrds = create.capt(captures.1, traps = traps.1, mrds.locs = mrds.locs.1)
#mask.1 = mask.m[[1]]



################################################################
#multi-session example test

o = capture.fun(capt = capt.m)

dims = o$dims

data.full = o$data.capt

bucket_info = o$bucket_info


###########################################
data.traps = trap.fun(traps = traps.m, dims = dims)

data.full = merge(data.full, data.traps, by = c('session', 'trap'), all = TRUE)

###################################################################################
sound.speed = 331
o.mask = mask.fun(mask = mask.m, dims = dims, 
                  data.traps = data.traps, data.full = data.full,
                  local = TRUE, bucket_info = bucket_info, sound.speed = sound.speed)

dims[["n.masks"]] = o.mask$n.masks
data.full = o.mask$data.full
bucket_info = o.mask$bucket_info
#keep the old version of output temporarily, might be useful later
mask = o.mask$mask
#dists = o.mask$dists
A = o.mask$A
buffer = o.mask$buffer
#all.which.local = o.mask$all.which.local
#all.n.local = o.mask$all.n.local
data.dists.thetas = o.mask$data.dists.thetas
data.ID_mask = o.mask$data.ID_mask
data.mask = o.mask$data.mask

#####################################################################################

#be careful 'D' should be the last one, because the default start value of D relies on others
fulllist.par = c('g0', 'sigma', 'lambda0', 'z', 'shape.1', 
                 'shape.2', 'shape', 'scale', 'b0.ss', 'b1.ss',
                 'b2.ss', 'sigma.ss', 'kappa', 'alpha', 'sigma.toa',
                 "sigma.b0.ss", 'D')

#requirement of mask level data input is, it must have the same structure as "mask"
#for each session, it must have two columns 'x' and 'y' which are exactly the same with "mask"
#and any columns include the information as explanatory variables
#if the user only need 'x' and 'y' as explanatory variables, this 'mask' level data input is not required
ihd.covariates = vector('list', dims$n.sessions)
for(i in 1:dims$n.sessions){
  ihd.covariates[[i]] = data.frame(var1 = rnorm(dims$n.masks[i]))
  ihd.covariates[[i]] = cbind(mask[[i]], ihd.covariates[[i]])
}

#"model" here do not accept interaction between variables from 'mask' level and other levels
#do not accept offset term in formula

par.extend = list(data = list(session = data.frame(session = 1:4,
                                        weather = c("rain", "rain", "cloudy", "cloudy")),
                              trap = data.frame(session = rep(1:4, dims$n.traps),
                                        trap = c(1:3,1:3,1:8,1:3), 
                                        brand = rep(c("a","b"), c(9,8))),
                              ID = data.frame(session = c(1,1,3,3,3,3), ID = c(1,2,1,2,3,4), 
                                        loc = c('high land', 'high land', 'forest', 'forest', 
                                        'lake', 'sea')),
                              mask = ihd.covariates), 
                  model = list(sigma = ~weather*brand, g0 = ~brand + loc + var1, 
                               D = ~te(x, y) + var1), 
                  link = list(sigma = 'log'),
                  scale = list(g0 = FALSE))

#par.extend = NULL

o.par.extend = par.extend.fun(par.extend = par.extend, data.full = data.full, data.mask = data.mask, 
                              dims = dims, fulllist.par = fulllist.par)

data.full = o.par.extend$data.full
data.mask = o.par.extend$data.mask
data.par = o.par.extend$data.par
name.extend.par = o.par.extend$name.extend.par

#####################################################################################

sv = list(g0 = 0.5, b1.ss = 1.5)
fix = list(sigma = 1, sigma.b0.ss = 1)

ss.opts = list(cutoff = 4.2)
o.ss = ss.fun(ss.opts = ss.opts, data.full = data.full, dims = dims, 
              bucket_info = bucket_info, sv = sv, fix = fix)

data.full = o.ss$data.full
dims = o.ss$dims
bucket_info = o.ss$bucket_info
ss.opts = o.ss$ss.opts

######################################################################

survey.length = NULL
cue.rates = NULL
o.CR_SL = CR_SL(cue.rates = cue.rates, survey.length = survey.length,
                bucket_info = bucket_info, dims = dims)

survey.length = o.CR_SL$survey.length
bucket_info = o.CR_SL$bucket_info
mu.rates = o.CR_SL$mu.rates


detfn = NULL
bounds = NULL

o.param = param.detfn.fun(sv = sv, fix = fix, bounds = bounds, name.extend.par = name.extend.par,
                          detfn = detfn, data.full = data.full, data.mask = data.mask, 
                          data.par = data.par, ss.opts = ss.opts, bucket_info = bucket_info,
                          fulllist.par = fulllist.par, A = A, buffer = buffer,
                          survey.length = survey.length, dims = dims)

DX.full = o.param$design.matrices.full
DX.mask = o.param$design.matrices.mask
data.par = o.param$data.par
sv.input = o.param$sv.input
fix.input = o.param$fix.input
bounds.input = o.param$bounds.input
name.extend.par = o.param$name.extend.par
detfn = o.param$detfn
param.og = o.param$param.og


#####################################################################
#####################################################################

library(TMB)
library(MASS)

data.full = numeric_ID(data.full)
data.ID_mask = numeric_ID(data.ID_mask)


data.full = sort.data(data.full, "data.full")
data.mask = sort.data(data.mask, "data.mask")
data.ID_mask = sort.data(data.ID_mask, "data.ID_mask")
data.dists.thetas = sort.data(data.dists.thetas, "data.dists.thetas")


data <- list(n_sessions = dims$n.sessions,
             n_animals = dims$n.animals,
             n_IDs = dims$n.IDs,
             n_traps = dims$n.traps,
             n_masks = dims$n.masks,
             
             #code of detfn_index:
             #1:hn, 2:hhn, 3:hr, 4:th, 5:lth, 6:ss, 7:ss_dir, 8:ss_het
             detfn_index = numeric_detfn(detfn),
             
             #simply use the row index of 'data.par' as the numeric subtitution of parameters
             #1:g0, 2:sigma; 3:lambda0, 4:z, 5:shape.1, 6:shape.2, 7:shpe, 8:scale, 9:b0.ss,
             #10:b1.ss, 11:b2.ss, 12:sigma.ss, 13:kappa, 14:alpha, 15:sigma.toa,
             #16:sigma.b0.ss, 17:D
             param.og = (1:nrow(data.par))[which(data.par[['par']] %in% param.og)],
             par_n_col = as.matrix(data.par[, c('n_col_full', 'n_col_mask')]),
             par_link = numeric_link(data.par[['link']]),
             
             is_animalID = as.numeric("animal_ID" %in% colnames(data.full)),
             is_ss_origin = as.numeric(all("ss" %in% bucket_info,
                                          !"ss.het" %in% bucket_info,
                                          !"ss.dir" %in% bucket_info)),
             is_ss_het = as.numeric("ss.het" %in% bucket_info),
             is_ss_dir = as.numeric("ss.dir" %in% bucket_info),
             is_bearing = as.numeric("bearing" %in% bucket_info),
             is_toa = as.numeric("toa" %in% bucket_info),
             is_dist = as.numeric("dist" %in% bucket_info),
             is_local = as.numeric("local" %in% bucket_info),
             is_freqs = as.numeric("freqs" %in% bucket_info),
             
             capt_bin = ifelse(is.na(data.full$bincapt), 0, data.full$bincapt),
             capt_bearing = ifelse(is.na(data.full$bearing), 0, data.full$bearing),
             capt_dist = ifelse(is.na(data.full$dist), 0, data.full$dist),
             capt_ss = ifelse(is.na(data.full$ss), 0, data.full$ss),
             capt_toa = ifelse(is.na(data.full$toa), 0, data.full$toa),
             
             dx = data.dists.thetas$dx,
             theta = data.dists.thetas$theta,
             
             index_local = data.ID_mask$local,
             toa_ssq = data.ID_mask$toa_ssq
             )

for(i in fulllist.par){
  data[[paste0(i, "_x")]] = DX.full[[i]]
  data[[paste0(i, "_x_mask")]] = DX.mask[[i]]
}


parameters <- list(b0=0, b1 = 0, b2 = 0, b3 = 0, logv = 0, logsigma2_a = 0, logitpho = 0, u = numeric(length(Y)))
