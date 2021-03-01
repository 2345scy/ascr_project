library(devtools)
load_all("G:/work/ben/ascr")
library(ascr)
library(mgcv)
setwd('C:/Users/haoli/Documents/GitHub/ascr_project')
source('functions_in_fit_ascr_new.r')
source('create.capt.r')
#source('function in progress.r')
source('support_functions.r')
source('test_data_preparation.r')


test_data("simple fitting -- half normal")



capt_input = create.capt(captures, traps = traps)

o = capture.fun(capt = capt_input)
dims = o$dims
data.full = o$data.capt
bucket_info = o$bucket_info

data.traps = trap.fun(traps = traps, dims = dims)
data.full = merge(data.full, data.traps, by = c('session', 'trap'), all = TRUE)


o.mask = mask.fun(mask = masks, dims = dims, 
                  data.traps = data.traps, data.full = data.full,
                  local = local_input, bucket_info = bucket_info, sound.speed = sound.speed)

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

fulllist.par = c('g0', 'sigma', 'lambda0', 'z', 'shape.1', 
                 'shape.2', 'shape', 'scale', 'b0.ss', 'b1.ss',
                 'b2.ss', 'sigma.ss', 'kappa', 'alpha', 'sigma.toa',
                 "sigma.b0.ss", 'D')
o.par.extend = par.extend.fun(par.extend = par.extend, data.full = data.full, data.mask = data.mask, 
                              dims = dims, fulllist.par = fulllist.par)

data.full = o.par.extend$data.full
data.mask = o.par.extend$data.mask
data.par = o.par.extend$data.par
name.extend.par = o.par.extend$name.extend.par

o.ss = ss.fun(ss.opts = ss.opts, data.full = data.full, data.ID_mask = data.ID_mask,
              dims = dims, bucket_info = bucket_info, sv = sv_list, fix = fix_list)

data.full = o.ss$data.full
data.ID_mask = o.ss$data.ID_mask
dims = o.ss$dims
bucket_info = o.ss$bucket_info
ss.opts = o.ss$ss.opts

o.CR_SL = CR_SL(cue.rates = cue.rates, survey.length = survey.length,
                bucket_info = bucket_info, dims = dims)

survey.length = o.CR_SL$survey.length
bucket_info = o.CR_SL$bucket_info
mu.rates = o.CR_SL$mu.rates

o.param = param.detfn.fun(sv = sv_list, fix = fix_list, bounds = bounds, name.extend.par = name.extend.par,
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


data.full = numeric_ID(data.full)
data.ID_mask = numeric_ID(data.ID_mask)

data.full = sort.data(data.full, "data.full")
data.mask = sort.data(data.mask, "data.mask")
data.ID_mask = sort.data(data.ID_mask, "data.ID_mask")
data.dists.thetas = sort.data(data.dists.thetas, "data.dists.thetas")

tem = subset(data.full, !is.na(ID))

if("animal_ID" %in% colnames(data.full)){
  dims$n.calls.each.animal = aggregate(tem$ID, list(session = tem$session, animal_ID = tem$animal_ID),
                                       function(x) length(unique(x)))$x
} else dims$n.calls.each.animal = 0

dims$n.detection = cal_n_det(data.full)

if(is.null(ss.opts)){
  n_dir_quadpoints = 0
  n_het_source_quadpoints = 0
  het_source_nodes = 0
  het_source_weights = 0
  cutoff = 0
} else {
  n_dir_quadpoints = ss.opts$n.dir.quadpoints
  n_het_source_quadpoints = ss.opts$n.het.source.quadpoints
  het_source_nodes = ss.opts$het.source.nodes
  het_source_weights = ss.opts$het.source.weights
  cutoff = ss.opts$cutoff
}


#####################################################################

library(TMB)
library(MASS)



data <- list(n_sessions = dims$n.sessions,
             n_animals = dims$n.animals,
             n_IDs = dims$n.IDs,
             #because in the data.full, if one session has no detection
             #it will still record n.trap's rows, a special n.id for
             #looking up the index of row in data.full is necessary
             n_IDs_for_datafull = ifelse(dims$n.IDs == 0, 1, dims$n.IDs),
             n_traps = dims$n.traps,
             n_masks = dims$n.masks,
             n_detection = dims$n.detection,
             n_calls_each_animal = dims$n.calls.each.animal,
             
             nrow_data_full = nrow(data.full),
             nrow_data_mask = nrow(data.mask),
             nrow_dx = nrow(data.dists.thetas),
             nrow_id_mask = nrow(data.ID_mask),
             
             A = A,
             survey_length = survey.length,
             sound_speed = sound.speed,
             
             #code of het_method:
             #1:NULL, 2:GH, 3:rect
             het_method = numeric_het_method(ss.opts$het.source.method),
             n_dir_quadpoints = n_dir_quadpoints,
             n_het_source_quadpoints = n_het_source_quadpoints,
             het_source_nodes = het_source_nodes,
             het_source_weights = het_source_weights,
             cutoff = cutoff,
             #code of detfn_index:
             #1:hn, 2:hhn, 3:hr, 4:th, 5:lth, 6:ss, 7:ss_dir, 8:ss_het
             #detfn_index = numeric_detfn(detfn),
             #temporarily for test purpose~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             detfn_index = numeric_detfn("hn"),
             
             
             #simply use the row index of 'data.par' as the numeric subtitution of parameters
             #1:g0, 2:sigma; 3:lambda0, 4:z, 5:shape.1, 6:shape.2, 7:shpe, 8:scale, 9:b0.ss,
             #10:b1.ss, 11:b2.ss, 12:sigma.ss, 13:kappa, 14:alpha, 15:sigma.toa,
             #16:sigma.b0.ss, 17:D
             #remeber in TMB, index begins from 0, so all of the index above should minus 1
             
             
             #param_og = (1:nrow(data.par))[which(data.par[['par']] %in% param.og)],
             #temporarily for test purpose~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             param_og = c(1, 2, 17),
             
             par_n_col = as.matrix(data.par[, c('n_col_full', 'n_col_mask')]),
             #link index: 1 - identity, 2 - log, 3 - logit
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

#to avoid "." in .cpp file
fulllist.par.4cpp = gsub("\\.", "_", fulllist.par)


parameters = vector('list', length(fulllist.par))
names(parameters) = fulllist.par.4cpp

for(i in 1:length(fulllist.par)){
  name.r = fulllist.par[i]
  name.cpp = fulllist.par.4cpp[i]
  data[[paste0(name.cpp, "_DX")]] = ifelse(is.na(DX.full[[name.r]]), 0, DX.full[[name.r]])
  #column in the design matrix for intercept should be always 1
  data[[paste0(name.cpp, "_DX")]][,1] = ifelse(data[[paste0(name.cpp, "_DX")]][,1] == 0, 1, 
                                               data[[paste0(name.cpp, "_DX")]][,1])
  data[[paste0(name.cpp, "_DX_mask")]] = ifelse(is.na(DX.mask[[name.r]]), 0, DX.mask[[name.r]])
  data[[paste0(name.cpp, "_bound")]] = bounds.input[[name.r]]
  parameters[[name.cpp]] = sv.input[[name.r]]
}

parameters$u = numeric(sum(dims$n.IDs * dims$n.masks * dims$n.traps))

#set the "map" argument for fixed parameters

name.fixed.par = names(fix.input)[!sapply(fix.input, is.null)]


if(length(name.fixed.par) == 0){
  map = list()
} else {
  map = vector('list', length(name.fixed.par))
  names(map) = name.fixed.par
  for(i in name.fixed.par) map[[i]] = NA
}




compile("fit_ascr.cpp")
dyn.load(dynlib("fit_ascr"))

obj <- MakeADFun(data = data, parameters = parameters, DLL="fit_ascr")
obj$hessian <- TRUE
opt = nlminb(obj$par, obj$fn, obj$gr)
o = sdreport(obj)
summary(o, "report")
