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
test_data("bearing fitting")

#the par.extend below is for "Multi-session models"
#par.extend = list(data = list(session = data.frame(session = 1:2, weather = c('rain', 'sunny')),
#                              trap = data.frame(session = rep(1:2, c(6,3)), trap = c(1:6, 1:3), 
#                                                brand = sample(c('a','b','c'), size = 9, replace = T))),
#                  model = list(g0 = "~weather + brand"))

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
ss.link = ss.opts$ss.link

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


###########################################uid stuffs

tem = extract_unique_id(data.full, dims)
#these two dimensions are a little bit ambiguous
#n.id.uid means the number of IDs under each unique capture history, this is a vector that
#simply combines these numbers for all unique capture histories in all sessions in a row
#for example, if there are 3 sessions and 1st session has 3 unique capture histories, and 2nd has no detection
#and the 3rd session has 2 unique capture histories
#then this n.id.uid will be a vector with length of 3 + 0 + 2 = 5

dims$n.id.uid = tem$n_id_uid

#n.uids is the number of unique capture histories in each session, it is a vector
#with length of n.sessions

dims$n.uids = tem$n_uids

#each unique capture history
data_u_bin = tem$data_u_bin
uid_dup_data_full = tem$uid_dup_data_full

DX.full.uid = lapply(DX.full, function(x) {x = x[!uid_dup_data_full,,drop = FALSE]; x = x[order(data_u_bin$session, data_u_bin$u_id, data_u_bin$trap),,drop = FALSE]; return(x)})

data_u_bin = sort.data(data_u_bin, "data_u_bin")

#this one is like the number of detection for each uid, like 'n.detection', this is a literately vector
#if 3 sessions, and their n.uid is 3,0,2, then this "n.detection.uid" is vector with length of 5
dims$n.detection.uid = cal_n_det(data_u_bin, is_uid = TRUE)

#this is a literately vector as well, "n.detection.uid" must be used to extract the index of traps for one uid
index_traps_uid = aggregate(data_u_bin$bincapt, list(session = data_u_bin$session, 
                                                     u_id = data_u_bin$u_id), function(x) which(x == 1))
index_traps_uid = sort.data(index_traps_uid, "index_traps_uid")
index_traps_uid = do.call('c', index_traps_uid$x)


#ID and unique capture history match table
u_id_match = tem$u_id_match
u_id_match = sort.data(u_id_match, "u_id_match")


###########################################end of uid stuffs
#on Fri, pack all of uid stuffs to a function
#and arrange the input for DX.full.uid in TMB


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

#dist is a little special because it could appear in log(capt_dist) in the likelihood
#capt_dist = ifelse(is.na(data.full$dist), 0, data.full$dist)
#capt_dist = ifelse(capt_dist == 0, 1e-20, capt_dist)

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
             n_detection_uid = dims$n.detection.uid,
             n_calls_each_animal = dims$n.calls.each.animal,
             n_uid_session = dims$n.uids,
             n_ids_each_uid = dims$n.id.uid,
             index_traps_uid = index_traps_uid, 
             
             
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
             detfn_index = numeric_detfn(detfn),
             
             
             #simply use the row index of 'data.par' as the numeric subtitution of parameters
             #1:g0, 2:sigma; 3:lambda0, 4:z, 5:shape.1, 6:shape.2, 7:shpe, 8:scale, 9:b0.ss,
             #10:b1.ss, 11:b2.ss, 12:sigma.ss, 13:kappa, 14:alpha, 15:sigma.toa,
             #16:sigma.b0.ss, 17:D
             #remeber in TMB, index begins from 0, so all of the index above should minus 1
             
             
             param_og = (1:nrow(data.par))[which(data.par[['par']] %in% param.og)],
             
             par_n_col = as.matrix(data.par[, c('n_col_full', 'n_col_mask')]),
             #link index: 1 - identity, 2 - log, 3 - logit
             par_link = numeric_link(data.par[['link']]),
             #link index: 1 - identity, 2 - log, 4 - spherical
             ss_link = numeric_link(ss.link),
             
             is_animalID = as.numeric("animal_ID" %in% colnames(data.full)),
             is_ss = as.numeric("ss" %in% bucket_info),
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
             
             u_id_match = u_id_match$ID,
             capt_bin_uid = data_u_bin$bincapt,
             
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
param.og.4cpp = gsub("\\.", "_", param.og)

parameters = vector('list', length(fulllist.par) + 1)
names(parameters) = c(fulllist.par.4cpp, "u")

for(i in 1:length(fulllist.par)){
  name.r = fulllist.par[i]
  name.cpp = fulllist.par.4cpp[i]
  data[[paste0(name.cpp, "_DX")]] = ifelse(is.na(DX.full[[name.r]]), 0, DX.full[[name.r]])
  #column in the design matrix for intercept should be always 1
  data[[paste0(name.cpp, "_DX")]][,1] = ifelse(data[[paste0(name.cpp, "_DX")]][,1] == 0, 1, 
                                               data[[paste0(name.cpp, "_DX")]][,1])
  data[[paste0(name.cpp, "_DX_mask")]] = ifelse(is.na(DX.mask[[name.r]]), 0, DX.mask[[name.r]])
  data[[paste0(name.cpp, "_DX_uid")]] = ifelse(is.na(DX.full.uid[[name.r]]), 0, DX.full.uid[[name.r]])
  data[[paste0(name.cpp, "_bound")]] = bounds.input[[name.r]]
  parameters[[name.cpp]] = sv.input[[name.r]]
}

parameters$u = numeric(sum(dims$n.IDs))

#set the "map" argument for fixed parameters
#obtain the fixed parameter that assigned by the user
name.fixed.par = names(fix.input)[!sapply(fix.input, is.null)]
name.fixed.par.4cpp = gsub("\\.", "_", name.fixed.par)

#combine this with the parameters that will not be used in this model
name.fixed.par.4cpp = c(name.fixed.par.4cpp, fulllist.par.4cpp[!(fulllist.par.4cpp %in% param.og.4cpp)])

#in case there is any duplication, delete them
name.fixed.par.4cpp = unique(name.fixed.par.4cpp)

if(!("ss.het" %in% bucket_info)){
  name.fixed.par.4cpp = c(name.fixed.par.4cpp, "u")
}

map = vector('list', length(name.fixed.par.4cpp))
names(map) = name.fixed.par.4cpp

for(i in 1:length(name.fixed.par.4cpp)){
  par_name = name.fixed.par.4cpp[i]
  map[[par_name]] = factor(rep(NA, length(parameters[[par_name]])))
}




compile("sp_bearing.cpp")
dyn.load(dynlib("sp_bearing"))

system.time({
  obj <- MakeADFun(data = data, parameters = parameters[param.og], map = list(g0 = factor(NA)), DLL="sp_bearing")

  obj$hessian <- TRUE
  opt = nlminb(obj$par, obj$fn, obj$gr)

  o = sdreport(obj)
})

summary(o)
opt$objective
