get_D_tmb = function(fit){
  
  #estimated population density for each mask is recorded in fit$D.mask
  #extract the values from it and make it into a usable format data frame
  #with session and mask indices
  
  dims = fit$output.tmb$dims
  output = data.frame(session = rep(1:dims$n.sessions, dims$n.masks))
  tem = vector('list', dims$n.sessions)
  for(s in 1:dims$n.sessions){
    tem[[s]] = data.frame(mask = 1:dims$n.masks[s], D_est = fit$D.mask[[s]])
  }
  tem = do.call('rbind', tem)
  output = cbind(output, tem)
  row.names(output) = NULL
  data_mask = get_data_mask(fit)
  output = merge(output, data_mask, by = c('session', 'mask'))
  output = output[order(output$session, output$mask),]
  return(output)
}

get_data_full = function(fit){
  output = fit$output.tmb$data.full
  return(output)
}

get_data_mask = function(fit){
  output = fit$output.tmb$data.mask
  return(output)
}

get_data_trap = function(fit){
  output = fit$output.tmb$data.traps
  return(output)
}

get_par_extend = function(fit){
  output = fit$args$par.extend
  return(output)
}

get_par_extend_name = function(fit){
  output = fit$output.tmb$param.extend
  return(output)
}

get_area_unit_tmb = function(fit){
  output = fit$output.tmb$area_unit
  return(output)
}

get_dims_tmb = function(fit){
  output = fit$output.tmb$dims
  return(output)
}

get_cue_rates = function(fit){
  output = fit$output.tmb$avg_cue_rates
  return(output)
}

get_survey_length = function(fit){
  output = fit$args$survey.length
  return(output)
}

get_dist_theta = function(fit){
  output = fit$output.tmb$data.dists.thetas
  return(output)
}

get_detfn = function(fit){
  output = fit$output.tmb$detfn
  return(output)
}

get_DX_full = function(fit){
  output = fit$output.tmb$DX$DX_non_mask
  return(output)
}

get_DX_mask = function(fit){
  output = fit$output.tmb$DX$DX_mask
  return(output)
}

get_data_param = function(fit){
  output = fit$output.tmb$param.info.table
  return(output)
}

get_coef = function(fit){
  output = fit$output.tmb$coef_link
  return(output)
}

get_param_og = function(fit){
  output = fit$output.tmb$param.og
  return(output)
}

get_ss_link = function(fit){
  output = fit$output.tmb$ss.link
  return(output)
}

get_cutoff = function(fit){
  output = fit$output.tmb$cutoff
  return(output)
}

get_buffer = function(fit){
  output = numeric(fit$n.sessions)
  for(s in 1:fit$n.sessions){
    buffer = attr(fit$args$mask[[s]], 'buffer')
    if(!is.null(buffer))
    output[s] = buffer
  }
  return(output)
}

######################################################################################
#get detection history, the key component of simulation
get_det_history = function(fit, data_density){
  source('support_functions.r', local = TRUE)
  source('detfn_tmb.r', local = TRUE)
  data_full = get_data_full(fit)
  data_mask = get_data_mask(fit)
  data_trap = get_data_trap(fit)
  #get design matrices for non-mask level and mask level
  #indicated as "full" and "mask" respectively
  DX_full = get_DX_full(fit)
  DX_mask = get_DX_mask(fit)
  #get the small data frame which records the link function
  #for each parameter and the numbers of covariates for each
  #parameter in non-mask level (in "DX_full") and mask level
  #(in "DX_mask")
  param_table = get_data_param(fit)
  #get all estimated values for each parameter's covariates
  param = get_coef(fit)
  det_fn = get_detfn(fit)
  #get original parameters' names for detection function
  #and extra information (will be used in another function)
  #since "D" is not used in detection function nor extra info
  #remove it
  param_name = get_param_og(fit)
  param_name = param_name[-which(param_name == 'D')]
  
  dims = get_dims_tmb(fit)
  #because in the data.full, if one session has no detection
  #it will still record n.trap's rows, a special n.id is necessary
  n_ID_data_full = ifelse(dims$n.IDs == 0, 1, dims$n.IDs)
  
  #create a new data.full but this time the columns are
  #back-transformed values for each parameters
  data_full_param_value = data_full[, c('session', 'ID', 'trap')]
  data_mask_param_value = data_mask[, c('session', 'mask')]
  
  #get the data frame records distance and angles between each mask and trap
  data_dist_theta = get_dist_theta(fit)
  
  for(i in param_name){
    tem = subset(param_table, param_table$par == i)
    link = tem$link
    #n_col essentially means number of covariates
    n_col_full = tem$n_col_full
    n_col_mask = tem$n_col_mask
    
    #get the design matrices
    tem_DX_full = DX_full[[i]]
    tem_DX_mask = DX_mask[[i]]
    
    #get the values for covariates and make it into a matrix with dim of k*1
    tem_coef_full = matrix(param[[i]][1:n_col_full], ncol = 1)
    
    #design matrix multiple parameter vector for non-mask level to get the
    #non-mask level part of parameters values 
    data_full_param_value[[i]] = tem_DX_full %*% tem_coef_full
    #since 'ID' is not extendable, shrink this data set
    tem = data_full_param_value[!duplicated(data_full_param_value[, c('session', 'trap')]),
                                c('session', 'trap', i)]
    #merge the data with distance and parameters values into one data set
    data_dist_theta = merge(data_dist_theta, tem, by = c('session', 'trap'))
    
    #if there is any mask level covariates, do the same operation with non-mask level
    #covariates to obtain the parameters values in mask level part
    if(n_col_mask > 0){
      tem_coef_mask = matrix(param[[i]][(n_col_full + 1) : (n_col_full + n_col_mask)], 
                             ncol = 1)
      data_mask_param_value[[paste0(i, '_mask')]] = tem_DX_mask %*% tem_coef_mask
      
      tem = data_mask_param_value[, c('session', 'mask', paste0(i, '_mask'))]
      tem = merge(data_dist_theta, data_mask_param_value, by = c('session', 'mask'))
      #combine the two parts to get the final parameters values
      data_dist_theta[[i]] = tem[[i]] + tem[[paste0(i, '_mask')]]
    }
    
    #based on link function, back transform the values of parameters
    if(link == 'log') data_dist_theta[[i]] = exp(data_dist_theta[[i]])
    if(link == 'logit') data_dist_theta[[i]] = exp(data_dist_theta[[i]]) / (exp(data_dist_theta[[i]]) + 1)
  }
  
  #merge the distances, parameters values and simulated number of calls together
  data_dist_theta = merge(data_dist_theta, data_density[,c('session', 'mask', 'n_calls')],
                          by = c('session', 'mask'))
  
  #shrink the data set by picking out the masks with no simulated call
  data_capt = subset(data_dist_theta, data_dist_theta$n_calls != 0)
  
  #sort it firstly, in order to assign ID to each call after the "split" operation much easier
  data_capt = data_capt[order(data_capt$session, data_capt$mask, data_capt$trap),]
  
  #for each n_call > 1, we need to split it to several identical individual call
  tem_gt1 = subset(data_capt, data_capt$n_calls > 1)
  tem_eq1 = subset(data_capt, data_capt$n_calls == 1)
  
  u_n_calls = unique(tem_gt1$n_calls)
  
  tem_gt1_list = vector('list', length(u_n_calls))
  
  for(i in 1:length(u_n_calls)){
    #extract the part with n_calls equals this unique number of calls
    tem = subset(tem_gt1, tem_gt1$n_calls == u_n_calls[i])
    #repeat this part u_n_calls[i] times
    tem_list = vector('list', u_n_calls[i])
    for(j in 1:u_n_calls[i]) tem_list[[j]] = tem
    tem_gt1_list[[i]] = do.call('rbind', tem_list)
  }
  
  tem_gt1 = do.call('rbind', tem_gt1_list)
  
  data_capt = rbind(tem_gt1, tem_eq1)
  #since all calls become individual calls, assign n_call to be 1
  data_capt[['n_calls']] = 1
  
  
  #remove useless objects
  rm(tem_gt1)
  rm(tem_eq1)
  rm(tem)
  rm(tem_list)
  rm(tem_gt1_list)
  
  #since data is well sorted, we could easily assign ID to each call
  data_capt[['ID']] = 0
  for(s in 1:dims$n.sessions){
    index_s = which(data_capt[['session']] == s)
    data_capt[['ID']][index_s] = rep(seq(length(index_s) / dims$n.traps[s]), each = dims$n.traps[s])
  }
  
  #merge the coordination of each trap and each mask to this data frame
  data_capt = merge(data_capt, data_trap, by = c('session', 'trap'), all.x = TRUE)
  
  data_capt = merge(data_capt, data_mask, by = c('session', 'mask'), all.x = TRUE)
  

  tem_list = vector('list', dims$n.sessions)
  
  for(s in 1:dims$n.sessions){
    tem_capt = subset(data_capt, session == s)
    tem_mask = subset(data_mask, session == s)
    
    #calculate the width of each mask
    h_x = unique(diff(sort(unique(tem_mask$x))))
    if(length(h_x) != 1) stop('masks are not evenly splited.')
    #calculate the height of each mask
    h_y = unique(diff(sort(unique(tem_mask$y))))
    if(length(h_y) != 1) stop('masks are not evenly splited.')
    
    #randomly select an coordination from the mask for each ID
    tem_capt$x = tem_capt$x + runif(nrow(tem_capt), -0.5 * h_x, 0.5 * h_x)
    tem_capt$y = tem_capt$y + runif(nrow(tem_capt), -0.5 * h_y, 0.5 * h_y)
    
    #re-calculate distance and theta
    tem_capt$dx = sqrt((tem_capt$trap_x - tem_capt$x) ^ 2 + 
                          (tem_capt$trap_y - tem_capt$y) ^ 2)
    
    tem_capt$theta = bearings_by_vec(tem_capt$trap_x, tem_capt$trap_y,
                                      tem_capt$x, tem_capt$y)
    
    #about the theta part, there is one point needs attention. For e.g.,
    #when the original theta is 1.469693, then the new theta after taking
    #a random location within a mask is 4.612040, although they seems 
    #to be different, but if we take tan(theta), we could see they are
    #quite similar, one is 9.857147 and the other one is 9.931751
    
    tem_list[[s]] = tem_capt
  }
  
  data_capt = do.call('rbind', tem_list)
  
  
  dx = data_capt[['dx']]
  
  #create variables for each parameter will be used in det_fn part
  det_par = vector('list', length(param_name))
  names(det_par) = param_name
  for(i in param_name) det_par[[i]] = data_capt[[i]]
  #calculate detection probability
  if(det_fn != 'ss'){
    data_capt[['det_prob']] = det_prob(det_fn, det_par, dx)
  } else {
    ss.link = get_ss_link(fit)
    cut_off = get_cutoff(fit)
    #ss is a little different, if the simulated ss > cut off, the detection probability is 1
    #and 0 otherwise
    mu = det_prob(det_fn, det_par, dx, ss.link)
    data_capt[['ss']] = rnorm(nrow(data_capt), mu, det_par[["sigma.ss"]])
    data_capt[['det_prob']] = ifelse(data_capt[['ss']] > cut_off, 1, 0)
  }
  
  #remove used variables to release RAM
  rm(det_par)
  
  data_capt[['n_det']] = 0
  
  #at least we need some detection history, otherwise we cannot proceed with any further procedure
  while(sum(data_capt[['n_det']]) == 0){
    data_capt[['n_det']] = rbinom(nrow(data_capt), data_capt[['n_calls']], data_capt[['det_prob']])
  }
  return(data_capt)
}

get_sound_speed = function(fit){
  output = fit$output.tmb$sound.speed
  return(output)
}

get_fit_type = function(fit){
  output = fit$fit.types
  return(output)
}


get_extra_info = function(fit, data_capt){
  fit_types = get_fit_type(fit)
  #check which rows have detection, only simulate for those rows
  index_det = which(data_capt[['n_det']] > 0)
  n_rows = length(index_det)
  if(fit_types['toa']){
    sound_speed = get_sound_speed(fit)
    data_capt[['toa']] = 0
    data_capt[['toa']][index_det] = rnorm(n_rows, data_capt[['dx']][index_det] / sound_speed, 
                                          data_capt[['sigma.toa']][index_det])
  }
  
  
  if(fit_types['dist']){
    data_capt[['dist']] = 0
    shape = data_capt[['alpha']][index_det]
    scale = data_capt[['dx']][index_det] / shape
    data_capt[['dist']][index_det] = rgamma(n_rows, shape = shape, scale = scale)
  }
  
  if(fit_types['bearing']){
    
    data_capt[['bearing']] = 0
    #function 'rvm' does not support vector computation
    if(length(get_coef(fit)$kappa) > 1){
      #if kappa is extended, do a loop to each row with a detection
      for(i in index_det){
        data_capt[['bearing']][i] = CircStats::rvm(1, data_capt[['theta']][i], data_capt[['kappa']][i])
      }
    } else {
      #if kappa is not extended, just take the first one as its value
      k = data_capt[['kappa']][1]
      data_capt[['bearing']][index_det] = (data_capt[['theta']][index_det] + 
                                             CircStats::rvm(n_rows, 0, k)) %% (2*pi)
    }
    
  }
  
  return(data_capt)
}

get_gam = function(fit, param, which.level){
  output = fit$output.tmb$gam_output[[param]]
  output = output[[which.level]]
  return(output)
}

get_DX_new_gam = function(mod, newdata){
  newdata = newdata[, names(mod$var.summary), drop = FALSE]
  newdata$gam.resp = 1
  olddata = mod$mf
  dat = rbind(olddata, newdata)
  output = gam(mod$formula, data = dat, fit = FALSE)$X
  output = output[(nrow(olddata) + 1):nrow(dat), , drop = FALSE]
  return(output)
}
