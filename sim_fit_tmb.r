simcapt_from_fit_tmb = function(fit){
  source('get_funs.r', local = T)
  #get pop density for each mask
  data_density = get_D_tmb(fit)
  #area for each mask
  area_unit = get_area_unit_tmb(fit)


  dims = get_dims_tmb(fit)
  
  #calculate lambda of each mask for poisson distribution
  data_density$lambda_unit = data_density$D_est * rep(area_unit, dims$n.masks)
  
  
  data_density$n_animals = 0
  data_density$n_calls = 0
  
  for(s in 1:dims$n.sessions){
    index_data_density_session = which(data_density$session == s)
    #get the lambda for this session, which sums up the lambda for each mask in this session
    lambda_unit_session = data_density$lambda_unit[index_data_density_session]
    lambda_session = sum(lambda_unit_session)
    #sample a random number from Poisson distribution with lambda_session
    pop_session = rpois(1, lambda_session)
    
    #allocate this number into each mask by a multinomial distribution
    #with (lambdas of each mask / lambda_session) as their probabilities
    
    #actually the operation above with Poisson(lambda_session) first and
    #multinomial later is equivalent to sample a number for each mask with
    #Poisson(lambda_each_mask)
    
    tem = as.vector(rmultinom(1, pop_session, prob = lambda_unit_session))
    data_density[index_data_density_session, 'n_animals'] = tem
    
  }
  
  #simulate capture history based on the simulated density and detection function
  data_capt = sim_det_history(fit, data_density)
  
  #simulate extra info (bearing, dist, ss, toa)
  data_capt = sim_extra_info(fit, data_capt)
  
  #take out the observations with no detection
  o = subset(data_capt, data_capt$n_det > 0)
  
  #pack all result into a data frame
  output = data.frame(session = o$session, ID = o$ID, occasion = 1, trap = o$trap)
  fit_types = get_fit_type(fit)
  if(fit_types['toa']){
    output[['toa']] = o$toa
  }
  
  
  if(fit_types['ss']){
    output[['ss']] = o$ss
  }
  
  if(fit_types['dist']){
    output[['dist']] = o$dist
  }
  
  if(fit_types['bearing']){
    output[['bearing']] = o$bearing
  }
  
  output = output[order(output$session, output$ID, output$trap),]
  
  return(output)
  
}

sim_result_agg = function(sim_result){
  tem = sim_result[[1]]
  n.sim = length(sim_result)
  
  output = vector('list', length(tem))
  names(output) = names(tem)
  
  for(i in names(output)){
    output[[i]] = vector('list', length(tem[[i]]))
    names(output[[i]]) = names(tem[[i]])
    for(j in names(output[[i]])) {
      output[[i]][[j]] = numeric(n.sim)
      for(n in 1:n.sim) output[[i]][[j]][n] = sim_result[[n]][[i]][j]
    }
  }
  
  return(output)
}

#input is the aggregated result
sim_result_plot = function(sim_agg, coef_real, pre_dir){
  for(i in names(coef_real)){
    tem_sim = sim_agg[[i]]
    tem_real = coef_real[[i]]
    for(j in names(tem_real)){
      sim_vec = tem_sim[[j]]
      real_val = tem_real[j]
      dir = paste0(pre_dir, "/", j, ".jpeg")
      jpeg(filename = dir, width = 1920, height = 1080)
      hist(sim_vec, main = paste0('simulation of ', j), xlab = j)
      abline(v = real_val, col = 2)
      dev.off()
    }
  }
  return(0)
}
