coef.ascr_tmb = function(object, types = 'linked', pars = NULL, new_covariates = NULL, ...){
  source('get_funs.r', local = TRUE)
  source('support_functions.r', local = TRUE)
  if(any(!types %in% c('all', 'fitted', 'linked', 'derived'))) stop("Argument 'types' must be a subset of {'fitted', 'linked', 'derived', 'all'}.")
  
  if(!is.null(new_covariates) & (!"fitted" %in% types)) types = c(types, 'fitted')
  
  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }
  
  is.fitted = 'fitted' %in% types
  is.linked = 'linked' %in% types
  is.derived = 'derived' %in% types
  
  
  #extract some key information for the object
  param_values_og = get_coef(object)
  name_og = get_param_og(object)
  name_extend = get_par_extend_name(object)
  df_param = get_data_param(object)
  esa = get_esa(object)$value
  
  #check which parameter will be displayed
  if(is.null(pars)){
    pars = name_og
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters' name in this model.")
  }
  
  if(is.null(name_extend) & !is.null(new_covariates)){
    warning('No parameter is extended, argument "new_covariates" will be ignored.')
    new_covariates = NULL
  }
  
  output = NULL
  
  for(i in pars){
    values_link = param_values_og[[i]]
    if(!i %in% name_extend){
      names(values_link) = paste(i, "link", sep = "_")
    } else {
      names(values_link) = gsub(" _ ", "\\.", names(values_link))
      names(values_link) = paste(names(values_link), "link", sep = "_")
    }
    
    
    
    if(is.linked){
      output = c(output, values_link)
    }
    
    if(is.fitted){
      par_info = subset(df_param, par == i)
      link = par_info$link
      if(i %in% name_extend){
        if(is.null(new_covariates)){
          #if there is no new covariates assigned and the parameter is extended, then regards all of the
          #relevant beta as identity linked parameters
          values_fitted = values_link
          names(values_fitted) = gsub("_link", "", names(values_fitted))
          link = 'identity'
        } else {
          if(nrow(new_covariates) != 1){
            stop('Argument "new_covariates" can only accept 1 row.')
          } else {
            gam = get_gam(object, i)
            tem = try({values_fitted = get_extended_par_value(gam, par_info$n_col_full, par_info$n_col_mask, values_link, new_covariates)})
            if(is(tem, 'try-error')) stop('Please make sure all covariates needed for assigned "par" are provided. Defaulty all parameters are assigned as "par".')
            names(values_fitted) = i
          }
        }
      } else {
        values_fitted = values_link
        names(values_fitted) = i
      }
      
      values_fitted = unlink.fun(link = link, value = values_fitted)
      
      output = c(output, values_fitted)
    }
    
  }
  
  
  if(is.derived){
    value_esa = esa$value
    names(value_esa) = paste('esa', seq(length(value_esa)), sep = '.')
    output = c(output, value_esa)
  }
  
  class(output) = 'coef_ascr_tmb'
  
  return(output)
}

print.coef_ascr_tmb = function(x){
  output = as.matrix(x)
  colnames(output) = "Est"
  print(output)
}



vcov.ascr_tmb = function(object, types = "linked", pars = NULL, new_covariates = NULL, ...){
  source('get_funs.r', local = TRUE)
  source('support_functions.r', local = TRUE)
  if(any(!types %in% c('all', 'fitted', 'linked', 'derived'))) stop("Argument 'types' must be a subset of {'fitted', 'linked', 'derived', 'all'}.")
  
  if(!is.null(new_covariates) & (!"fitted" %in% types)) types = c(types, 'fitted')
  
  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }
  
  #'og' below means original output from the object
  param_values_og = get_coef(object)
  cov_og = object$vcov
  name_dim_og = colnames(cov_og)
  name_dim = gsub(" _ ", "\\.",name_dim_og)
  which.derived = which(substr(name_dim, 1, 3) == 'esa')
  name_dim = name_dim[-which.derived]
  name_dim_og = name_dim_og[-which.derived]
  cov_derived = cov_og[which.derived, which.derived, drop = FALSE]
  cov_linked = cov_og[-which.derived, -which.derived, drop = FALSE]
  dimnames(cov_linked) = list(name_dim, name_dim)
  name_extend = get_par_extend_name(object)
  
  
  
  if(is.null(name_extend) & !is.null(new_covariates)){
    warning('No parameter is extended, argument "new_covariates" will be ignored.')
    new_covariates = NULL
  }
  
  #the original name is needed, in "name_og", "D.(Intercept)", "D.x1", "D.xxx" will be regarded as 'D' only
  #and the same to the other extendable parameters
  
  #the 'og' here means the original parameters name
  name_og = ori_name(name_dim)
  
  #extract the link functions for all element in the linked covariance matrix
  #if new_covariates is provided, this "link_funs" will be useless as cov matrix of fitted and
  #cov matrix of linked will have different dimentions
  df_param = get_data_param(object)
  link_funs = character(length(name_og))
  #extracted linked estimated values, the reason why not doing do.call('c', param_values_og)
  #is that, the fixed parameters could screw this up
  param_values = numeric(length(name_og))
  
  for(i in 1:length(name_dim)){
    param_values[i] = param_values_og[[name_og[i]]][name_dim_og[i]]
    if(!name_og[i] %in% name_extend){
      name_dim[i] = name_og[i]
      link_funs[i] = df_param[which(df_param$par == name_og[i]), 'link']
    } else {
      link_funs[i] = 'identity'
    }
  }
  
  
  #check which parameter will be displayed
  if(is.null(pars)){
    pars = unique(name_og)
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters, which is not fixed, in this model.")
  }
  
  #obtain the indices for the assigned "pars"
  index_par = which(name_og %in% pars)
  name_og = name_og[index_par]
  cov_linked = cov_linked[index_par, index_par, drop = FALSE]
  link_funs = link_funs[index_par]
  param_values = param_values[index_par]
  
  output = vector('list', length(types))
  names(output) = types
  
  for(type in types){

    
    if(type == "linked"){
      output[[type]] = cov_linked
      tem = paste(name_dim, 'link', sep = '_')
      dimnames(output[[type]]) = list(tem, tem)
    } else if(type == 'derived'){
      output[[type]] = cov_derived
    } else if(type == 'fitted'){
      #if 'new_covariates' is not provided, just keep the extended covariates as parameters with identity link function
      if(is.null(new_covariates)){
        output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, link_funs = link_funs)
        dimnames(output[[type]]) = list(name_dim, name_dim)
      } else {
        
        #standardize new_covariates first
        new_covariates = object$scale.covs(new_covariates)
        
        gam.output = get_gam(object)
        tem = try({output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, new_covariates = new_covariates, pars = pars,
                                                         name_og = name_og, name_extend = name_extend, df_param = df_param,
                                                         gam.output = gam.output)})
        if(is(tem, 'try-error')) stop('Please make sure all covariates needed for assigned "par" are provided. Defaulty all parameters are assigned as "par".')
        dimnames(output[[type]]) = list(pars, pars)
      }
    }
  }
  
  #if only one type be assigned, change the output from a list to a matrix only
  if(length(types) == 1) output = output[[types]]
  
  return(output)
  
}



stdEr.ascr_tmb = function(object, types = 'linked', pars = NULL, new_covariates = NULL, ...){
  output = vcov(object = object, types = types, pars = pars, new_covariates = new_covariates)
  
  if(is(output, 'list')){
    for(i in names(output)) output[[i]] = sqrt(diag(output[[i]]))
    class(output) = c('std_ascr_tmb', 'list')
  } else {
    output = sqrt(diag(output))
    class(output) = 'std_ascr_tmb'
  }
  
  return(output)
}

print.std_ascr_tmb = function(x){
  if(is(x, 'list')){
    for(i in names(x)){
      x[[i]] = as.matrix(x[[i]])
      colnames(x[[i]]) = 'Std'
    }
    class(x) = 'list'
  } else {
    x = as.matrix(x)
    colnames(x) = 'Std'
  }
  
  print(x)
}




confint.ascr_tmb = function(object, types = 'linked', level = 0.95, method = 'default',
                            linked = FALSE, pars = NULL, new_covariates = NULL,
                            qqplot = FALSE, ask = FALSE, ...){
  source("support_functions.r", local = TRUE)
  if(method != 'default') stop('Please apply bootstrap to the output of the model before using other methods.')
  
  if(any(!types %in% c('all', 'fitted', 'linked', 'derived'))) stop("Argument 'types' must be a subset of {'fitted', 'linked', 'derived', 'all'}.")
  
  if(!is.null(new_covariates) & (!"fitted" %in% types)) types = c(types, 'fitted')
  if(linked & (!"fitted" %in% types)) types = c(types, 'fitted')
  
  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }
  
  if(level < 0.5) stop('argument "level" cannot be above 0.5')
  p_upper = 0.5 + 0.5 * level
  p_lower = 0.5 - 0.5 * level
  col_name = paste(c(p_lower, p_upper) * 100, "%", sep = " ")
  q_upper = qnorm(p_upper)
  q_lower = qnorm(p_lower)
  
  #if linked = FALSE, use coef(fit, 'fitted') + 1.96 * stdEr(fit, 'fitted')
  #if linked = TRUE, use unlink.fun(coef(fit, 'linked') + 1.96 * stdEr(fit, 'linked'))

  if('linked' %in% types | linked){
    df_est = vector_to_df(coef(object, types = 'linked', pars = pars))
    df_std = vector_to_df(stdEr(object, types = 'linked', pars = pars))
    df_linked = merge(df_est, df_std, by = 'name', all.x = TRUE)
    colnames(df_linked) = c('par', 'est', 'std')
    df_linked$std = ifelse(is.na(df_linked$std), 0, df_linked$std)
    df_linked$lower = df_linked$est + q_lower * df_linked$std
    df_linked$upper = df_linked$est + q_upper * df_linked$std
  }
  
  output = vector('list', length(types))
  names(output) = types
  
  for(type in types){
    if(type == 'linked'){
      output[[type]] = as.matrix(df_linked[, c('lower', 'upper')])
      rownames(output[[type]]) = df_linked$par
      colnames(output[[type]]) = col_name
    } else if(type == 'fitted'){
      if(linked){
        tem = as.matrix(df_linked[, c('lower', 'upper')])
        
      } else {
        
      }
    } else {
      
    }
  }
  
  #if only one type be assigned, change the output from a list to a matrix only
  if(length(types) == 1) output = output[[types]]
  
  return(output)
}
