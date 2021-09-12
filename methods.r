coef.ascr_tmb <- function(object, pars = "fitted", ...){
  if ("all" %in% pars){
    pars <- c("fitted", "derived", "linked")
  }
  if (any(pars == "esa")){
    pars <- pars[-which(pars == "esa")]
    pars <- c(pars, paste("esa", 1:object$n.sessions, sep = "."))
  }
  par.names <- names(object$coefficients)
  if (!all(pars %in% c("fitted", "derived", "linked", "esa", par.names))){
    stop("Argument 'pars' must either contain a vector of parameter names, or a subset of \"fitted\", \"derived\", \"linked\", and \"all\".")
  }
  if (any(c("fitted", "derived", "linked") %in% pars)){
    which.linked <- grep("_link", par.names)
    which.D = which(par.names == "D")
    linked <- object$coefficients[c(which.linked, which.D)]
    #what is this "Da"?
    which.derived <- which(substr(par.names, 1, 3) == "esa" | par.names == "Da")
    derived <- object$coefficients[c(which.derived, which.D)]
    fitted <- object$coefficients[-c(which.linked, which.derived)]
    out <- mget(pars)
    names(out) <- NULL
    out <- c(out, recursive = TRUE)
    if (!object$fit.ihd){
      out <- out[c("D", names(out)[!(names(out) %in% c("D.(Intercept)", "D"))])]
    }
  } else {
    out <- object$coefficients[pars]
  }
  out
}




vcov.ascr_tmb <- function(object, types = "linked", pars = NULL, new_covariates = NULL, ...){
  source('get_funs.r', local = TRUE)
  source('support_functions.r', local = TRUE)
  if(!types %in% c('all', 'fitted', 'linked', 'derived')) stop("Argument 'types' must be a subset of 'fitted', 'linked', 'derived', 'all'.")
  
  if(!is.null(new_covariates) & (!"fitted" %in% types)) types = c(types, 'fitted')
  
  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }
  param_values_og = get_coef(object)
  cov_og = object$vcov
  name_dim_og = colnames(cov_og)
  name_dim = gsub(" _ ", "\\.",name_dim_og)
  which.derived = which(substr(name_dim, 1, 3) == 'esa')
  name_dim = name_dim[-which.derived]
  name_dim_og = name_dim_og[-which.derived]
  cov_derived = cov_og[which.derived, which.derived, drop = FALSE]
  cov_linked = cov_og[-which.derived, -which.derived, drop = FALSE]
  
  name_extend = get_par_extend_name(object)
  
  
  
  if(is.null(name_extend) & !is.null(new_covariates)){
    warning('No parameter is extended, argument "new_covariates" will be ignored.')
    new_covariates = NULL
  }
  
  #the original name is needed, in "name_og", "D.(Intercept)", "D.x1", "D.xxx" will be regarded as 'D' only
  #and the same to the other extendable parameters
  name_og = character(length(name_dim))
  
  #extract the link functions for all element in the linked covariance matrix
  #if new_covariates is provided, this "link_funs" will be useless as cov matrix of fitted and
  #cov matrix of linked will have different dimentions
  df_param = get_data_param(object)
  link_funs = character(length(name_og))
  #extracted linked estimated values, the reason why not doing do.call('c', param_values_og)
  #is that, the fixed parameters could screw this up
  param_values = numeric(length(name_og))
  
  for(i in 1:length(name_dim)){
    tem = gsub("\\..+$", "", name_dim[i])
    name_og[i] = tem
    param_values[i] = param_values_og[[tem]][name_dim_og[i]]
    if(!tem %in% name_extend){
      name_dim[i] = tem
      link_funs[i] = df_param[which(df_param$par == tem), 'link']
    } else {
      link_funs[i] = 'identity'
    }
  }
  
  dimnames(cov_linked) = list(name_dim, name_dim)
  
  #check which parameter will be displayed
  if(is.null(pars)){
    pars = unique(name_og)
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters, which is not fixed, in this model.")
  }
  
  #obtain the indices for the assigned "pars"
  index_par = which(name_og %in% pars)
  cov_linked = cov_linked[index_par, index_par, drop = FALSE]
  link_funs = link_funs[index_par]
  param_values = param_values[index_par]
  
  output = vector('list', length(types))
  names(output) = types
  
  for(type in types){

    
    if(type == "linked"){
      output[[type]] = cov_linked
    } else if(type == 'derived'){
      output[[type]] = cov_derived
    } else if(type == 'fitted'){
      #if 'new_covariates' is not provided, just keep the extended covariates as parameters with identity link function
      if(is.null(new_covariates)){
        output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, link_funs = link_funs)
        dimnames(output[[type]]) = dimnames(cov_linked)
      } else {
        gam.output = get_gam(object)
        output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, new_covariates = new_covariates,
                                               name_og = name_og, name_extend = name_extend, df_param = df_param,
                                               gam.output = gam.output)
        dimnames(output[[type]]) = list(pars, pars)
      }
    }
  }
  
  #if only one type be assigned, change the output from a list to a matrix only
  if(length(types) == 1) output = output[[types]]
  
  return(output)
  
}
