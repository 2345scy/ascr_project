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




vcov.ascr_tmb <- function(object, types = "linked", pars = NULL, newdata = NULL, ...){
  source('get_funs.r')
  if(types %in% c('all', 'fitted', 'linked', 'derived')) stop("Argument 'types' must be a subset of 'fitted', 'linked', 'derived', 'all'.")
  
  if(!is.null(newdata) & (!"fitted" %in% types)) types = c(types, 'fitted')
  
  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }
  
  cov_og = object$vcov
  name_dim = gsub(" _ ", "\\.",colnames(object$vcov))
  which.derived = which(substr(name_dim, 1, 3) == 'esa')
  name_dim = name_dim[-which.derived]
  cov_derived = cov_og[which.derived, which.derived, drop = FALSE]
  cov_linked = cov_og[-which.derived, -which.derived, drop = FALSE]
  
  name_extend = get_par_extend_name(object)
  #the original name is needed, in "name_og", "D.(Intercept)", "D.x1", "D.xxx" will be regarded as 'D' only
  #and the same to the other extendable parameters
  name_og = character(length(name_dim))
  
  #extract the link functions for all element in the covariance matrix
  df_param = get_data_param(object)
  link_funs = character(length(name_og))
  
  
  for(i in 1:length(name_dim)){
    tem = gsub("\\..+$", "", name_dim[i])
    name_og[i] = tem
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
  
  
  output = vector('list', length(types))
  names(output) = types
  
  for(type in types){
    
    if(type == "linked"){
      index_par = which(name_og %in% pars)
      output[type] = cov_linked[index_par, index_par, drop = FALSE]
    } else if(type == 'derived'){
      output[type] = cov_derived
    } else if(type == 'fitted'){
      #if 'newdata' is not provided, just keep the extended covariates as parameters with identity link function
      if(is.null(newdata)){
        #extract the covariance matrix for 'linked' and link functions for the assigned parameters only
        index_par = which(name_og %in% pars)
        link_funs = link_funs[index_par]
        cov_linked = cov_linked[index_par, index_par, drop = FALSE]
        
        #to be continued...
        
        
      } else {
        
      }
    }
  }
  
  #if only one type be assigned, change the output from a list to a matrix only
  if(length(types) == 1) output = output[[types]]
  
  return(output)
  
}