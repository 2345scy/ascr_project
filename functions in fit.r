
capture.fun = function(capt){
  all.types <- c("bearing", "dist", "ss", "toa")

  if("bincapt" %in% names(capt)){
    n.sessions = 1
    is.animalID = "animal_ID" %in% colnames(capt$bincapt)
    n.traps = ncol(capt$bincapt)
    #if "is.animal_ID", 2 columns were used for "animal_ID" and "ID", so deduct 2
    if(is.animalID) {
      n.traps = n.traps - 2
      n.animals = length(unique(capt$bincapt$animal_ID))
    } else {
      n.animals = NULL
    }
    n.IDs = nrow(capt$bincapt)
    info.types <- all.types[all.types %in% names(capt)]
    is.mrds = "mrds" %in% names(capt)
    
  } else {
    n.sessions = length(capt)
    is.animalID = "animal_ID" %in% colnames(capt[[1]]$bincapt)
    n.traps = sapply(capt, function(x) ncol(x$bincapt))
    if(is.animalID) {
      n.traps = n.traps - 2
      n.animals = sapply(capt, function(x) length(unique(x$bincapt$animal_ID)))
    } else {
      n.animals = NULL
    }
    n.IDs = sapply(capt, function(x) nrow(x$bincapt))
    info.types <- all.types[all.types %in% names(capt[[1]])]
    is.mrds = "mrds" %in% names(capt[[1]])
    
  }
  
  tem.data.capt = vector('list', n.sessions)
  for(i in 1:n.sessions){
    if(n.sessions == 1){
      tem = capt
    } else {
      tem = capt[[i]]
    }
    
    number.row = n.IDs[i] * n.traps[i]
    tem.df = data.frame(session = rep(i, number.row), ID = character(number.row))
    if(is.animalID) tem.df$animal_ID = character(number.row)
    for(j in c("trap", "bincapt", all.types, "mrds_x", "mrds_y")) tem.df[[j]] = numeric(number.row)
    
    if(number.row > 0){
      if(is.animalID){
        tem.df$animal_ID = rep(tem$bincapt$animal_ID, n.traps[i])
        tem.df$ID = rep(tem$bincapt$ID, n.traps[i])
        tem$bincapt = as.matrix(tem$bincapt[, 1:n.traps[i]])
      } else {
        tem.df$ID = rep(rownames(tem$bincapt), n.traps[i])
      }
      
      tem.df$trap = rep(1:n.traps[i], each = n.IDs[i])
      for(k in c('bincapt', info.types)) tem.df[[k]] = as.vector(tem[[k]])
      
      if(is.mrds){
        tem$mrds = as.data.frame(tem$mrds, stringsAsFactors = FALSE)
        if(!"ID" %in% colnames(tem$mrds)) tem$mrds$ID = as.character(1:nrow(tem$mrds))
        if("animal_ID" %in% colnames(tem$mrds)){
          tem.mrds = merge(tem.df, tem$mrds, by = c("animal_ID","ID"))
        } else {
          tem.mrds = merge(tem.df, tem$mrds, by = "ID")
        }
        
        tem.df$mrds_x = tem.mrds$V1
        tem.df$mrds_y = tem.mrds$V2
      }
      
    }
    tem.data.capt[[i]] = tem.df
  }
  data.capt = do.call("rbind", tem.data.capt)
  data.capt = data.capt[order(data.capt$session, data.capt$ID, data.capt$trap),]
  return(list(data.capt = data.capt, dims = list(n.sessions = n.sessions, n.traps = n.traps,
                                                 n.IDs = n.IDs, n.animals = n.animals),
              info.types = c(info.types, "mrds"[is.mrds])))
  
}




#################################################################################################
#requirement of "traps" is the same as "create.capt" function. It is possible to add any features
#refer to traps only into "traps", but the first 2 columns must be coordinates
#"dims" is one of the output component of "capture.fun"
trap.fun = function(traps, dims){
  stopifnot(any(class(traps) == 'list', is.data.frame(traps), is.matrix(traps)))
  n.sessions = dims$n.sessions
  n.traps = dims$n.traps
  #make traps to a list with length == n.sessions
  if(any(class(traps) == "list" & length(traps) == 1, is.data.frame(traps), is.matrix(traps))){
    tem.traps = vector('list', n.sessions)
    for(i in 1:n.sessions){
      if(class(traps) == 'list') tem.traps[[i]] = traps[[1]]
      if(is.data.frame(traps) | is.matrix(traps)) tem.traps[[i]] = traps
    }
    traps = tem.traps
  } else {
    stopifnot(length(traps) == n.sessions) 
  }
  stopifnot(all(sapply(traps, nrow) == n.traps), all(sapply(traps, ncol) == 2))
  
  
  for(i in 1:n.sessions) {
    traps[[i]] = as.data.frame(traps[[i]])
    colnames(traps[[i]]) = c("trap_x", "trap_y")
    traps[[i]]$session = i
    traps[[i]]$trap = 1:nrow(traps[[i]])
  }
  return(do.call("rbind", traps))
}

###############################################################################


par.extend.fun = function(par.extend = par.extend, data.full = data.full){
  namelist.par.extend = c('g0', 'sigma', 'lambda0', 'zeta', 'shape1', 
                          'shape2', 'scale', 'ss.beta0', 'ss.beta1',
                          'ss.beta2', 'ss.sigma', 'kappa', 'alpha', 'sigma.toa')
  if(!is.null(par.extend)){
    stopifnot('data' %in% names(par.extend))
    input_data = par.extend$data
    par.extend.name = names(par.extend)[-which(names(par.extend) == 'data')]
    if(length(par.extend.name) == 0) stop("please provide at least one parameter for extension")
    stopifnot(all(par.extend.name %in% namelist.par.extend))
    
    
    #######################################################
    #check input
    stopifnot(all(sapply(input_data,class) == 'data.frame') & class(input_data) == 'list')
    if(!all(names(input_data) %in% c('session', 'trap', 'ID', 'animal_ID'))){
      stop('Currently only "session", "trap", "ID", "animal_ID" are acceptable data input')
    }
    #check the component of data
    for(k in names(input_data)){
      if(!all(c("session", k) %in% colnames(input_data[[k]]))){
        tem = ifelse(k == 'session', " ", paste0(" and '", k, "' "))
        msg = paste0("In the data of ", k, ", columns of 'session'", tem, "must be provided")
        stop(msg)
      }
    }
    #check other components one by one
    for(i in par.extend.name){
      input = par.extend[[i]]
      stopifnot(class(input) == 'list')
      if(!all(c("model", "link") %in% names(input))){
        stop("input for each extended parameter must contains 'model' and 'link'")
      }
      if(!input$link %in% c('identical', 'log', 'logit')){
        stop('Currently only "identical", "log", "logit" are supported as link function for extended parameters')
      }
      
      
      
      #check whether the variables in the 'model' argument are all provided by the data
      formula.check = strsplit(unique(trimws(strsplit(input$model, "\\+|\\*|\\:")[[1]])), "_")
      for(k in 1:length(formula.check)){
        stopifnot(length(formula.check[[k]]) > 1)
        tem = formula.check[[k]]
        names(formula.check)[k] = tem[1]
        formula.check[[k]] = tem[2:length(tem)]
      }
      
      if(!all(names(formula.check) %in% names(input_data))){
        stop("please provide data for all information included in the argument 'model'")
      }
      
      for(k in names(formula.check)){
        if(!all(formula.check[[k]] %in% colnames(input_data[[k]]))){
          msg = paste0("Please provide variable(s) of '", paste(formula.check[[k]], collapse = "', '"), 
                       "' in the data of '",k, "'")
          stop(msg)
        }
      }
      
    }
    
    #################################################
    #deal with data first
    
    is.animalID = "animal_ID" %in% colnames(data.full)
    tem.data = data.full[,c("session", "animal_ID"[is.animalID], 'trap', 'ID')]
    input_data = par.extend$data
    
    for(k in names(input_data)){
      tem = which(!names(input_data[[k]]) %in% c('session', k))
      names(input_data[[k]])[tem] = paste0(k,'_', names(input_data[[k]])[tem])
      tem.data = merge(tem.data, input_data[[k]], by = c('session', k[!k == 'session']), all = TRUE)
      
    }
    tem.data$pseudo_y = 1
    
    df.for.link = data.frame(name = par.extend.name, link = character(length(par.extend.name)),
                             stringsAsFactors = FALSE)
    
    for(i in par.extend.name){
      input = par.extend[[i]]
      #extract the columns index which will be used for this extended parameter
      tem.index = which(colnames(tem.data) %in% unique(trimws(strsplit(input$model, '\\+|\\*|\\:')[[1]])))
      #get the subset which will be used for this extended parameter
      tem = subset(tem.data,apply(!is.na(tem.data[,tem.index]), 1, all))
      foo = paste0('pseudo_y~', input$model)
      X = gam(as.formula(foo), data = tem, fit = FALSE)$X
      X = as.data.frame(X)
      colnames(X) = paste0(i, '_', colnames(X))
      X = cbind(tem[, c('session', 'trap', 'animal_ID'[is.animalID], 'ID')], X)
      data.full = merge(data.full, X, by = c('session', k[!k == 'session'], 'trap', 'ID'), all.x = TRUE)
      df.for.link[which(df.for.link$name == i), "link"] = input$link
      
    }
    
    par.std.name = namelist.par.extend[!namelist.par.extend %in% par.extend.name]
    
    for(i in par.std.name){
      data.full[[paste0(i, '_(Intercept)')]] = 1
    }

  }
  else{
    for(i in namelist.par.extend){
      data.full[[paste0(i, '_(Intercept)')]] = 1
    }
    
    df.for.link = data.frame(name = character(0), link = character(0), 
                             stringsAsFactors = FALSE)
  }
  
  return(list(data.full = data.full, link = df.for.link))
}


###################################################################################
mask.fun = function(mask, dims, data.traps){
  
  stopifnot(any(class(mask) == 'list', is.data.frame(mask), is.matrix(mask)))
  n.sessions = dims$n.sessions
  n.traps = dims$n.traps
  
  if(any(class(mask) == "list" & length(mask) == 1, is.data.frame(mask), is.matrix(mask))){
    tem.mask = vector('list', n.sessions)
    for(i in 1:n.sessions){
      if(class(mask) == 'list') tem.mask[[i]] = mask[[1]]
      if(is.data.frame(mask) | is.matrix(mask)) tem.mask[[i]] = mask
    }
    mask = tem.mask
  } else {
    stopifnot(length(mask) == n.sessions) 
  }
  
  n.masks <- sapply(mask, nrow)
  
  A <- numeric(n.sessions)
  buffer <- numeric(n.sessions)
  dists <- vector(mode = "list", length = n.sessions)
  for (i in 1:n.sessions){
    traps = data.traps[data.traps$session == i, c("trap_x", "trap_y")]
    
    if(ncol(mask[[i]]) != 2) {
      stop("each element of 'mask' must be a matrix or data frame with 2 columns")
    }
    
    colnames(mask[[i]]) = c('x', 'y')
    
    mask[[i]] <- as.matrix(mask[[i]])
    dists[[i]] <- distances(as.matrix(traps), mask[[i]])
    ## Filling in area and buffer, if missing.
    if (is.null(attr(mask[[i]], "area"))){
      mask.dists <- distances(mask[[i]], mask[[i]])
      A[i] <- min(mask.dists[mask.dists > 0])^2/10000
    } else {
      A[i] <- attr(mask[[i]], "area")
    }
    if (is.null(attr(mask[[i]], "buffer"))){
      buffer[i] <- max(apply(dists[[i]], 2, min))
    } else {
      buffer[i] <- attr(mask[[i]], "buffer")
    }
    attr(mask[[i]], "area") <- A[i]
    attr(mask[[i]], "buffer") <- buffer[i]
  }
  
  
  return(list(mask = mask, dists = dists, n.masks = n.masks, A = A, buffer = buffer))
}
