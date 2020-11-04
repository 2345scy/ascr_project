
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
    bucket_info <- all.types[all.types %in% names(capt)]
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
    bucket_info <- all.types[all.types %in% names(capt[[1]])]
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
      for(k in c('bincapt', bucket_info)) tem.df[[k]] = as.vector(tem[[k]])
      
      if(is.mrds){
        tem$mrds = as.data.frame(tem$mrds, stringsAsFactors = FALSE)
        if(!"ID" %in% colnames(tem$mrds)) tem$mrds$ID = as.character(1:nrow(tem$mrds))
        if("animal_ID" %in% colnames(tem$mrds)){
          tem.mrds = merge(tem.df, tem$mrds, by = c("animal_ID","ID"))
        } else {
          tem.mrds = merge(tem.df, tem$mrds, by = "ID")
        }
        tem.df = tem.df[order(tem.df$session, tem.df$ID, tem.df$trap),]
        tem.mrds = tem.mrds[order(tem.mrds$session, tem.mrds$ID, tem.mrds$trap),]
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
              bucket_info = c(bucket_info, "mrds"[is.mrds])))
  
}


#################################################################################################

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


###################################################################################

mask.fun = function(mask, dims, data.traps, data.full, bucket_info, local){
  
  stopifnot(any(class(mask) == 'list', is.data.frame(mask), is.matrix(mask)))
  n.sessions = dims$n.sessions
  n.traps = dims$n.traps
  n.IDs = dims$n.IDs
  
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
  tem.df = vector('list', n.sessions)
  
  for (i in 1:n.sessions){
    traps = data.traps[data.traps$session == i, c("trap_x", "trap_y")]
    
    if(ncol(mask[[i]]) != 2) {
      stop("each element of 'mask' must be a matrix or data frame with 2 columns")
    }
    
    colnames(mask[[i]]) = c('x', 'y')
    mask[[i]] <- as.matrix(mask[[i]])
    dists[[i]] <- distances(as.matrix(traps), mask[[i]])
    
    n.rows = n.traps[i] * n.masks[i]
    tem.dist = data.frame(mask = rep(1:n.masks[i], each = n.traps[i]), 
                          trap = rep(1:n.traps[i], n.masks[i]),
                          dx = as.vector(dists[[i]]))
    
    
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
    attr(mask[[i]], "area") = A[i]
    attr(mask[[i]], "buffer") = buffer[i]
    
    tem.df[[i]] = as.data.frame(mask[[i]], stringsAsFactors = FALSE)
    tem.df[[i]]$mask = 1:nrow(mask[[i]])
    tem.df[[i]]$session = i
    tem.df[[i]] = merge(tem.df[[i]], tem.dist, by = 'mask')
  }
  
  data.mask = do.call('rbind', tem.df)
  
  ################################################################################  
  #deal with local argument
  #sort it first, because later, in aggregate() step, 
  #the result will be sorted automatically, however, u.id generated by unique() is not sorted
  is.animal_ID = "animal_ID" %in% colnames(data.full)
  if(is.animal_ID){
    data.full = data.full[order(data.full$session, data.full$animal_ID, 
                                data.full$ID, data.full$trap), ]
  } else {
    data.full = data.full[order(data.full$session, data.full$ID, data.full$trap), ]
  }
  
  all.which.local <- vector(mode = "list", length = n.sessions)
  all.n.local <- vector(mode = "list", length = n.sessions)
  
  
  if ("mrds" %in% bucket_info){
    local <- TRUE
    for (i in 1:n.sessions){
      
      if(n.IDs[i] > 0){
        mrds = data.full[data.full$session == i, c("mrds_x", "mrds_y")]
        mrds = mrds[!duplicated(mrds),]
        all.which.local[[i]] <- find.nearest.mask(as.matrix(mrds), mask[[i]])
        all.n.local[[i]] <- laply(all.which.local[[i]], length)
        all.which.local[[i]] <- c(all.which.local[[i]], recursive = TRUE)
      }
      
    }
  } else if (local){
    
    for (i in 1:n.sessions){
      if(n.IDs[i] > 0){
        if(is.animal_ID){
          bincapt = data.full[data.full$session == i, c("bincapt","animal_ID", "ID")]
          bincapt = aggregate(bincapt$bincapt, list(animal_ID = bincapt$animal_ID,
                                                    ID = bincapt$ID),
                              function(x) t(x))
          
          bincapt = bincapt[order(bincapt$animal_ID, bincapt$ID),]
          bincapt = bincapt[, -c(1,2)]
        } else {
          bincapt = data.full[data.full$session == i, c("bincapt","ID")]
          bincapt = aggregate(bincapt$bincapt, list(ID = bincapt$ID),
                              function(x) t(x))
          bincapt = bincapt[order(bincapt$ID),]
          bincapt = bincapt[, -1]
        }
        
        all.which.local[[i]] <- find_local(bincapt, dists[[i]], buffer[i])
        all.n.local[[i]] <- laply(all.which.local[[i]], length)
        all.which.local[[i]] <- c(all.which.local[[i]], recursive = TRUE)
      } 
    }
  } else {
    for (i in 1:n.sessions){
      if(n.IDs[i] > 0){
        all.n.local[[i]] <- rep(1, n.IDs[i])
        all.which.local[[i]] <- rep(0, n.IDs[i])
      }
    }
  }
  
  
  #prepare a data set for ID*mask since "local" or not is a ID*mask level variable
  tem.df = vector('list', n.sessions)
  u.id = vector('list', n.sessions)
  #to avoid confusing, u.id is animal_ID-ID if is.animal_ID, and the "ID" column is
  #temporarily to be animal_ID-ID as well until in the later step where we restore it
  for(i in 1:n.sessions){
    n.rows = n.IDs[i] * n.masks[i]
    if(n.rows > 0){
      tem = subset(data.full, session == i)
      if(is.animal_ID){
        u.id[[i]] = unique(paste0(tem[['animal_ID']], '-', tem[['ID']]))
        tem.df[[i]] = data.frame(session = rep(i, n.rows),
                                 ID = rep(u.id[[i]], each = n.masks[i]),
                                 mask = rep(1:n.masks[i], n.IDs[i]))
      } else {
        u.id[[i]] = unique(tem[["ID"]])
        tem.df[[i]] = data.frame(session = rep(i, n.rows),
                                 ID = rep(u.id[[i]], each = n.masks[i]),
                                 mask = rep(1:n.masks[i], n.IDs[i]))
      }
    }
  }
  
  data.IDxmask = do.call('rbind', tem.df)
  
  tem.df = vector('list', n.sessions)
  for(i in 1:n.sessions){
    if(n.IDs[i] > 0){
      if(local){
        
        tem.df[[i]] = data.frame(session = rep(i, sum(all.n.local[[i]])),
                                 ID = rep(u.id[[i]], all.n.local[[i]]),
                                 mask = all.which.local[[i]],
                                 local = rep(1, sum(all.n.local[[i]])),
                                 stringsAsFactors = FALSE)
      } else {
        tem.df[[i]] = data.frame(session = i, local = 1)
      }
    }
  }
  
  tem = do.call('rbind', tem.df)
  
 
  data.local = merge(data.IDxmask, tem, 
                     by=c('session', 'ID'[local], 'mask'[local]), all = TRUE)
  
  if(is.animal_ID){
    data.local[['animal_ID']] = gsub("-(.+)", "", data.local[['ID']])
    data.local[['ID']] = gsub("^(.+)-", "", data.local[["ID"]])
  }
  
  tem = merge(data.local, data.mask, by = c('session', 'mask'), all = TRUE)
  
  tem$local = ifelse(is.na(tem$local), 0, tem$local)
  

  data.full = merge(data.full, tem, by = c('session', 'animal_ID'[is.animal_ID], 
                                             'ID', 'trap'), all = TRUE)

  
  
  
  return(list(mask = mask, dists = dists, n.masks = n.masks, A = A, buffer = buffer,
              local = local, all.which.local=all.which.local, all.n.local= all.n.local,
              data.full = data.full))
}



###############################################################################



par.extend.fun = function(par.extend = par.extend, data.full = data.full, dims = dims,
                          fulllist.par = fulllist.par){
  namelist.par.extend = fulllist.par[-which(fulllist.par == 'sigma.b0.ss')]
  fulllist.par.extend = fulllist.par
  
  n.sessions = dims$n.sessions
  n.IDs = dims$n.IDs
  n.traps = dims$n.traps
  n.masks = dims$n.masks
  
  is.animal_ID = "animal_ID" %in% colnames(data.full)
  
  if(!is.null(par.extend)){
    if(!all(names(par.extend) %in% c('data', 'model', 'link', 'scale'))) {
      stop("'par.extend' only accepts 'data', 'model', 'link' and 'scale' as input")
    }
    if(!all(c('data', 'model') %in% names(par.extend))){
      stop("'data' and 'model' must be provided")
    }
    
    #check 'data' component first
    input_data = par.extend$data
    if(!all(names(input_data) %in% c('session', 'trap', 'ID', 'mask'))){
      stop("only 'session', 'trap', 'ID' or 'mask' level data could be used as input")
    }
    
    data.par.extend = data.frame(session = 1:n.sessions)
    
    #check 'session' level data
    df.s = input_data$session
    if(!is.null(df.s)){
      stopifnot(class(df.s) == 'data.frame')
      name.s = colnames(df.s)
      stopifnot('session' %in% name.s)
      stopifnot(length(name.s) > 1)
      if(any(duplicated(df.s$session))) stop('duplicated "session" input')
      stopifnot(df.s$session %in% unique(data.full$session))
      var.s = name.s[-which(name.s == 'session')]
      data.par.extend = merge(data.par.extend, df.s, by = 'session', all = TRUE)
    } else {
      var.s = NULL
    }
    
    #check 'trap' level data
    df.t = input_data$trap
    if(!is.null(df.t)){
      stopifnot(class(df.t) == 'data.frame')
      name.t = colnames(df.t)
      stopifnot(all(c('session', 'trap') %in% name.t))
      stopifnot(length(name.t) > 2)
      if(any(duplicated(df.t[, c('session', 'trap')]))) stop('duplicated "trap" input')
      stopifnot(all(paste(df.t$session, df.t$trap, sep = '-') %in% 
                      unique(paste(data.full$session, data.full$trap, sep = '-'))))
      var.t = name.t[-which(name.t == 'session'| name.t == 'trap')]
      data.par.extend = merge(data.par.extend, df.t, by = 'session', all = TRUE)
    } else {
      var.t = NULL
    }
    
    #check 'ID' level data
    df.I = input_data$ID
    if(!is.null(df.I)){
      stopifnot(class(df.I) == 'data.frame')
      name.I = colnames(df.I)
      stopifnot(all(c('session', 'animal_ID'[is.animal_ID], 'ID') %in% name.I))
      stopifnot(length(name.I) > (2 + as.numeric(is.animal_ID)))
      if(any(duplicated(df.I[, c('session', 'animal_ID'[is.animal_ID], 'ID')]))) stop('duplicated "ID" input')
      
      if(is.animal_ID){
        stopifnot(all(paste(df.I$session, df.I$animal_ID, df.I$ID, sep = '-') %in% 
                  unique(paste(data.full$session, data.full$animal_ID, data.full$ID, sep = '-'))))
        var.I = name.I[-which(name.I == 'session'| name.I == 'ID' | name.I == 'animal_ID')]
      } else {
        stopifnot(all(paste(df.I$session, df.I$ID, sep = '-') %in% 
                        unique(paste(data.full$session, data.full$ID, sep = '-'))))
        var.I = name.I[-which(name.I == 'session'| name.I == 'ID')]
      }

      data.par.extend = merge(data.par.extend, df.I, by ='session', all = TRUE)
    } else {
      var.I = NULL
    }
    
    #check 'mask' level data
    df.m = input_data$mask
    if(!is.null(df.m)){
      stopifnot(is.list(df.m))
      if(class(df.m) == 'list'){
        stopifnot(length(df.m) == n.sessions)
        stopifnot(sapply(df.m, nrow) == n.masks)
        for(i in 1:n.sessions){
          df.m[[i]][['session']] = i
          df.m[[i]][['mask']] = 1:n.masks[i]
        }
        df.m = do.call('rbind', df.m)
      } else {
        stopifnot(n.sessions == 1)
        stopifnot(nrow(df.m) == n.masks)
        df.m[['session']] = 1
        df.m[['mask']] = 1:n.masks
      }
      name.m = colnames(df.m)
      var.m = name.m[-which(name.m == 'session'| name.m == 'mask')]
      data.par.extend = merge(data.par.extend, df.m, by = 'session', all.x = TRUE)
    } else {
      var.m = NULL
    }
    
    #check whether there is duplicated column names across these  data sets
    
    var.ex = c(var.s, var.t, var.I, var.m)
    
    if(any(duplicated(var.ex))) {
      stop('duplicated column name across different level data sets')
    }
    
    #check 'model'
    
    name.extend.par = names(par.extend$model)
    if(!all(name.extend.par %in% namelist.par.extend)){
      stop('one or more parameters indicated in "model" are not supported to be extended currently')
    }
    
    foo = vector('list', length(name.extend.par))
    
    
    #deal the extend parameters one by one
    for(i in 1:length(name.extend.par)){
      
      #scale default is set as TRUE
      if(is.null(par.extend$scale[[name.extend.par[i]]])){
        is.scale = TRUE
      } else {
        is.scale = par.extend$scale[[name.extend.par[i]]]
      }
      
      
      foo[[i]] = as.formula(paste(c('pseudo_y',
                                    as.character(par.extend$model[[name.extend.par[i]]])), collapse = ""))
      tem.var = all.vars(foo[[i]][[3]])
      if(length(tem.var) == 0) stop('please input a valid formula except "~ 1"')
      tem = as.data.frame(lapply(data.par.extend[, tem.var], is.na), stringsAsFactors = FALSE)
      index = apply(tem, 1, any)
      tem.data = subset(data.par.extend, !index)
      tem.data[['pseudo_y']] = 1
      
      
      #here the name of 'scale' is not checked since it is unnecessary
      for(j in tem.var){
        if(is.scale){
          if(class(tem.data[[j]]) == 'numeric'){
            mu = mean(tem.data[[j]])
            std = sd(tem.data[[j]])
            tem.data[[j]] = (tem.data[[j]] - mu)/std
          }
        }
      }
      
      design.matrix = gam(foo[[i]], data = tem.data, fit = FALSE)$X
      design.matrix = as.data.frame(design.matrix, stringsAsFactors = FALSE)
      colnames(design.matrix) = paste(name.extend.par[i], "_", colnames(design.matrix), sep = " ")
      design.matrix[['session']] = tem.data[['session']]
      if(!is.null(tem.data[['trap']])) design.matrix[['trap']] = tem.data[['trap']]
      if(!is.null(tem.data[['animal_ID']])) design.matrix[['animal_ID']] = tem.data[['animal_ID']]
      if(!is.null(tem.data[['ID']])) design.matrix[['ID']] = tem.data[['ID']]
      if(!is.null(tem.data[['mask']])) design.matrix[['mask']] = tem.data[['mask']]
      data.par.extend = merge(data.par.extend, design.matrix, by = c('session', 
                                                                     'trap'[!is.null(design.matrix[['trap']])],
                                                                     'animal_ID'[!is.null(design.matrix[['animal_ID']])],
                                                                     'ID'[!is.null(design.matrix[['ID']])],
                                                                     'mask'[!is.null(design.matrix[['mask']])]), all = TRUE)
      
    }
    
    data.par.extend = data.par.extend[, -which(colnames(data.par.extend) %in% var.ex)]
    
    data.full = merge(data.full, data.par.extend, by = c('session',
                                                         'trap'[!is.null(data.par.extend[['trap']])],
                                                         'animal_ID'[!is.null(data.par.extend[['animal_ID']])],
                                                         'ID'[!is.null(data.par.extend[['ID']])],
                                                         'mask'[!is.null(data.par.extend[['mask']])]), all = TRUE)
    
    
    for(i in fulllist.par.extend){
      if(!i %in% name.extend.par){
        data.full[[paste0(i, ' _ (Intercept)')]] = 1
      }
    }
  } else {
    for(i in fulllist.par.extend){
      data.full[[paste0(i, ' _ (Intercept)')]] = 1
    }
  }
  
  
  
  #finally, record the link, here we only check whether the name is in the full list of parameters,
  #so technically, user could assign link function to any parameter here, but not recommended
  df.link = data.frame(name = character(0), link = character(0), 
                           stringsAsFactors = FALSE)
  
  
  link = par.extend$link
  
  if(!is.null(link)){
    if(!all(names(link) %in% fulllist.par.extend)){
      warning("one or more parameters denoted in 'link' are ignored as they are invalid.")
    }
  }
  
  
  j = 1 
  for(i in fulllist.par.extend){
    df.link[j, 1] = i
    if(i %in% names(link)){
      df.link[j, 2] = link[[i]]
    } else {
      if(i == "D") df.link[j, 2] = 'log'
      if(i == "g0") df.link[j, 2] = 'logit'
      if(i == "sigma") df.link[j, 2] = 'log'
      if(i == "lambda0") df.link[j, 2] = 'log'
      if(i == "z") df.link[j, 2] = 'log'
      if(i == "shape.1") df.link[j, 2] = 'log'
      if(i == "shape.2") df.link[j, 2] = 'identity'
      if(i == "shape") df.link[j, 2] = 'identity'
      if(i == "scale") df.link[j, 2] = 'log'
      if(i == "b0.ss") df.link[j, 2] = 'log'
      if(i == "b1.ss") df.link[j, 2] = 'log'
      if(i == "b2.ss") df.link[j, 2] = 'log'
      if(i == "sigma.ss") df.link[j, 2] = 'log'
      if(i == "kappa") df.link[j, 2] = 'log'
      if(i == "alpha") df.link[j, 2] = 'log'
      if(i == "sigma.toa") df.link[j, 2] = 'log'
      if(i == "sigma.b0.ss") df.link[j, 2] = 'log'
    }
    j = j + 1
  }
  
  
  
  return(list(data.full = data.full, link = df.link, name.extend.par = name.extend.par))
}



#############################################################################
ss.fun = function(ss.opts, data.full, bucket_info, dims, sv, fix, df.link){
  cutoff <- ss.opts[["cutoff"]]
  ss.link <- ss.opts[["ss.link"]]
  directional <- ss.opts[["directional"]]
  het.source <- ss.opts[["het.source"]]
  het.source.method <- ss.opts[["het.source.method"]]
  n.dir.quadpoints <- ss.opts[["n.dir.quadpoints"]]
  n.het.source.quadpoints <- ss.opts[["n.het.source.quadpoints"]]
  

  fit.ss = "ss" %in% bucket_info
  

  if(fit.ss){
    if(is.null(ss.opts)) stop("Argument 'ss.opts' is missing.")
    if(is.null(cutoff)) stop("The 'cutoff' component of 'ss.opts' must be specified.")
    ## Warning for unexpected component names.
    if (!all(names(ss.opts) %in% c("cutoff", "het.source", "het.source.method", 
                                   "n.het.source.quadpoints", "directional", 
                                   "n.dir.quadpoints", "ss.link"))){
      warning("Components of 'ss.opts' may only consist of \"cutoff\", \"het.source\", 
            \"het.source.method\", \"n.het.source.quadpoints\", \"directional\",  
            \"n.dir.quadpoints\", \"ss.link\"; 
            others are being ignored.")
    }
    
    ## Setting default values for het.source and directional, ss.link.
    
    ## By default, directional calling model is only used if b2.ss appears in sv or fix.
    if (is.null(directional)){
      if (is.null(sv$b2.ss) & is.null(fix$b2.ss)){
        directional <- FALSE
      } else {
        directional <- TRUE
      }
    }
    
    #make detailed and final judgment about ss.het, and settle some basic settings
    if (!directional){
      if (!is.null(sv$b2.ss) | !is.null(fix$b2.ss)){
        warning("As the 'directional' component of 'ss.opts' is FALSE, 
        the values of parameter b2.ss in 'sv' and 'fix' are being ignored")
      }
      n.dir.quadpoints = NULL
    } else {
      if(!is.null(fix$b2.ss)){
        if(fix$b2.ss != 0) {
          bucket_info = c(bucket_info, "ss.dir")
        } else {
          warning("'fix$b2.ss' is zero, 
                  all corresponding parameters of 'directional' are ignored")
          directional = FALSE
          n.dir.quadpoints = NULL
        }
      } else {
        bucket_info = c(bucket_info, "ss.dir")
      }
    }
    
    #after settling down all basic info of ss.dir, set default n.dir.quadpoints
    if (directional){
      if (is.null(n.dir.quadpoints)){
        n.dir.quadpoints <- 8
      }
    } else {
      n.dir.quadpoints <- 1
    }
    
    ## By default, heterogeneity source strength model is only
    ## used if sigma.b0.ss appears in sv or fix, or het.source.method is specified
    if (is.null(het.source)){
      if (all(is.null(sv$sigma.b0.ss),
              is.null(fix$sigma.b0.ss),
              is.null(het.source.method))){
        het.source <- FALSE
      } else {
        het.source <- TRUE
      }
    }

    #make detailed and final judgment about ss.het, and settle some basic settings
    if (het.source){
      if (is.null(het.source.method)){
        het.source.method <- "GH"
      }
      if(!is.null(fix$sigma.b0.ss)){
        if(fix$sigma.b0.ss != 0) {
          bucket_info = c(bucket_info, "ss.het")
        } else {
          warning("'fix$sigma.b0.ss' is zero, 
                  all corresponding parameters of 'het.source' are ignored")
          het.source = FALSE
          het.source.method = NULL
          n.het.source.quadpoints = NULL
        }
      } else {
        bucket_info = c(bucket_info, "ss.het")
      }
    } else {
      ## Fixing sigma.b0.ss to 0 if a heterogeneous source
      ## strength model is not being used.
      if (!is.null(sv$sigma.b0.ss) | !is.null(fix$sigma.b0.ss)){
        warning("As the 'het.source' component of 'ss.opts' is FALSE, 
        the values of the parameter sigma.b0.ss in 'sv' and 'fix' are being ignored")
      }
      het.source.method = NULL
      n.het.source.quadpoints = NULL
    }
    
    #after settling down basic settings, assign default value to n.het...
    if (het.source){
      if (is.null(n.het.source.quadpoints)){
        n.het.source.quadpoints <- 15
      }
    } else {
      n.het.source.quadpoints <- 1
    }
    
    ## Setting het.source.gh.
    if (het.source){
      het.source.gh <- het.source.method == "GH"
    } else {
      het.source.gh <- FALSE
    }
    
    het.source.nodes <- 0
    het.source.weights <- 0
    if(het.source){
      if(het.source.method == "GH"){
        GHd <- gaussHermiteData(n.het.source.quadpoints)
        het.source.nodes <- GHd$x
        het.source.weights <- GHd$w
      }
    }
    
    
    ##################################################################
    #simple check and simple argument adjustment
    if(all(c("ss.dir", "ss.het") %in% bucket_info)){
      stop("Fitting of models with both heterogeneity in source signal strength 
       and directional calling is not yet implemented.")
    }
    
    #set default for ss.link
    if (is.null(ss.link)){
      ss.link <- "identity"
    } else if (!(ss.link %in% c("identity", "log", "spherical"))){
      stop("Component 'ss.link' in 'ss.opts' must be \"identity\", \"log\", or \"spherical\".")
    }
    
    ## Error thrown for model with heterogeneity in source strength and a log-link function.
    if ("ss.het" %in% bucket_info){
      if (ss.link != "identity"){
        stop("Fitting of signal strength models with heterogeneity in source 
         signal strength is only implemented with an identity link function.")
      }
    }
    

    
    
    #data.full$ss = NA observations are generated by the previous par.extend()
    #because the session in which no call is detected still contains useful
    #information such as trap level covariates (e.g. detector's brand) and 
    #session level covariates (e.g. weather).
    #in this case, they are merged into the data.full without bincapt nor ss etc.
    #so we generated NA in that merge step.
    #if one call is detected by trap 1 but not trap 2, "ss" will be recorded as 0
    data.ss = subset(data.full, !is.na(ss))
    data.ss.na = subset(data.full, is.na(ss))
    
    
    #adjust all capture history to zero if "ss" < cutoff, and then remove the observations if
    #there is no capture history for that "ID"
    rem = data.ss$ss < cutoff
    for(i in c('bincapt', 'ss', 'bearing', 'toa', 'dist')){
      data.ss[[i]] = ifelse(rem, 0, data.ss[[i]])
    }
    
    is.animal_ID = !is.null(dims$n.animals)
    if(is.animal_ID){
      keep = aggregate(data.ss$bincapt, list(session = data.ss$session, 
                                             animal_ID = data.ss$animal_ID, 
                                             ID = data.ss$ID), sum)
    } else {
      keep = aggregate(data.ss$bincapt, list(session = data.ss$session, 
                                             ID = data.ss$ID), sum)
    }
    
    keep = subset(keep, x > 0)
    colnames(keep)[3] = "keep"
    data.ss = merge(data.ss, keep, by = c("session", "animal_ID"[is.animal_ID], "ID"), all = TRUE)
    data.ss = subset(data.ss, !is.na(data.ss$keep))
    data.ss = data.ss[,-which(colnames(data.ss) == "keep")]
    data.full = rbind(data.ss, data.ss.na)
    
    
    #refresh "dims"
    if(is.animal_ID){
      new.n.animal_IDs = aggregate(data.ss$animal_ID, list(session = data.ss$session), 
                                   function(x) length(unique(x)))
      new.n.IDs = aggregate(paste0(data.ss$animal_ID, "-", data.ss$ID), list(session = data.ss$session), 
                            function(x) length(unique(x)))
      dims$n.animals[new.n.animal_IDs$session] = new.n.animal_IDs$x
      
    } else {
      new.n.IDs = aggregate(data.ss$ID, list(session = data.ss$session), 
                            function(x) length(unique(x)))
    }
    
    n.removed = sum(dims$n.IDs) - sum(new.n.IDs$x)
    if(n.removed > 0){
      message(n.removed, " capture history entries have no received signal", 
              " strengths above the cutoff and have therefore been removed.\n", 
              sep = "")
    }
    
    dims$n.IDs[new.n.IDs$session] = new.n.IDs$x
    
    #modified ss.opts, contains much condensed and regular information
    ss.opts = list(het.source.method = het.source.method,
                   n.dir.quadpoints= n.dir.quadpoints,
                   n.het.source.quadpoints = n.het.source.quadpoints,
                   het.source.nodes = het.source.nodes,
                   het.source.weights = het.source.weights,
                   ss.link = ss.link)
    
    
  } else {
    if (!is.null(ss.opts)){
      warning("Argument 'ss.opts' is being ignored as a signal 
            strength model is not being fitted.")
    }
    
    ss.opts <- NULL
    
    if(any(!is.null(sv$b0.ss),
           !is.null(sv$b1.ss),
           !is.null(sv$b2.ss),
           !is.null(sv$sigma.ss),
           !is.null(sv$sigma.b0.ss),
           !is.null(fix$b0.ss),
           !is.null(fix$b1.ss),
           !is.null(fix$b2.ss),
           !is.null(fix$sigma.ss),
           !is.null(fix$sigma.b0.ss))){
      warning(paste0("Not 'ss' model, coresponding parameters specified in", 
                     "'sv' or 'fix' will be ignored"))
    }
  }
  
  return(list(data.full = data.full, dims = dims,
              ss.opts = ss.opts, link = df.link, bucket_info = bucket_info))
}

