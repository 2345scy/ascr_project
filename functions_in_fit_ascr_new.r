
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
      n.animals = 0
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
      n.animals = numeric(n.sessions)
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
  data.capt = sort.data(data.capt, "data.full")
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

mask.fun = function(mask, dims, data.traps, data.full, bucket_info, local,
                    sound.speed = sound.speed){
  
  stopifnot(any(class(mask) == 'list', is.data.frame(mask), is.matrix(mask)))
  n.sessions = dims$n.sessions
  n.traps = dims$n.traps
  n.IDs = dims$n.IDs
  is.animal_ID = "animal_ID" %in% colnames(data.full)
  
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
  #dists records the distances between every mask to every trap
  #thetas here records the angle of the line between every mask and trap
  #in the original package, it is called "bearings" but there is already 
  #a function occupies this name, so I changed it
  dists <- vector(mode = "list", length = n.sessions)
  thetas = vector("list", length = n.sessions)
  tem.df = vector('list', n.sessions)
  
  for (i in 1:n.sessions){
    traps = data.traps[data.traps$session == i, c("trap_x", "trap_y")]
    
    if(ncol(mask[[i]]) != 2) {
      stop("each element of 'mask' must be a matrix or data frame with 2 columns")
    }
    
    colnames(mask[[i]]) = c('x', 'y')
    mask[[i]] <- as.matrix(mask[[i]])
    dists[[i]] <- distances(as.matrix(traps), mask[[i]])
    #based on the main procedure, at this stage, we don't know whether it is "ss.dir",
    #since it should be determined at later stage when deal with "ss"
    #so we generate bearings anyway.
    thetas[[i]] = bearings(as.matrix(traps), mask[[i]])

    
    n.rows = n.traps[i] * n.masks[i]
    tem.dist.theta = data.frame(mask = rep(1:n.masks[i], each = n.traps[i]), 
                                trap = rep(1:n.traps[i], n.masks[i]),
                                dx = as.vector(dists[[i]]), 
                                theta = as.vector(thetas[[i]]))
    
    
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
    tem.df[[i]] = merge(tem.df[[i]], tem.dist.theta, by = 'mask')
  }
  
  data.dists.thetas = do.call('rbind', tem.df)

  
  data.mask = data.dists.thetas[, c("session", 'mask', 'x', 'y'), drop = FALSE]
  data.mask = data.mask[!duplicated(data.mask[, c('session', 'mask')]), , drop = F]
  data.mask = sort.data(data.mask, "data.mask")
  data.dists.thetas = data.dists.thetas[, -which(colnames(data.dists.thetas) 
                                                 %in% c('x', 'y'))]
  data.dists.thetas = sort.data(data.dists.thetas, "data.dists.thetas")
  
  ################################################################################  
  #deal with local argument
  #sort it first, because later, in aggregate() step, 
  #the result will be sorted automatically, however, u.id generated by unique() is not sorted
  data.full = sort.data(data.full, "data.full")
  
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
  
  #prepare a data set for ID*mask

  data.ID_mask = vector('list', n.sessions)
  u.id = vector('list', n.sessions)
  
  if(local) {
    bucket_info = c(bucket_info, "local")
    tem.df = vector('list', n.sessions)
  }
  

  #to avoid confusing, u.id is animal_ID-ID if is.animal_ID, and the "ID" column is
  #temporarily to be animal_ID-ID as well until in the later step where we restore it
  for(i in 1:n.sessions){
    if(n.IDs[i] > 0){
      tem = subset(data.full, session == i)
      #be careful that 'unique' will return a vector that not been sorted if the
      #original vector is not sorted. but here data.full has been sorted before,
      #so we could literally use it.
      if(is.animal_ID){
        u.id[[i]] = unique(paste0(tem[['animal_ID']], '-', tem[['ID']]))
      } else {
        u.id[[i]] = unique(tem[["ID"]])
      }
      
      if(local){
        tem.df[[i]] = data.frame(session = rep(i, sum(all.n.local[[i]])),
                                 ID = rep(u.id[[i]], all.n.local[[i]]),
                                 mask = all.which.local[[i]], local = 1,
                                 stringsAsFactors = F)
        data.ID_mask[[i]] = data.frame(session = rep(i, n.IDs[i] * n.masks[i]),
                                       ID = rep(u.id[[i]], each = n.masks[i]),
                                       mask = rep(1:n.masks[i], n.IDs[i]),
                                       stringsAsFactors = F)
      } else {
        data.ID_mask[[i]] = data.frame(session = rep(i, n.IDs[i] * n.masks[i]),
                                       ID = rep(u.id[[i]], each = n.masks[i]),
                                       mask = rep(1:n.masks[i], n.IDs[i]),
                                       local = 1,
                                       stringsAsFactors = F)
      }
      
      #assign 'toa_ssq' in this ID_mask data set
      #data.full must be sorted by 'session', 'ID', 'trap'
      #far before this step, it has already been sorted
      #we take 'order' so careful because 'aggregate' returns
      #results based on sorted 'list(ID = xx)' agrgument
      if('toa' %in% bucket_info){
        if(is.animal_ID){
          toa = aggregate(tem$toa, list(animal_ID = tem$animal_ID, ID = tem$ID),
                          function(x) t(x))$x
        } else {
          toa = aggregate(tem$toa, list(ID = tem$ID), function(x) t(x))$x
        }
        toa.ssq = make_toa_ssq(toa, dists[[i]], sound.speed)
        #based on the designed structure of 'data.ID_mask' and the structure
        #of toa.ssq, it need to be transposed
        data.ID_mask[[i]][['toa_ssq']] = as.vector(t(toa.ssq))
      } else {
        data.ID_mask[[i]][['toa_ssq']] = 0
      }
    }
  }
  
  data.ID_mask = do.call('rbind', data.ID_mask)
  
  if(local){
    tem.df = do.call('rbind', tem.df)
    data.ID_mask = merge(data.ID_mask, tem.df, by = c('session', 'ID', 'mask'), all = T)
    data.ID_mask[['local']] = ifelse(is.na(data.ID_mask[['local']]), 0, data.ID_mask[['local']])
  }

  
  if(is.animal_ID){
    data.ID_mask[['animal_ID']] = gsub("-(.+)", "", data.ID_mask[['ID']])
    data.ID_mask[['ID']] = gsub("^(.+)-", "", data.ID_mask[["ID"]])
    data.ID_mask = sort.data(data.ID_mask, "data.ID_mask")
  } else {
    data.ID_mask = sort.data(data.ID_mask, "data.ID_mask")
  }
  
  return(list(mask = mask, dists = dists, thetas = thetas, n.masks = n.masks, A = A, 
              buffer = buffer, all.which.local = all.which.local, all.n.local= all.n.local,
              data.full = data.full, data.dists.thetas = data.dists.thetas,
              data.mask = data.mask, data.ID_mask = data.ID_mask, bucket_info = bucket_info))
}



###############################################################################

par.extend.fun = function(par.extend = par.extend, data.full = data.full, data.mask = data.mask,
                          dims = dims, fulllist.par = fulllist.par){
  namelist.par.extend = fulllist.par[-which(fulllist.par == 'sigma.b0.ss')]
  
  n.sessions = dims$n.sessions
  n.IDs = dims$n.IDs
  n.traps = dims$n.traps
  n.masks = dims$n.masks
  
  is.animal_ID = "animal_ID" %in% colnames(data.full)
  
  #a list records the number of beta in both data.full and data.mask needed for each parameter
  param_distribution = vector('list', length(fulllist.par))
  names(param_distribution) = fulllist.par
  
  if(!is.null(par.extend)){
    if(!all(names(par.extend) %in% c('data', 'model', 'link', 'scale'))) {
      stop("'par.extend' only accepts 'data', 'model', 'link' and 'scale' as input.")
    }
    if(!all(c('data', 'model') %in% names(par.extend))){
      stop("'data' and 'model' must be provided.")
    }
    
    #check 'data' component first
    
    #temporarily not allow "ID" level extension
    #and "animal_ID" level extension will make it too complicated
    #maybe add it in the future
    
    input_data = par.extend$data
    if(!all(names(input_data) %in% c('session', 'trap', 
                                     #'ID', 
                                     'mask'))){
      #stop("only 'session', 'trap', 'ID' or 'mask' level data could be used as input.")
      stop("only 'session', 'trap', or 'mask' level data could be used as input.")
    }
    
    data.par.non.mask = data.frame(session = 1:n.sessions)
    
    #check 'session' level data
    df.s = input_data$session
    if(!is.null(df.s)){
      stopifnot(class(df.s) == 'data.frame')
      name.s = colnames(df.s)
      stopifnot('session' %in% name.s)
      stopifnot(length(name.s) > 1)
      if(any(duplicated(df.s$session))) stop('duplicated "session" input.')
      stopifnot(df.s$session %in% unique(data.full$session))
      var.s = name.s[-which(name.s == 'session')]
      data.par.non.mask = merge(data.par.non.mask, df.s, by = 'session', all = TRUE)
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
      if(any(duplicated(df.t[, c('session', 'trap')]))) stop('duplicated "trap" input.')
      stopifnot(all(paste(df.t$session, df.t$trap, sep = '-') %in% 
                      unique(paste(data.full$session, data.full$trap, sep = '-'))))
      var.t = name.t[-which(name.t == 'session'| name.t == 'trap')]
      data.par.non.mask = merge(data.par.non.mask, df.t, by = 'session', all = TRUE)
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
      if(any(duplicated(df.I[, c('session', 'animal_ID'[is.animal_ID], 'ID')]))) stop('duplicated "ID" input.')
      
      if(is.animal_ID){
        stopifnot(all(paste(df.I$session, df.I$animal_ID, df.I$ID, sep = '-') %in% 
                        unique(paste(data.full$session, data.full$animal_ID, data.full$ID, sep = '-'))))
        var.I = name.I[-which(name.I == 'session'| name.I == 'ID' | name.I == 'animal_ID')]
      } else {
        stopifnot(all(paste(df.I$session, df.I$ID, sep = '-') %in% 
                        unique(paste(data.full$session, data.full$ID, sep = '-'))))
        var.I = name.I[-which(name.I == 'session'| name.I == 'ID')]
      }
      
      data.par.non.mask = merge(data.par.non.mask, df.I, by ='session', all = TRUE)
    } else {
      var.I = NULL
    }
    
    if(any(c('x', 'y') %in% c(var.I, var.s, var.t))){
      stop("'x' and 'y' are reserved column names, please choose another name.")
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
      data.par.mask = df.m
    } else {
      var.m = c('x', 'y')
      data.par.mask = data.mask
    }
    
    #check whether there is duplicated column names across these  data sets
    
    var.ex = c(var.s, var.t, var.I, var.m)
    
    if(any(duplicated(var.ex))) {
      stop('duplicated column name across different level data sets.')
    }
    
    #check 'model'
    
    name.extend.par = names(par.extend$model)
    if(!all(name.extend.par %in% namelist.par.extend)){
      stop('one or more parameters indicated in "model" are not supported to be extended currently.')
    }
    
    #deal the extend parameters one by one
    for(i in name.extend.par){
      
      #scale default is set as TRUE
      if(is.null(par.extend$scale[[i]])){
        is.scale = TRUE
      } else {
        is.scale = par.extend$scale[[i]]
      }
      
      
      foo = as.formula(paste(c('pseudo_y', as.character(par.extend$model[[i]])), collapse = ""))
      foo_var = all.vars(foo[[3]])
      if(length(foo_var) == 0) stop('please input a valid formula except "~ 1"')
      if(!all(foo_var %in% var.ex)) stop(paste0('one or more specifed explaintory varrable in "model"',
                                                ' is not provided in "data".'))
      
      #sigma_toa is not trap extendable
      if(i == 'sigma_toa' & any(foo_var %in% var.t)) stop("'sigma_toa' is not extendable in 'trap' level.")
      
      tem = formula_separate(foo, var.m)
      
      #"full" here does not mean full information, it does not contain mask information
      
      foo_full = tem$foo_full
      foo_mask = tem$foo_mask
      
      foo_full_var = all.vars(foo_full[[3]])
      foo_mask_var = all.vars(foo_mask[[3]])
      
      if(length(foo_full_var) > 0){
        tem = as.data.frame(lapply(data.par.non.mask[, foo_full_var, drop = FALSE], is.na), stringsAsFactors = FALSE)
        index = apply(tem, 1, any)
        tem.data = subset(data.par.non.mask, !index)
        tem.data[['pseudo_y']] = 1
        
        for(j in foo_full_var){
          if(is.scale){
            if(class(tem.data[[j]]) == 'numeric'){
              mu = mean(tem.data[[j]])
              std = sd(tem.data[[j]])
              if(std == 0) stop(paste0("input of ", j, " is constant."))
              tem.data[[j]] = (tem.data[[j]] - mu)/std
            }
          }
        }
        
        design.matrix = gam(foo_full, data = tem.data, fit = FALSE)$X
        n_col_full = ncol(design.matrix)
        design.matrix = as.data.frame(design.matrix, stringsAsFactors = FALSE)
        colnames(design.matrix) = paste(i, "_", colnames(design.matrix), sep = " ")
        design.matrix[['session']] = tem.data[['session']]
        if(!is.null(tem.data[['trap']])) design.matrix[['trap']] = tem.data[['trap']]
        if(!is.null(tem.data[['animal_ID']])) design.matrix[['animal_ID']] = tem.data[['animal_ID']]
        if(!is.null(tem.data[['ID']])) design.matrix[['ID']] = tem.data[['ID']]
        data.par.non.mask = merge(data.par.non.mask, design.matrix, 
                                  by = c('session', 'trap'[!is.null(design.matrix[['trap']])],
                                         'animal_ID'[!is.null(design.matrix[['animal_ID']])],
                                         'ID'[!is.null(design.matrix[['ID']])]), all = TRUE)
        
      } else {
        #no matter whether this parameter is extended under other levels, the intercept always
        #in this data set. On the other hand, no matter it is extended under 'mask' level,
        #there is no intercept column in data.mask
        data.par.non.mask[[paste0(i, ' _ (Intercept)')]] = 1
        n_col_full = 1
      }
      
      
      if(length(foo_mask_var) > 0){
        tem = as.data.frame(lapply(data.par.mask[, foo_mask_var, drop = FALSE], is.na), stringsAsFactors = FALSE)
        index = apply(tem, 1, any)
        tem.data = subset(data.par.mask, !index)
        tem.data[['pseudo_y']] = 1
        
        for(j in foo_mask_var){
          if(is.scale){
            if(class(tem.data[[j]]) == 'numeric'){
              mu = mean(tem.data[[j]])
              std = sd(tem.data[[j]])
              if(std == 0) stop(paste0("input of ", j, " is constant."))
              tem.data[[j]] = (tem.data[[j]] - mu)/std
            }
          }
        }
        
        #delete the first column, which is intercept
        design.matrix = gam(foo_mask, data = tem.data, fit = FALSE)$X[, -1, drop = FALSE]
        n_col_mask = ncol(design.matrix)
        design.matrix = as.data.frame(design.matrix, stringsAsFactors = FALSE)
        colnames(design.matrix) = paste(i, "_", colnames(design.matrix), sep = " ")
        design.matrix[['session']] = tem.data[['session']]
        design.matrix[['mask']] = tem.data[['mask']]
        data.par.mask = merge(data.par.mask, design.matrix, 
                              by = c('session', 'mask'), all = TRUE)
      } else {
        n_col_mask = 0
      }
      
      param_distribution[[i]] = c(n_col_full, n_col_mask)
    }
    
    
    data.par.non.mask = data.par.non.mask[, -which(colnames(data.par.non.mask) %in% var.ex)]
    data.par.mask = data.par.mask[, -which(colnames(data.par.mask) %in% var.ex)]
    data.full = merge(data.full, data.par.non.mask, by = c('session',
                                                           'trap'[!is.null(data.par.non.mask[['trap']])],
                                                           'animal_ID'[!is.null(data.par.non.mask[['animal_ID']])],
                                                           'ID'[!is.null(data.par.non.mask[['ID']])]), all = TRUE)
    data.mask = merge(data.mask, data.par.mask, by = c('session', 'mask'))
    
    for(i in fulllist.par){
      if(!(i %in% name.extend.par)){
        data.full[[paste0(i, ' _ (Intercept)')]] = 1
        param_distribution[[i]] = c(1, 0)
      }
    }
  } else {
    for(i in fulllist.par){
      data.full[[paste0(i, ' _ (Intercept)')]] = 1
      param_distribution[[i]] = c(1, 0)
    }
    name.extend.par = NULL
  }
  
  #finally, record the link, here we only check whether the name is in the full list of parameters,
  #so technically, user could assign link function to any parameter here, but not recommended
  df.par = data.frame(par = character(0), link = character(0), n_col_full = numeric(0),
                      n_col_mask = numeric(0), stringsAsFactors = FALSE)
  
  link = par.extend$link
  
  if(!is.null(link)){
    if(!all(names(link) %in% fulllist.par)){
      warning("one or more parameters denoted in 'link' are ignored as they are invalid.")
    }
  }
  
  j = 1
  for(i in fulllist.par){
    df.par[j, 'par'] = i
    df.par[j, 'n_col_full'] = param_distribution[[i]][1]
    df.par[j, 'n_col_mask'] = param_distribution[[i]][2]
    if(i %in% names(link)){
      df.par[j, 'link'] = link[[i]]
    } else {
      if(i == "D") df.par[j, 'link'] = 'log'
      if(i == "g0") df.par[j, 'link'] = 'logit'
      if(i == "sigma") df.par[j, 'link'] = 'log'
      if(i == "lambda0") df.par[j, 'link'] = 'log'
      if(i == "z") df.par[j, 'link'] = 'log'
      if(i == "shape.1") df.par[j, 'link'] = 'log'
      if(i == "shape.2") df.par[j, 'link'] = 'identity'
      if(i == "shape") df.par[j, 'link'] = 'identity'
      if(i == "scale") df.par[j, 'link'] = 'log'
      if(i == "b0.ss") df.par[j, 'link'] = 'log'
      if(i == "b1.ss") df.par[j, 'link'] = 'log'
      if(i == "b2.ss") df.par[j, 'link'] = 'log'
      if(i == "sigma.ss") df.par[j, 'link'] = 'log'
      if(i == "kappa") df.par[j, 'link'] = 'log'
      if(i == "alpha") df.par[j, 'link'] = 'log'
      if(i == "sigma.toa") df.par[j, 'link'] = 'log'
      if(i == "sigma.b0.ss") df.par[j, 'link'] = 'log'
    }
    j = j + 1
  }
  
  data.full = sort.data(data.full, "data.full")
  
  data.mask = sort.data(data.mask, "data.mask")
  
  return(list(data.full = data.full, data.mask = data.mask, data.par = df.par,
              name.extend.par = name.extend.par))
}



#############################################################################

ss.fun = function(ss.opts, data.full, data.ID_mask, bucket_info, dims, sv, fix){
  cutoff <- ss.opts[["cutoff"]]
  ss.link <- ss.opts[["ss.link"]]
  directional <- ss.opts[["directional"]]
  het.source <- ss.opts[["het.source"]]
  het.source.method <- ss.opts[["het.source.method"]]
  n.dir.quadpoints <- ss.opts[["n.dir.quadpoints"]]
  n.het.source.quadpoints <- ss.opts[["n.het.source.quadpoints"]]
  
  is.animal_ID = "animal_ID" %in% colnames(data.full)
  
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
    
    #make detailed and final judgment about ss.dir, and settle some basic settings
    if (!directional){
      if (!is.null(sv$b2.ss) | !is.null(fix$b2.ss)){
        warning("Since the 'directional' component of 'ss.opts' is FALSE, 
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
      } else {
        if(het.source.method %in% c("GH", "rect")) stop(paste0("Only 'GH' and 'rect' are ",
                                                               "supported as het.source.method."))
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
    if(het.source){
      if(is.null(n.het.source.quadpoints)){
        n.het.source.quadpoints <- 15
      }
    } else {
      n.het.source.quadpoints <- 1
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
    
    if('ss.dir' %in% bucket_info){
      if(!('toa' %in% bucket_info)){
        stop('"toa" is required for directional signal strength model.')
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
    
    
    if(is.animal_ID){
      keep = aggregate(data.ss$bincapt, list(session = data.ss$session, 
                                             animal_ID = data.ss$animal_ID, 
                                             ID = data.ss$ID), sum)
    } else {
      keep = aggregate(data.ss$bincapt, list(session = data.ss$session, 
                                             ID = data.ss$ID), sum)
    }
    
    keep = subset(keep, x > 0)
    if(nrow(keep) == 0) stop("None of signal received is greater than the cut off threshold.")
    colnames(keep)[3] = "keep"
    data.ss = merge(data.ss, keep, by = c("session", "animal_ID"[is.animal_ID], "ID"), all = TRUE)
    data.ss = subset(data.ss, !is.na(data.ss$keep))
    
    if(nrow(data.ss) == 0) stop('no signal strength is greater than the cutoff.')
    
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
      #a logical value to record whether there is any ID been removed
      ID.removed = TRUE
    } else ID.removed = FALSE
    
    dims$n.IDs[new.n.IDs$session] = new.n.IDs$x
    
    if(ID.removed){
      #when any ID is deleted due to the cut off of ss
      if("animal_ID" %in% colnames(data.full)){
        data.ID_mask = subset(data.ID_mask, paste(session, animal_ID, ID, sep = "-") %in%
                                unique(paste(data.full$session, data.full$animal_ID,data.full$ID, sep = "-")))
      } else {
        data.ID_mask = subset(data.ID_mask, paste(session, ID, sep = "-") %in%
                                unique(paste(data.full$session, data.full$ID, sep = "-")))
      }
    }
    
    #modified ss.opts, contains much condensed and regular information
    ss.opts = list(het.source.method = het.source.method,
                   n.dir.quadpoints = n.dir.quadpoints,
                   n.het.source.quadpoints = n.het.source.quadpoints,
                   het.source.nodes = het.source.nodes,
                   het.source.weights = het.source.weights,
                   ss.link = ss.link,
                   cutoff = cutoff)
    
    
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
  
  data.full = sort.data(data.full, "data.full")

  
  return(list(data.full = data.full, dims = dims, data.ID_mask = data.ID_mask,
              ss.opts = ss.opts, bucket_info = bucket_info))
}

############################################################################
CR_SL = function(cue.rates = cue.rates, survey.length = survey.length,
                 bucket_info = bucket_info, dims = dims){
  n.sessions = dims$n.sessions
  
  if (is.null(survey.length)){
    survey.length <- rep(1, n.sessions)
    if (!is.null(cue.rates)){
      stop("The use of `cue.rates' without `survey.length' is no longer supported.
           Please provide `survey.length', and ensure `cue.rates' is measured in the same time units.")
    }
  } else {
    if (length(survey.length) != n.sessions){
      stop("The argument `survey.length' must have a value for each session.")
    }
  }
  
  if (!is.null(cue.rates)){
    bucket_info = c(bucket_info, "freqs")
    mu.rates <- mean(cue.rates)
  } else {
    mu.rates = 1
  }
  
  return(list(bucket_info = bucket_info, mu.rates = mu.rates, survey.length = survey.length))
  
}


############################################################################

param.detfn.fun = function(sv = sv, fix = fix, bounds = bounds, name.extend.par = name.extend.par,
                           detfn = detfn, data.full = data.full, data.mask = data.mask, 
                           data.par = data.par, ss.opts = ss.opts, bucket_info = bucket_info,
                           fulllist.par = fulllist.par, A = A, buffer = buffer,
                           survey.length = survey.length, dims = dims){
  #check bounds, sv, fix
  if(!is.null(bounds)){
    if(class(bounds) != 'list') stop ('"bounds" must be a list.')
    if(!all(names(bounds) %in% fulllist.par)) stop('one or more element in "bounds" is not valid.')
    if(!all(sapply(bounds, function(x) nrow(as.matrix(x))) == 2)){
      stop('each element of "bounds" should be a vector with two values, or a matrix with two rows.')
    } 
  } else {
    bounds = vector('list', 0)
  }
  
  if(!is.null(sv)){
    if(class(sv) != 'list') stop('"sv" must be a list.')
    if(!all(names(sv) %in% fulllist.par)) stop('one or more element in "sv" is not valid.')
  } else {
    sv = vector('list', 0)
  }
  
  if(!is.null(fix)){
    if(class(fix) != 'list') stop('"fix" must be a list.')
    if(!all(names(fix) %in% fulllist.par)) stop('one or more element in "fix" is not valid.')
  } else {
    fix = vector('list', 0)
  }
  
  
  
  #detfn is the basis of choosing parameters, settle this first
  
  if(is.null(detfn)){
    #the detfn in 'ss' is a little different from original, may be needed
    #to be restored as original at very end of 'fit' function
    #ss_het only support ss.link == 'identity', but this has been checked in ss.fun
    if("ss" %in% bucket_info){
      if("ss.het" %in% bucket_info){
        detfn = 'ss_het'
      } else if("ss.dir" %in% bucket_info){
        detfn = 'ss_dir'
      } else {
        detfn = 'ss'
      }
    } else {
      detfn = "hn"
    }
  } else {
    if("ss" %in% bucket_info){
      if(detfn != 'ss') {
        warning('signal strenth is provided, therefore the specified "detfn" is ignored')
      }
      
      if("ss.het" %in% bucket_info){
        detfn = 'ss_het'
      } else if("ss.dir" %in% bucket_info){
        detfn = 'ss_dir'
      } else {
        detfn = 'ss'
      }
      
    } else {
      if(detfn == 'ss'){
        warning('signal strenth is not provided, therefore the specified "detfn" is ignored')
        detfn = 'hn'
      } else {
        if(!detfn %in% c('hn', 'hhn', 'hr', 'th', 'lth')){
          warning('invalid "detfn" is specified, the default will be used instead')
          detfn = 'hn'
        }
      }
    }
  }
  
  #then, decide the original names of the parameters will be input into TMB
  if('ss' %in% bucket_info){
    param.og = c('b0.ss', 'b1.ss', 'sigma.ss')
    if('ss.dir' %in% bucket_info){
      param.og = c(param.og, 'b2.ss')
    }
    
    if('ss.het' %in% bucket_info){
      param.og = c(param.og, 'sigma.b0.ss')
    }
    
  } else {
    if(detfn == 'hn'){
      param.og = c('g0', 'sigma')
    }
    if(detfn == 'hhn'){
      param.og = c('lambda0', 'sigma')
    }
    if(detfn == 'hr'){
      param.og = c('g0', 'sigma', 'z')
    }
    if(detfn == 'th'){
      param.og = c('shape', 'scale')
    }
    if(detfn == 'lth'){
      param.og = c('shape.1', 'shape.2', 'scale')
    }
  }
  
  if('bearing' %in% bucket_info){
    param.og = c(param.og, 'kappa')
  }
  
  if('dist' %in% bucket_info){
    param.og = c(param.og, 'alpha')
  }
  
  if('toa' %in% bucket_info){
    param.og = c(param.og, 'sigma.toa')
  }
  
  param.og = c(param.og, "D")
  
  #generate all DATA_MATRIX for every each parameter
  #and in the loop, we also do some basic conflict checking
  data.full = sort.data(data.full, "data.full")
  
  
  design.matrices.full = vector('list', length(fulllist.par))
  names(design.matrices.full) = fulllist.par
  
  tem = strsplit(colnames(data.full), " _ ")
  tem = lapply(tem, function(x) x[[1]])
  clean.names.full = do.call('c', tem)
  
  design.matrices.mask = vector('list', length(fulllist.par))
  names(design.matrices.mask) = fulllist.par
  
  tem = strsplit(colnames(data.mask), " _ ")
  tem = lapply(tem, function(x) x[[1]])
  clean.names.mask = do.call('c', tem)
  
  
  for(i in fulllist.par){
    index_par = which(data.par$par == i)
    
    if(data.par[index_par, "n_col_mask"] > 0){
      design.matrices.mask[[i]] = as.matrix(data.mask[, which(clean.names.mask == i), drop = FALSE])
    } else {
      design.matrices.mask[[i]] = matrix(0, ncol = 2, nrow = 2) 
    }
    
    design.matrices.full[[i]] = as.matrix(data.full[, which(clean.names.full == i), drop = FALSE])
    
    if(i %in% param.og){
      #deal with "fix" first
      if(i %in% name.extend.par){
        #a valid extended parameter should not be fixed
        if(!is.null(fix[[i]])){
          warning(paste0("parameter ", i, " is going to be extended, therefore, it cannot be fixed."))
          fix[[i]] = NULL
        }
      } else {
        #if this valid parameter is not expended, check whether its fixed value is in the range
        if(!is.null(fix[[i]])){
          tem = param.range.validate(i, fix[[i]])
          if(!tem) stop(paste0("fixed value for parameter ", i, " is out of valid range."))
        }
      }
      #then deal with "sv" and "bounds"
      #if a parameter is extended, the start value will be its intercept beta's start value
      #the bounds of an extended parameter will be assigned to its intercept as well
      if(!is.null(bounds[[i]])){
        #the lower and upper bounds all should be in the valid range of this parameter
        tem1 = param.range.validate(i, ifelse(bounds[[i]][1] == 0, 1e-15, bounds[[i]][1]))
        tem2 = param.range.validate(i, ifelse(bounds[[i]][2] == 0, 1e-15, bounds[[i]][2]))
        if(!all(tem1, tem2)) stop(paste0("bounds for parameter ", i, " is out of valid range."))
        if(bounds[[i]][2] <= bounds[[i]][1]) stop(paste0("invalid range for parameter ", i))
        
        if(!is.null(sv[[i]])){
          tem = param.range.validate(i, sv[[i]])
          if(!tem) stop(paste0("start value for parameter ", i, " is out of valid range."))
          
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            stop(paste0("the start value of ", i, " is out of the specified corresponding bounds."))
          }
        } else {
          #if sv is not assigned, assign a default value
          sv[[i]] = default.sv(i, info = list(data.full = data.full, data.mask = data.mask,
                                              buffer = buffer, A = A, ss.opts = ss.opts, 
                                              detfn = detfn, dims = dims, 
                                              survey.length = survey.length, param.og = param.og,
                                              sv = sv))
          
          #since the bound is provided, we need the default start value is located within the bound
          #if the default start value is out of the specified bounds, re-assign the start value
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            sv[[i]] = 0.5 * (bounds[[i]][1] + bounds[[i]][2])
          }
        }
        
        #whether the fixed value is inside the assigned bound
        #since we have check the validation of fixed value, no need to check it here like sv
        #if bound is NULL, no need to check fixed value against the defaulted assigned bound
        if(!is.null(fix[[i]])){
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            stop(paste0("the start value of ", i, " is out of the specified corresponding bounds."))
          }
        }
        
      } else {
        bounds[[i]] = default.bounds(i)
        
        if(!is.null(sv[[i]])){
          tem = param.range.validate(i, sv[[i]])
          if(!tem) stop(paste0("start value for parameter ", i, " is out of valid range."))
          #I'm not sure if the assigned start value be out of the default bound, what we should do
          #discuss with Ben later.
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            warning(paste0("specified start value of '", i, "' is extreme, a carefull checking ",
                           "is recommended."))
            bounds[[i]][1] = -2 * abs(sv[[i]])
            bounds[[i]][2] = 2 * abs(sv[[i]])
          }
          
        } else {
          #if sv is not assigned, assign a default value

          sv[[i]] = default.sv(i, info = list(data.full = data.full, data.mask = data.mask,
                                              buffer = buffer, A = A, ss.opts = ss.opts, 
                                              detfn = detfn, dims = dims, 
                                              survey.length = survey.length, param.og = param.og,
                                              sv = sv))
        }
      }
    } else {
      #a start value is needed for every parameter regardless it is used or not
      if(!is.null(sv[[i]])){
        warning(paste0("parameter ", i, " assigned in argument 'sv' will be ignored because ",
                       "the detect function is incompatible"))
      }
      sv[[i]] = 0.5
      
      #if this parameter is not used in the model and in the "fix", remove it from "fix"
      #by assigning NULL to it
      if(!is.null(fix[[i]])){
        warning(paste0("parameter ", i, " assigned in argument 'fix' will be ignored because ",
                       "the detect function is incompatible"))
        fix[[i]] = NULL
      }
      
      if(!is.null(bounds[[i]])){
        warning(paste0("parameter ", i, " assigned in argument 'bounds' will be ignored because ",
                       "the detect function is incompatible"))
      }
      
      bounds[[i]] = c(0.3, 0.7)
      
      #if this not used parameter specified to be extended, everything refer to extension will be ignored
      if(i %in% name.extend.par){
        warning(paste0("parameter ", i, " specified to be extended, however, this parameter is not used ",
                       "in the model, everything about extension will be ignored"))
        #correct all of the relevent information
        design.matrices.full[[i]] = design.matrices.full[[i]][, 1, drop = FALSE]
        design.matrices.mask[[i]] = matrix(0, ncol = 2, nrow = 2)
        data.par[index_par, 'n_col_mask'] = 0
        data.par[index_par, 'n_col_full'] = 1
        name.extend.par = name.extend.par[-which(name.extend.par == i)]
      }
      
      
    }
  }
  
  
  #construct the parameters input, which must be transformed
  #they are not needed to be named, especially the betas for extended parameters
  
  #the sv.input contains vectors as element in the list
  #in the vector, the first 'n_col_full' numbers are for the relevant columns in
  #data.full, and the rest for the relevant columns in data.mask
  #the fit.input contains scalars as element in the list, because extended param
  #has no "fit"
  sv.input = vector('list', length(fulllist.par))
  names(sv.input) = fulllist.par
  fix.input = vector('list', length(fulllist.par))
  names(fix.input) = fulllist.par
  #the bounds.input contains matrices as each element in the list
  bounds.input = vector('list', length(fulllist.par))
  names(bounds.input) = fulllist.par
  
  
  for(i in fulllist.par){
    index_par = which(data.par$par == i)
    length_param = data.par[index_par, 'n_col_full'] + data.par[index_par, 'n_col_mask']
    link = data.par[index_par, 'link']
    
    
    #here is a small problem need attention, if the user specify a parameter to be fixed to zero,
    #and its link function is 'log', than, the real input of fixed value in TMB will be log(1e-20)
    #it is not exactly zero. But this problem could be solved by rounding the output to 7 digits,
    #so in the output, this parameter will be zero again.
    if(!is.null(fix[[i]])) {
      fix.input[[i]] = link.fun(link = link, value = fix[[i]])
      #if a parameter is fixed, we set its start value to its fixed value
      #since only not extended parameter could be fixed, so we don't need to
      #worry about its length, it is always a scaler
      sv.input[[i]] = fix.input[[i]]
    } else {
      sv.input[[i]] = numeric(length_param)
      #since the numeric(n) is defaults set to 0, we only need to assign the
      #first element, no matter it is extended or not
      sv.input[[i]][1] = link.fun(link = link, value = sv[[i]])
    }

    
    
    #set the default bounds of extended betas as the starter of our bounds matrix
    #and replace the first one (or the only one) by the assigned bounds' transformation
    if(link == 'identity'){
      bounds.input[[i]] = matrix(rep(c(-1e8, 1e8), length_param),
                                 nrow = 2)
    } else {
      bounds.input[[i]] = matrix(rep(c(-log(1e20), log(1e20)), length_param),
                                 nrow = 2)
    }
    
    bounds.input[[i]][1, 1] = link.fun(link = link, value = bounds[[i]][1])
    bounds.input[[i]][2, 1] = link.fun(link = link, value = bounds[[i]][2])
  }
  
  return(list(design.matrices.full = design.matrices.full, 
              design.matrices.mask = design.matrices.mask,
              name.extend.par = name.extend.par,
              data.par = data.par,
              sv.input = sv.input,
              fix.input = fix.input,
              bounds.input = bounds.input,
              detfn = detfn,
              param.og = param.og))
  
}


