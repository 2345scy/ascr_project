
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
  namelist.par.extend = c('g0', 'sigma', 'lambda0', 'z', 'shape.1', 
                          'shape.2', 'shape', 'scale', 'b0.ss', 'b1.ss',
                          'b2.ss', 'sigma.ss', 'kappa', 'alpha', 'sigma.toa')
  
  par.extend.name = character(0)
  
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
      if(!input$link %in% c('identity', 'log', 'logit')){
        stop('Currently only "identity", "log", "logit" are supported as link function for extended parameters')
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
  
  return(list(data.full = data.full, link = df.for.link, par.extend.name = par.extend.name))
}


###################################################################################
mask.fun = function(mask, dims, data.traps, data.full, info.types, local){
  
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
################################################################################  
  #deal with local argument (question remains, needs Ben's help)
  all.which.local <- vector(mode = "list", length = n.sessions)
  all.n.local <- vector(mode = "list", length = n.sessions)
  
  if ("mrds" %in% info.types){
    local <- TRUE
    for (i in 1:n.sessions){
      
      if(dims$n.IDs[i] > 0){
        mrds = data.full[data.full$session == i, c("mrds_x", "mrds_y")]
        mrds = mrds[!duplicated(mrds),]
        all.which.local[[i]] <- find.nearest.mask(as.matrix(mrds), mask.m[[i]])
        all.n.local[[i]] <- laply(all.which.local[[i]], length)
        all.which.local[[i]] <- c(all.which.local[[i]], recursive = TRUE)
      } else {
        
      }
      
    }
  } else if (local){
    is.animal_ID = "animal_ID" %in% colnames(data.full)
    for (i in 1:n.sessions){
      if(dims$n.IDs[i] > 0){
        if(is.animal_ID){
          bincapt = data.full[data.full$session == i, c("bincapt","animal_ID", "ID")]
          bincapt = aggregate(bincapt$bincapt, list(animal_ID = bincapt$animal_ID,
                                                    ID = bincapt$ID),
                              function(x) t(x))
          bincapt = bincapt[, -c(1,2)]
        } else {
          bincapt = data.full[data.full$session == i, c("bincapt","ID")]
          bincapt = aggregate(bincapt$bincapt, list(ID = bincapt$ID),
                              function(x) t(x))
          bincapt = bincapt[, -1]
        }
        
        all.which.local[[i]] <- find_local(bincapt, dists[[i]], buffer[i])
        all.n.local[[i]] <- laply(all.which.local[[i]], length)
        all.which.local[[i]] <- c(all.which.local[[i]], recursive = TRUE)
      } else {
        all.n.local[[i]] = NULL
        all.which.local[[i]] = NULL
      }
    }
  } else {
    for (i in 1:n.sessions){
      if(dims$n.IDs[i] > 0){
        all.n.local[[i]] <- rep(1, dims$n.IDs[i])
        all.which.local[[i]] <- rep(0, dims$n.IDs[i])
      } else {
        all.n.local[[i]] <- NULL
        all.which.local[[i]] <- NULL
      }
    }
  }
  
################################################################################ 
  
  return(list(mask = mask, dists = dists, n.masks = n.masks, A = A, buffer = buffer,
              local = local, all.which.local=all.which.local, all.n.local= all.n.local))
}


###########################################################################

ihd.fun = function(ihd.opts, mask, dims, info.types){
  n.sessions = dims$n.sessions
  n.masks = dims$n.masks
  
  if(!is.null(ihd.opts)){
    if(class(ihd.opts) != 'list'){
      stop("'ihd.opts' must be a list with 'model', 'covariates' and 'scale' (optional)")
    }
    
    if(!all(names(ihd.opts) %in% c('model', 'covariates', 'scale'))){
      stop("'ihd.opts' must be a list with 'model', 'covariates' and 'scale' (optional)")
    }
    
    if(is.null(ihd.opts$model)) stop("'model' is compulsory.")
    
    covariates = ihd.opts$covariates
    model = ihd.opts$model
    cov.scale = ihd.opts$scale
    
    if(is.null(cov.scale)) cov.scale = TRUE
    if(!is.formula(model) | length(as.character(model)) != 2){
      stop("'model' must be a formula begins with '~'")
    }
    
    if(!is.null(covariates)){
      stopifnot(any(class(covariates) == 'list', is.data.frame(covariates), 
                    is.matrix(covariates)))
      
      if(any(class(covariates) == "list" & length(covariates) == 1, 
             is.data.frame(covariates), is.matrix(covariates))){
        tem.covariates = vector('list', n.sessions)
        for(i in 1:n.sessions){
          if(class(covariates) == 'list') tem.covariates[[i]] = covariates[[1]]
          if(is.data.frame(covariates) | is.matrix(covariates)) tem.covariates[[i]] = covariates
        }
        covariates = tem.covariates
      } else {
        stopifnot(length(covariates) == n.sessions) 
      }
      
      if(!all(sapply(covariates, nrow) == n.masks)){
        stop("Each data frame of 'covariates' should provide covariate values at each mask point")
      }
      
      for(i in 1:n.sessions){
        covariates[[i]] = as.data.frame(covariates[[i]], stringsAsFactors = FALSE)
        if(any(c('x','y') %in% colnames(covariates[[i]]))) {
          stop("Covariate names 'x' and 'y' are reserved for x- and y-coordiates.")
        }
        covariates[[i]] = cbind(covariates[[i]], 
                                as.data.frame(mask[[i]], stringsAsFactors = FALSE))
      }
    } else {
      if(model!=~1){
        stop('please provide corresponding covariates')
      }
      covariates = vector('list', n.sessions)
      for(i in 1:n.sessions){
        covariates[[i]] = as.data.frame(mask[[i]], stringsAsFactors = FALSE)
      }
    }
    
    all.covariates = do.call('rbind', covariates)
    ################################################
    
    if(cov.scale){
      for(i in 1:ncol(all.covariates)){
        if(is.numeric(all.covariates[[i]])){
          tem = all.covariates[[i]]
          tem.mean = mean(x)
          tem.sd = sd(x)
          if(tem.sd!=0) all.covariates[[i]] = (tem - tem.mean)/tem.sd
        }
      }
    }
    
    #########################################
    all.covariates$gam.resp = 0
    
    if(model!=~1){
      model.formula = paste0("gam.resp", paste(as.character(model), collapse = ""))
      model.formula = as.formula(paste0(model.formula, "+te(x,y)"))
    } else {
      model.formula = as.formula("gam.resp~te(x,y)")
    }
    
    tem.fit = gam(model.formula, data = all.covariates, fit = FALSE)
    design.matrix = as.data.frame(tem.fit$X)
    colnames(design.matrix) = paste0("D_", colnames(design.matrix))
    
    which.session <- rep(1:n.sessions, times = n.masks)
    data.ihd = vector('list', n.sessions)
    for (i in 1:n.sessions){
      data.ihd[[i]] <- design.matrix[which.session == i, , drop = FALSE]
    }
    
    info.types = c(info.types, "ihd")
    
  } else {
    data.ihd = vector('list', n.sessions)
    for (i in 1:n.sessions){
      data.ihd[[i]] <- data.frame(tem = rep(1, n.masks[i]))
      colnames(data.ihd[[i]])[1] = "D_(Intercept)"
    }
  }
  
  return(data.ihd)
  
}


#############################################################################
ss.fun = function(ss.opts, data.full, info.types, dims, sv, fix, df.link){
  cutoff <- ss.opts$cutoff
  ss.link <- ss.opts$ss.link
  directional <- ss.opts$directional
  het.source <- ss.opts$het.source
  het.source.method <- ss.opts$het.source.method
  n.dir.quadpoints <- ss.opts$n.dir.quadpoints
  n.het.source.quadpoints <- ss.opts$n.het.source.quadpoints
  
  
  fit.ss = "ss" %in% info.types
  
  
  if(fit.ss){
    ## Warning for unexpected component names.
    if (!all(names(ss.opts) %in% c("cutoff", "het.source", "het.source.method", 
                                   "n.het.source.quadpoints", "directional", 
                                   "n.dir.quadpoints", "ss.link"))){
      warning("Components of 'ss.opts' may only consist of \"cutoff\", \"het.source\", 
            \"het.source.method\", \"n.het.source.quadpoints\", \"directional\",  
            \"n.dir.quadpoints\", \"ss.link\"; 
            others are being ignored.")
    }
    
    if(is.null(ss.opts)) stop("Argument 'ss.opts' is missing.")
    if(is.null(cutoff)) stop("The 'cutoff' component of 'ss.opts' must be specified.")
    
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
      sv$b2.ss <- 0
      fix$b2.ss <- 0
      n.dir.quadpoints = NULL
    } else {
      if(!is.null(fix$b2.ss)){
        if(fix$b2.ss != 0) {
          info.types = c(info.types, "ss.dir")
        } else {
          warning("'fix$b2.ss' is zero, 
                  all corresponding parameters of 'directional' are ignored")
          directional = FALSE
          sv$b2.ss <- 0
          n.dir.quadpoints = NULL
        }
      } else {
        info.types = c(info.types, "ss.dir")
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
          info.types = c(info.types, "ss.het")
        } else {
          warning("'fix$sigma.b0.ss' is zero, 
                  all corresponding parameters of 'het.source' are ignored")
          het.source = FALSE
          sv$sigma.b0.ss = 0
          het.source.method = NULL
          n.het.source.quadpoints = NULL
        }
      } else {
        info.types = c(info.types, "ss.het")
      }
    } else {
      ## Fixing sigma.b0.ss to 0 if a heterogeneous source
      ## strength model is not being used.
      if (!is.null(sv$sigma.b0.ss) | !is.null(fix$sigma.b0.ss)){
        warning("As the 'het.source' component of 'ss.opts' is FALSE, 
        the values of the parameter sigma.b0.ss in 'sv' and 'fix' are being ignored")
      }
      sv$sigma.b0.ss <- 0
      fix$sigma.b0.ss <- 0
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
        Hd <- gaussHermiteData(n.het.source.quadpoints)
        het.source.nodes <- GHd$x
        het.source.weights <- GHd$w
      }
    }
    
    
    ##################################################################
    #simple check and simple argument adjustment
    if(all(c("ss.dir", "ss.het") %in% info.types)){
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
    if ("ss.het" %in% info.types){
      if (ss.link != "identity"){
        stop("Fitting of signal strength models with heterogeneity in source 
         signal strength is only implemented with an identity link function.")
      }
    }
    
    #record ss.link
    tem = nrow(df.link)
    df.link[tem + 1, 1] = "ss"
    df.link[tem + 1, 2] = ss.link
    
    
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
      message(n.removed, " capture history entries have no received signal 
            strengths above the cutoff and have therefore been removed.\n", sep = "")
    }
    
    dims$n.IDs[new.n.IDs$session] = new.n.IDs$x
    
    #modified ss.opts, contains much condensed and regular information
    ss.opts = list(het.source.method = het.source.method,
                   n.dir.quadpoints= n.dir.quadpoints,
                   n.het.source.quadpoints = n.het.source.quadpoints,
                   het.source.nodes = het.source.nodes,
                   het.source.weights = het.source.weights)
    
    
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
    sv$b0.ss = 0
    sv$b1.ss = 0
    sv$b2.ss = 0
    sv$sigma.ss = 0
    sv$sigma.b0.ss = 0
    
    fix$b0.ss = 0
    fix$b1.ss = 0
    fix$b2.ss = 0
    fix$sigma.ss = 0
    fix$sigma.b0.ss = 0
  }
  
  return(list(data.full = data.full, dims = dims, sv = sv, fix = fix,
              ss.opts = ss.opts, link = df.link, info.types = info.types))
}

