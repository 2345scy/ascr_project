p.dot.defaultD = function(points = NULL, traps = NULL, detfn = NULL, ss_dir = NULL,
                          ss.link = NULL, pars = NULL, A = A, n.quadpoints = 8){
  
  dists <- distances(traps, points)
  
  if(ss_dir){
    n.traps <- nrow(traps)
    n.points <- nrow(points)
    dirs <- (0:(n.quadpoints - 1))*2*pi/n.quadpoints
    probs <- numeric(n.points)
    ## Integrating over all possible directions.
    ## TODO: Write all this in C++.
    for (i in 1:n.quadpoints){
      dir <- dirs[i]
      bearings <- bearings(traps, points)
      orientations <- abs(dir - bearings)
      for (j in 1:n.points){
        ## Probabilities of detection given orientation.
        o.prob <- numeric(n.traps)
        for (k in 1:n.traps){
          o.prob[k] <- calc.detfn(dists[k, j], detfn, pars, ss.link,
                                  orientations[k, j])
        }
        probs[j] <- probs[j] + (1/n.quadpoints)*(1 - prod(1 - o.prob))
      }
    }
    out <- probs
  } else {
    probs <- calc.detfn(dists, detfn, pars, ss.link)
    out <- aaply(probs, 2, function(x) 1 - prod(1 - x))
  }
  
  
  out <- A*sum(out)
  
  return(out)
}


formula_separate = function(foo, var.m){
  foo_vars = all.vars(foo[[3]])
  response_var = all.vars(foo[[2]])
  index = which(foo_vars %in% var.m)
  if(length(index) == 0){
    full_vars = foo_vars
    mask_vars = character(0)
  } else {
    full_vars = foo_vars[-index]
    mask_vars = foo_vars[index]
  }
  
  foo_terms = attr(terms(foo), 'factors')[-1,,drop = FALSE]
  
  row_names = as.list(rownames(foo_terms))
  col_names = as.list(colnames(foo_terms))
  #use col names to check whether there is interaction term
  #between variables from mask level and other levels
  for(i in 1:length(col_names)){
    tem.name = col_names[[i]]
    tem.foo = as.formula(paste0('psedudo_y~',tem.name))
    tem.foo_vars = all.vars(tem.foo[[3]])
    if(!(sum(tem.foo_vars %in% var.m) %in% c(0, length(tem.foo_vars)))){
      stop("interaction between variables from mask level and other level is not supported.")
    }
  }
  #use row names to separate this foo_terms matrix into two matrix
  #one for mask level and other one for other levels
  for(i in 1:length(row_names)){
    tem.name = row_names[[i]]
    tem.foo = as.formula(paste0('psedudo_y~',tem.name))
    tem.foo_vars = all.vars(tem.foo[[3]])
    row_names[[i]] = tem.foo_vars
  }
  
  len = lapply(row_names, length)
  foo_terms = foo_terms[rep(1:nrow(foo_terms), len),,drop = FALSE]
  row_names = do.call('c', row_names)
  rownames(foo_terms) = row_names
  
  n.var.full = sum(!(foo_vars %in% var.m))
  n.var.mask = sum((foo_vars %in% var.m))
  
  if(n.var.full > 0) {
    tem_full = vector('list', n.var.full)
    names(tem_full) = full_vars
  }
  
  if(n.var.mask > 0) {
    tem_mask = vector('list', n.var.mask) 
    names(tem_mask) = mask_vars
  }
  for(i in foo_vars){
    if(i %in% full_vars){
      tem_full[[i]] = foo_terms[which(row_names == i),,drop = FALSE]
      tem_full[[i]] = apply(tem_full[[i]], 2, function(x) sum(x) > 0)
      tem_full[[i]] = matrix(tem_full[[i]], ncol = ncol(foo_terms),
                             byrow = T)
      colnames(tem_full[[i]]) = colnames(foo_terms)
      rownames(tem_full[[i]]) = i
    } else {
      tem_mask[[i]] = foo_terms[which(row_names == i),,drop = FALSE]
      tem_mask[[i]] = apply(tem_mask[[i]], 2, function(x) sum(x) > 0)
      tem_mask[[i]] = matrix(tem_mask[[i]], ncol = ncol(foo_terms),
                             byrow = T)
      colnames(tem_mask[[i]]) = colnames(foo_terms)
      rownames(tem_mask[[i]]) = i
    }
  }
  
  if(n.var.full > 0){
    tem = apply(do.call('rbind', tem_full), 2, any)
    tem.terms = names(tem)[which(tem)]
    foo_full = as.formula(paste0(response_var, "~", paste(tem.terms, collapse = "+")))
  } else {
    foo_full = as.formula(paste0(response_var, "~ 1"))
  }
  
  if(n.var.mask > 0){
    tem = apply(do.call('rbind', tem_mask), 2, any)
    tem.terms = names(tem)[which(tem)]
    foo_mask = as.formula(paste0(response_var, "~", paste(tem.terms, collapse = "+")))
  } else {
    foo_mask = as.formula(paste0(response_var, "~ 1"))
  }
  
  return(list(foo_full = foo_full, foo_mask = foo_mask))
}


#if the value is valid, return TRUE

param.range.validate = function(param, value){
  if(param %in% c("sigma", "lambda0", "z", "shape.1","scale", "sigma.b0.ss", "kappa", "alpha", "sigma.toa",
                  "b0.ss", "b1.ss", "b2.ss", "sigma.ss", 'D')){
    if(value > 0){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if(param == 'g0'){
    if(value > 0 & value <= 1){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }
}


#default sv, this part needs read the original function carefully

default.sv = function(param, info){
  data.full = info[['data.full']]
  data.mask = info[['data.mask']]
  buffer = info[['buffer']]
  A = info[['A']]
  ss.opts = info[['ss.opts']]
  cutoff = ss.opts[['cutoff']]
  ss.link = ss.opts[["ss.link"]]
  detfn = info[['detfn']]
  dims = info[["dims"]]
  survey.length = info[["survey.length"]]
  param.og = info[['param.og']]
  sv = info[['sv']]
  
  
  
  if(param == 'sigma'){
    easy.out = mean(buffer) / 4
    is.animal_ID = "animal_ID" %in% colnames(data.full)
    same.traplocs = (nrow(data.full[!duplicated(data.full[,c('trap_x', 'trap_y')]), ]) == 1)
    if(same.traplocs){
      return(easy.out)
    } else {
      tem = data.full[!duplicated(data.full[,c("session", "animal_ID"[is.animal_ID], "ID", "trap")]), ]
      tem = subset(tem, !is.na(bincapt) & bincapt == 1)
      u.id = unique(tem[, c('session', "animal_ID"[is.animal_ID], "ID")])
      u.id = paste(u.id$session, u.id$animal_ID[is.animal_ID], u.id$ID, sep = "_")
      
      mean.dists = data.frame(avg_dist = numeric(0), weight = numeric(0))
      j = 1
      for(i in 1:length(u.id)){
        tem1 = subset(tem, paste(tem$session, tem$animal_ID[is.animal_ID], tem$ID, sep = "_") == u.id[i],
                      select = c("trap_x", "trap_y"))
        
        if(nrow(tem1) > 1){
          rc.dists <- distances(as.matrix(tem1), as.matrix(tem1))
          mean.dists[j, 1] = mean(rc.dists[rc.dists > 0])
          mean.dists[j, 2] = length(rc.dists[rc.dists > 0])
          j = j + 1
        }
      }
      
      if(nrow(mean.dists) > 0){
        out <- sum(mean.dists$avg_dist * mean.dists$weight)/sum(mean.dists$weight)
      } else {
        out = easy.out
      }
      return(out)
    }
  } else if(param == "D"){
    session_to_use = which(dims$n.IDs > 0)[1]
    tem.mask = subset(data.mask, session == session_to_use)
    mask = as.matrix(unique(tem.mask[, c('x', 'y')]))
    tem.trap = subset(data.full, session == session_to_use)
    traps = tem.trap[, c('trap', 'trap_x', 'trap_y')]
    traps = traps[!duplicated(traps[['trap']]), ]
    traps = as.matrix(traps[order(traps$trap), c('trap_x', 'trap_y')])
    survey.length = survey.length[session_to_use]
    A = A[session_to_use]
    
    detpar.names = param.og[which(!(param.og %in% c('kappa', 'alpha', 'sigma.toa', 'D')))]
    pars <- sv[detpar.names]
    
    if(detfn == 'ss_het') ss_het = TRUE else ss_het = FALSE
    if(detfn == 'ss_dir') ss_dir = TRUE else ss_dir = FALSE
    
    if(ss_het | ss_dir | detfn == 'ss'){
      detfn = 'ss'
      if(ss.link == 'log'){
        detfn = 'log.ss'
      } else if(ss.link == 'spherical'){
        detfn = 'spherical.ss'
      }
      pars[["sigma.b0.ss"]] <- 0
      if(!any(detpar.names == 'b2.ss')) pars[['b2.ss']] = 0
    }
    
    if (!is.null(cutoff)){
      pars$cutoff <- cutoff
    }
    
    esa <- p.dot.defaultD(points = mask, traps = traps, detfn = detfn, ss_dir = ss_dir,
                         ss.link = ss.link, pars = pars, A = A, n.quadpoints = 8)
    
    return(log(dims$n.IDs[session_to_use]/(esa*survey.length)))
    #this part in the original function is difficult, better to ask for help
  } else if(param == "g0"){
    return(0.95)
  } else if(param == "lambda0"){
    return(2)
  } else if(param == "z"){
    return(1)
  } else if(param == "sigma.toa"){
    return(0.0025)
  } else if(param == "kappa"){
    return(10)
  } else if(param == "b0.ss"){
    if(ss.link == 'log'){
      return(log(max(data.full[['ss']], na.rm = TRUE)))
    } else {
      return(max(data.full[['ss']], na.rm = TRUE))
    }
  } else if(param == "sigma.b0.ss"){
    if(ss.link == 'log'){
      return(0.01 * log(max(data.full[['ss']], na.rm = TRUE)))
    } else {
      return(0.01 * max(data.full[['ss']], na.rm = TRUE))
    }
  } else if(param == "b1.ss"){
    max.ss = max(data.full[["ss"]], na.rm = TRUE)
    out = (max.ss - cutoff)/(buffer/2)
    if(ss.link == "log"){
      out = out/max.ss
    }
    return(out)
  } else if(param == "b2.ss"){
    return(0.1)
  } else if(param == "sigma.ss"){
    ss_valid = data.full[['ss']][data.full[['ss']] >= cutoff & !is.na(data.full[['ss']])]
    return(sd(ss_valid))
  } else if(param == "alpha"){
    return(2)
  } else if(param == "shape"){
    return(2)
  } else if(param == "scale"){
    if(detfn == "th"){
      out = default.sv(param = "sigma", info = info)
    } else if(detfn == "lth"){
      sigma = default.sv(param = "sigma", info = info)
      shape.1 = default.sv(param = "shape.1", info = info)
      out = (log(shape.1) - log(shape.1 - 0.5)) / sigma
    } else {
      stop("Detection function not recognised.")
    }
    return(out)
  } else if(param == "shape.1"){
    return(0.809017)
  } else if(param == "shape.2"){
    sigma = default.sv(param = "sigma", info = info)
    shape.1 = default.sv(param = "shape.1", info = info)
    scale = default.sv(param = "scale", info = info)
    return(shape.1 + scale * sigma)
  }
}


#default bounds
#the default bounds of "D" is not specified in the origin
default.bounds = function(param){
  if(param == 'g0') {
    return(c(0, 1))
  } else if(param %in% c('sigma', 'lambda0', 'shape.1', 'scale', 'b0.ss', 
                         'b1.ss', 'b2.ss', 'sigma.b0.ss', 'sigma.ss',
                         'z', 'sigma.toa', 'alpha')){
    return(c(0, 1e8))
  } else if(param == 'kappa'){
    return(c(0, 700))
  } else if(param %in% c('shape', 'shape.2')){
    return(c(-100, 100))
  } else if(param == "D"){
    return(c(0, 1e20))
  }
}

#link function
link.fun = function(link, value){
  if(link == "identity") return(value)
  if(link == "log"){
    if(value == 0) value = value + 1e-20
    
    return(log(value))
  }
  
  if(link == "logit"){
    if(value == 0) value = value + 1e-15
    if(value == 1) value = value - 1e-15
    
    return(log(value/(1 - value)))
  }
}

sort.data = function(dat, name){
  is.animal_ID = "animal_ID" %in% colnames(dat)
  
  if(name == 'data.full'){
    if(is.animal_ID){
      dat = dat[order(dat$session, dat$animal_ID, 
                      dat$ID, dat$trap), ]
    } else {
      dat = dat[order(dat$session, dat$ID, dat$trap), ]
    }
  }
  
  if(name == "data.dists.thetas"){
    dat = dat[order(dat$session, dat$trap, dat$mask), ]
  }
  
  if(name == "data.ID_mask"){
    if(is.animal_ID){
      dat = dat[order(dat$session, dat$animal_ID, dat$ID, dat$mask),]
    } else {
      dat = dat[order(dat$session, dat$ID, dat$mask), ]
    }
  }
  
  if(name == "data.mask"){
    dat = dat[order(dat$session, dat$mask),]
  }
  
  return(dat)
}



numeric_ID = function(dat){
  is.animal_ID = "animal_ID" %in% colnames(dat)
  if(is.animal_ID){
    tem.dat.non.na = subset(dat, !is.na(ID))
    tem.dat.non.na$animal_ID = as.numeric(as.factor(tem.dat.non.na$animal_ID))
    
    tem.dat.na = subset(dat, is.na(ID))
    u.id = unique(tem.dat.non.na$animal_ID)
    tem = vector('list', length(u.id))
    for(k in 1:length(u.id)){
      tem[[k]] = subset(tem.dat.non.na, animal_ID== u.id[k])
      tem[[k]]$ID = as.numeric(as.factor(tem[[k]]$ID))
    }
  } else {
    tem.dat.non.na = subset(dat, !is.na(ID))
    tem.dat.na = subset(dat, is.na(ID))
    u.session = unique(tem.dat.non.na$session)
    tem = vector('list', length(u.session))
    for(k in 1:length(u.session)){
      tem[[k]] = subset(tem.dat.non.na, session == u.session[[k]])
      tem[[k]]$ID = as.numeric(as.factor(tem[[k]]$ID))
    }
  }
  tem.dat.non.na = do.call('rbind', tem)
  dat = rbind(tem.dat.na, tem.dat.non.na)
  return(dat)
}

numeric_link = function(link){
  ans = numeric(length(link))
  for(i in 1:length(link)){
    if(link[i] == 'identity') ans[i] = 1
    if(link[i] == 'log') ans[i] = 2
    if(link[i] == 'logit') ans[i] = 3
    if(link[i] == 'spherical') ans[i] = 4
  }
  return(ans)
}

numeric_detfn = function(detfn){
  if(detfn == 'hn') return(1)
  if(detfn == 'hhn') return(2)
  if(detfn == 'hr') return(3)
  if(detfn == 'th') return(4)
  if(detfn == 'lth') return(5)
  if(detfn == 'ss') return(6)
  if(detfn == 'ss_dir') return(7)
  if(detfn == 'ss_het') return(8)
}

numeric_het_method = function(het_method){
  if(is.null(het_method)) {
    return(1)
  } else {
    if(het_method == "GH") return(2)
    if(het_method == "rect") return(3)
  }
}

cal_n_det = function(data.full){
  if("animal_ID" %in% colnames(data.full)){
    tem = aggregate(data.full$bincapt, list(session = data.full$session,
                                            animal_ID = data.full$animal_ID,
                                            ID = data.full$ID), sum)
    tem = tem[order(tem$session, tem$animal_ID, tem$ID),]
    return(tem$x)
  } else {
    tem = aggregate(data.full$bincapt, list(session = data.full$session,
                                            ID = data.full$ID), sum)
    tem = tem[order(tem$session, tem$ID),]
    return(tem$x)
  }
  
}