###################################################################################################################
#sort out output for the function
out = vector('list', 36)

names(out) = c("fn", "coefficients", "coeflist", "se", "loglik", "maxgrad", "cor", "vcov", "npar", "npar_re",
               "npar_sdrpt", "npar_rep", "npar_total", "hes", "eratio", "D.mask", "mm.ihd", "args", "n.sessions",
               "fit.types", "infotypes", "detpars", "suppars", "D.betapars", "phases", "par.links", "par.unlinks",
               "fit.ihd", "re.detfn", "fit.freqs", "first.calls", "scale.covs", "model.formula", "fgam", "all.covariates", 
               "output.tmb")

#create output for TMB model
out[['output.tmb']] = vector('list', 7)
names(out[['output.tmb']]) = c('coef_link', 'se_link', 'esa', 'DX', 'param.og', 'param.extend', 'param.fix', 'param.info.table')

#give an index to each parameter, to make it easier to find it in "data.par"
par.id = 1:nrow(data.par)
names(par.id) = data.par$par

param.with.esa = c(param.og, 'esa')
param.og.4cpp = c(param.og.4cpp, 'esa')
o_value = vector('list', length(param.with.esa))
names(o_value) = param.with.esa
o_sd = vector('list', length(param.with.esa))
names(o_sd) = param.with.esa
for(i in 1:length(param.with.esa)){
  o_value[[param.with.esa[i]]] = o$value[which(names(o$value) == param.og.4cpp[i])]
  o_sd[[param.with.esa[i]]] = o$sd[which(names(o$value) == param.og.4cpp[i])]
}


############################################################################################################
#the 1st component: "fn"
out$fn = 'secr_tmb'


#############################################################################################################
#the 2nd component: "coefficients", there are "esa.[session]" not included yet
#this one can be generated from "name.fited.pars", "name.fixed.par", "name.extend.par"

#create a vector for parameter names used in various of output components
#adjust it to be the same with the output of current ADMB model

#here "name.fited.pars" is used instead of "param.with.esa" is to make sure the order
#of parameters are the same as TMB model, and "esa" is not in "param.og"
#if nothing goes wrong, 'name.fited.pars" shoule be exactly the same with "param.with.esa"
name.fited.pars = unique(names(o$value))
name_output = vector('list', length(name.fited.pars))
names(name_output) = name.fited.pars

#names for tmb output, "esa" would be totally separated from other parameters
name_output_tmb = name_output
name_output_tmb$esa = NULL

name_split_full = strsplit(colnames(data.full), " _ ")
name_split_mask = strsplit(colnames(data.mask), " _ ")

#out_coef is a temporary object, which will be assigned to out$coefficients later
out_coef = vector('list', length(name.fited.pars))
names(out_coef) = name.fited.pars
#take 'esa' out from tmb output
coef_link_tmb = out_coef
coef_link_tmb$esa = NULL

for(i in name.fited.pars){
  #esa is not in the original parameter list, deal with it separately
  if(i != 'esa'){
    #firstly, sort out names for output coefficients
    index_full = which(sapply(name_split_full, function(x) x[1] == i))
    index_mask = which(sapply(name_split_mask, function(x) x[1] == i))
    name_output[[i]] = c(colnames(data.full)[index_full], colnames(data.mask)[index_mask])
    name_output_tmb[[i]] = name_output[[i]]
    name_output[[i]] = gsub(" _ ", ".", name_output[[i]])
    #if i is "D", then keep the .(Intercept) no matter what model it is
    if(i != "D" & !i %in% name.extend.par){
      name_output[[i]] = gsub("\\.\\(Intercept\\)", "", name_output[[i]])
    }
    
    if(!i %in% name.fixed.par){
      #then, input the estimated values
      out_coef[[i]] = o_value[[i]]
      names(out_coef[[i]]) = paste(name_output[[i]], 'link', sep = "_")
      coef_link_tmb[[i]] = o_value[[i]]
      names(coef_link_tmb[[i]]) = name_output_tmb[[i]]
      
      if(i == "D"){
        #for "D", take it as extended no matter whether it actually is
        tem = out_coef[[i]]
      } else {
        if(!i %in% name.extend.par){
          tem = unlink.fun(data.par[par.id[i], 'link'], out_coef[[i]])
        } else {
          tem = out_coef[[i]]
        }
      }
      names(tem) = name_output[[i]]
      out_coef[[i]] = c(out_coef[[i]], tem)
      
      
      #for homo "D" models, there is a element in the vector named "D"
      if(i == 'D' & !i %in% name.extend.par){
        tem = unlink.fun(data.par[par.id[i], 'link'], tem)
        names(tem) = 'D'
        out_coef[[i]] = c(out_coef[[i]], tem)
      }
    } else {
      coef_link_tmb[[i]] = fix.input[[i]]
      names(coef_link_tmb[[i]]) = name_output_tmb[[i]]
    }
  } else {
    name_output[[i]] = paste('esa', 1:dims$n.sessions, sep = ".")
    out_coef[[i]] = o_value[[i]]
    names(out_coef[[i]]) = name_output[[i]]
  }
}

out$coefficients = do.call('c', out_coef)
names(out$coefficients) = do.call('c', lapply(out_coef, function(x) names(x)))
out[['output.tmb']][['coef_link']] = coef_link_tmb
#######################################################################################################
#the 3rd component: "coeflist"

out$coeflist = vector('list', 0)

D_pars = name_output$D
det_pars = name.fited.pars[which(!name.fited.pars %in% c('D', 'kappa', 'alpha', 'sigma.toa', 'esa'))]
sup_pars = name.fited.pars[which(name.fited.pars %in% c('kappa', 'alpha', 'sigma.toa'))]

tem = val(o_value[['D']])
for(i in 1:length(name_output$D)){
  tem_name = paste0('D_betapars_link[', i, ']')
  out$coeflist[[tem_name]] = tem[i]
}

for(i in 1:length(det_pars)){
  if(!det_pars[i] %in% name.fixed.par){
    tem = val(o_value[[det_pars[i]]])
    tem_name = paste0('detpars_link[', i, ']')
    out$coeflist[[tem_name]] = tem
  }
}

if(length(sup_pars) > 0){
  for(i in 1:length(sup_pars)){
    if(!sup_pars[i] %in% name.fixed.par){
      tem = val(o_value[[sup_pars[i]]])
      tem_name = paste0('suppars_link[', i, ']')
      out$coeflist[[tem_name]] = tem
    }
  }
}


for(i in 1:dims$n.sessions){
  tem = val(o_value[['esa']][i])
  tem_name = paste0('esa[', i, ']')
  out$coeflist[[tem_name]] = tem
}

#deal with "D" and "D.mask[s]"
if(!'D' %in% name.extend.par){
  out$coeflist[['D']] = unlink.fun(data.par[par.id['D'], 'link'], val(o_value[['D']]))
  for(i in 1:dims$n.sessions){
    tem_name = paste0('D_mask[', i, ']')
    out$coeflist[[tem_name]] = rep(out$coeflist[['D']], dims$n.masks[i])
  }
} else {
  tem_n_col = data.par[par.id['D'], c('n_col_full', 'n_col_mask')]
  tem_full = DX.full[['D']] %*% val(o_value[['D']])[1:tem_n_col[1, 1]]
  if(tem_n_col[1, 2] > 0){
    tem_mask = DX.mask[['D']] %*% val(o_value[['D']])[(tem_n_col[1, 1] + 1):sum(tem_n_col)]
  } else {
    tem_mask = rep(0, sum(dims$n.masks))
  }
  
  for(i in 1:dims$n.sessions){
    tem_name = paste0('D_mask[', i, ']')
    out$coeflist[[tem_name]] = tem_full[which(data.full$session ==i)[1]] + 
      tem_mask[which(data.mask$session == i)]
    out$coeflist[[tem_name]] = unlink.fun(data.par[par.id['D'], 'link'],
                                          out$coeflist[[tem_name]])
  }
}



####################################################################################################
#the 4th component: "se"

out_se = vector('list', length(name.fited.pars))
names(out_se) = name.fited.pars

#in tmb output, not include 'esa'
se_tmb = out_se
se_tmb$esa = NULL

for(i in name.fited.pars){
  if(i != 'esa'){
    if(!i %in% name.fixed.par){
      #then, input the estimated values
      out_se[[i]] = o_sd[[i]]
      names(out_se[[i]]) = paste(name_output[[i]], 'link', sep = "_")
      
      se_tmb[[i]] = o_sd[[i]]
      names(se_tmb[[i]]) = name_output_tmb[[i]]
      if(i == "D"){
        #for "D", take it as extended no matter whether it actually is
        tem = out_se[[i]]
      } else {
        if(!i %in% name.extend.par){
          tem = delta.fun(data.par[par.id[i], 'link'], out_se[[i]], o_value[[i]])
        } else {
          tem = out_se[[i]]
        }
      }
      names(tem) = name_output[[i]]
      out_se[[i]] = c(out_se[[i]], tem)
      
      
      #for homo "D" models, there is a element in the vector named "D"
      if(i == 'D' & !i %in% name.extend.par){
        tem = delta.fun(data.par[par.id[i], 'link'], tem, o_value[['D']])
        names(tem) = 'D'
        out_se[[i]] = c(out_se[[i]], tem)
      }
    } else {
      se_tmb[[i]] = NA
      names(se_tmb[[i]]) = name_output_tmb[[i]]
    }
  } else {
    out_se[[i]] = o_sd[[i]]
    names(out_se[[i]]) = name_output[[i]]
  }
  
}

out$se = do.call('c', out_se)
names(out$se) = do.call('c', lapply(out_se, function(x) names(x)))
out[['output.tmb']][['se_link']] = se_tmb
#####################################################################################################
#insert the necessary TMB output
tem = vector('list', 2)
names(tem) = c('DX_non_mask', 'DX_mask')
tem[['DX_non_mask']] = DX.full
tem[['DX_mask']] = DX.mask

out[['output.tmb']][['DX']] = tem

tem = vector('list', 2)
names(tem) = c('value', 'std')
tem[['value']] = val(o_value$esa)
tem[['std']] = o_sd$esa

out[['esa']] = tem

out[['output.tmb']][['param.og']] = param.og
out[['output.tmb']][['param.extend']] = name.extend.par
out[['output.tmb']][['param.fix']] = name.fixed.par
out[['output.tmb']][['param.info.table']] = data.par


######################################################################################################
#the 5th component: "loglik"
out$loglik = -1 * opt$objective

#####################################################################################################
#the 6th component: "maxgrad"
out$maxgrad = max(abs(o$gradient.fixed))

#####################################################################################################
#the 7th and 8th component: "cor" and "vcov", after communication it is fine to
#use the output from TMB directly.
tem = name_output_tmb

for(i in names(tem)){
  if(i %in% name.fixed.par) tem[i] = list(NULL)
}
tem = do.call('c', tem)
out$vcov = o$cov.fixed
dimnames(out$vcov) = list(tem, tem)

#convert the covariance matrix to correlation matrix
tem = sqrt(diag(out$vcov))
out$cor = out$vcov / outer(tem, tem)


######################################################################################################
#The 9th component: "npar"

out$npar = sum(sapply(name_output, length))

######################################################################################################
#From the 10th component to the 14th component, not sure what they are, skip them temporarily




######################################################################################################
#The 15th component: 'eratio'
#obtain all estimated values excluding 'esa'
tem = o$value[which(names(o$value) != 'esa')]
out$eratio = max(abs(o$gradient.fixed) / abs(tem))


######################################################################################################
#The 16th component: "D.mask"
out$D.mask = vector('list', dims$n.sessions)

for(s in 1:dims$n.sessions){
  tem = paste0("D_mask[", s, "]")
  out$D.mask[[s]] = out$coeflist[[tem]]
  tem = stringr::str_pad(as.character(1:dims$n.masks[s]), nchar(as.character(dims$n.masks[s])), 
                         pad = '0', side = 'left')
  names(out$D.mask[[s]]) = paste(paste0('D_mask[', s, ']'), tem, sep = ".")
}



########################################################################################################
#The 17th component: "mm.ihd"

out$mm.ihd = vector('list', dims$n.sessions)

for(s in 1:dims$n.sessions){
  #firstly got the first index of this session in data.full
  tem = min(which(data.full$session == s))
  #then, got the data from DX.full, which denotes "design matrix for data.full"
  tem = DX.full$D[tem,]
  #finally, generate a matrix with the same length of n.masks and contains the information
  #of "D" under "session - level"
  tem = matrix(rep(tem, dims$n.masks[s]), nrow = dims$n.masks[s], byrow = T)
  
  if(data.par[par.id['D'], 'n_col_mask'] > 0){
    out$mm.ihd[[s]] = cbind(tem, DX.mask$D[data.mask$session == s,])
  } else {
    out$mm.ihd[[s]] = tem
  }
  colnames(out$mm.ihd[[s]]) = gsub("D\\.", "", name_output$D)
  rownames(out$mm.ihd[[s]]) = NULL
}


######################################################################################################
#The 18th component: "args"

out$args = arg.input

######################################################################################################
#the 19th component: 'n.sessions'
out$n.sessions = dims$n.sessions
#####################################################################################################
#the 20th component: 'fit.types'
tem = logical(5)
names(tem) = c('bearing', 'dist', 'ss', 'toa', 'mrds')

tem['bearing'] <- ('bearing' %in% bucket_info)
tem['dist'] <- ('dist' %in% bucket_info)
tem['ss'] <- ('ss' %in% bucket_info)
tem['toa'] <- ('toa' %in% bucket_info)
tem['mrds'] <- ('mrds' %in% bucket_info)

out$fit.types = tem

#####################################################################################################
#the 21st component: 'infotypes'

out$infotypes = bucket_info[which(bucket_info %in% c('bearing', 'dist', 'ss', 'toa', 'mrds'))]


#####################################################################################################
#The 22nd component: 'detpars'
out$detpars = param.og[which(!param.og %in% c('kappa', 'alpha', 'sigma.toa', 'D'))]


######################################################################################################
#the 23rd component: 'suppars'
out$suppars = param.og[which(param.og %in% c('kappa', 'alpha', 'sigma.toa'))]
if(length(out$suppars) == 0) out$suppars = NULL

######################################################################################################
#the 24th component: 'D.betapars'
out$D.betapars = name_output$D

######################################################################################################
#the 25th component: 'phases'
all.pars = c(out$detpars, out$suppars, out$D.betapars)

out$phases = ifelse(all.pars %in% names(out$args$fix), -1, 1)
names(out$phases) = all.pars
out$phases = as.list(out$phases)

#####################################################################################################
#the 26th component: 'par.links'
#the 27th component: 'par.unlinks'

out$par.links = vector('list', length(all.pars))
names(out$par.links) = all.pars

out$par.unlinks = vector('list', length(all.pars))
names(out$par.unlinks) = all.pars

for(i in c(out$detpars, out$suppars)){
  if(!i %in% name.extend.par){
    if(data.par[par.id[i], 'link'] == 'log'){
      out$par.links[[i]] = cus_log
      out$par.unlinks[[i]] = cus_log_unlink
    } else if(data.par[par.id[i], 'link'] == 'logit'){
      out$par.links[[i]] = cus_logit
      out$par.unlinks[[i]] = cus_logit_unlink
    } else {
      out$par.links[[i]] = cus_identity
      out$par.unlinks[[i]] = cus_identity
    }
  } else {
    out$par.link[[i]] = cus_identity
    out$par.unlinks[[i]] = cus_identity
  }
}

for(i in out$D.betapars){
  out$par.links[[i]] = cus_identity
  out$par.unlinks[[i]] = cus_identity
}


######################################################################################################
#the 28th component: 'fit.ihd'
#only if the model is 'D' extended and the extension is "marks" extension only,
#in other words, the model which is equivalent to the incumbent "inhomo" model,
#then this 'ihd' is set to TRUE

if('D' %in% name.extend.par & data.par[par.id['D'], 'n_col_full'] == 1){
  out$fit.ihd = TRUE
} else {
  out$fit.ihd = FALSE
}


#####################################################################################################
#the 29th component: 're.detfn', not sure what this is
#in all of the tested sample data, this is FALSE, so set it as FALSE temporarily

out$re.detfn = FALSE


#####################################################################################################
#the 30th component: 'fit.freqs'
#this part is not included in the TMB template yet

out$fit.freqs = !is.null(cue.rates)

######################################################################################################
#the 31st component: 'first.calls'
#based on the instruction, set it to FALSE

out$first.calls = FALSE


#######################################################################################################
#the 32nd component: 'scale.covs'
#from the code, it seems this is a function receive data.frame as input and output a data.frame as well
#not 100% sure

out$scale.covs = scale.covs

#######################################################################################################
#the 33rd component: 'model.formula'
if('D' %in% name.extend.par){
  out$model.formula = as.formula(paste(c('gam.resp', as.character(par.extend$model$D)), collapse = ''))
} else {
  out$model.formula = as.formula('gam.resp~1')
}

########################################################################################################
#the 34th component: 'fgam'
out$fgam = fgam

########################################################################################################
#the 35th component: 'all.covariates'
if(is(out$args$mask, 'list')){
  out$all.covariates = do.call('rbind', out$args$mask)
} else {
  out$all.covariates = out$args$mask
}

if(is.scale){
  for(i in 1:ncol(out$all.covariates)){
    tem = out$all.covariates[, i]
    out$all.covariates[, i] = (tem - mean(tem)) / sd(tem)
  }
}

if(ncol(fgam$X) > 1){
  out$all.covariates = cbind(out$all.covariates, fgam$X[,-1])
}
out$all.covariates = as.data.frame(out$all.covariates)
