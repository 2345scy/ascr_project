test_data = function(model_type){
  context("Testing model fits")
  if(model_type == "simple fitting -- half normal"){
    #recover the example data to original format
    captures <<- data.frame(session = rep(1, 127 * 6), ID = rep(1:127, 6), occasion = rep(1, 127 * 6),
                          trap = rep(1:6, each = 127), stringsAsFactors = FALSE)
    captures <<- captures[as.logical(as.vector(example.data$capt[["bincapt"]])),]
    
    traps <<- example.data$traps
    masks <<- example.data$mask
    
    par.extend <<- NULL
    ss.opts <<- NULL
    survey.length <<- NULL
    cue.rates <<- NULL
    detfn <<- NULL
    bounds <<- NULL
    sv_list <<- NULL
    fix_list <<- NULL
    local_input <<- FALSE
    sound.speed <<- 331
  } else if (model_type == "joint bearing/dist fitting"){
    #detfn = "hn"
    captures = data.frame(session = rep(1, 127 * 6), ID = rep(1:127, 6), occasion = rep(1, 127 * 6),
                            trap = rep(1:6, each = 127), stringsAsFactors = FALSE)
    captures$bearing = as.vector(example.data$capt[["bearing"]])
    captures$dist = as.vector(example.data$capt[['dist']])
    captures <<- captures[as.logical(as.vector(example.data$capt[["bincapt"]])),]
    traps <<- example.data$traps
    masks <<- example.data$mask
    par.extend <<- NULL
    ss.opts <<- NULL
    survey.length <<- NULL
    cue.rates <<- NULL
    detfn <<- NULL
    bounds <<- NULL
    sv_list <<- NULL
    fix_list <<- list(g0 = 1)
    local_input <<- FALSE
    sound.speed <<- 331
  } else if (model_type == "simple fitting -- hazard halfnormal"){
    captures <<- data.frame(session = rep(1, 127 * 6), ID = rep(1:127, 6), occasion = rep(1, 127 * 6),
                            trap = rep(1:6, each = 127), stringsAsFactors = FALSE)
    captures <<- captures[as.logical(as.vector(example.data$capt[["bincapt"]])),]
    
    traps <<- example.data$traps
    masks <<- example.data$mask
    
    par.extend <<- NULL
    ss.opts <<- NULL
    survey.length <<- NULL
    cue.rates <<- NULL
    detfn <<- "hhn"
    bounds <<- NULL
    sv_list <<- NULL
    fix_list <<- NULL
    local_input <<- FALSE
    sound.speed <<- 331
  } else if (model_type == "simple fitting -- hazard rate"){
    captures <<- data.frame(session = rep(1, 127 * 6), ID = rep(1:127, 6), occasion = rep(1, 127 * 6),
                            trap = rep(1:6, each = 127), stringsAsFactors = FALSE)
    captures <<- captures[as.logical(as.vector(example.data$capt[["bincapt"]])),]
    
    traps <<- example.data$traps
    masks <<- example.data$mask
    
    par.extend <<- NULL
    ss.opts <<- NULL
    survey.length <<- NULL
    cue.rates <<- NULL
    detfn <<- "hr"
    bounds <<- NULL
    sv_list <<- list(z = 5)
    fix_list <<- NULL
    local_input <<- FALSE
    sound.speed <<- 331
  } else if (model_type == "bearing fitting"){
    #detfn = "hn"
    captures = data.frame(session = rep(1, 127 * 6), ID = rep(1:127, 6), occasion = rep(1, 127 * 6),
                          trap = rep(1:6, each = 127), stringsAsFactors = FALSE)
    captures$bearing = as.vector(example.data$capt[["bearing"]])
    captures <<- captures[as.logical(as.vector(example.data$capt[["bincapt"]])),]
    traps <<- example.data$traps
    masks <<- example.data$mask
    par.extend <<- NULL
    ss.opts <<- NULL
    survey.length <<- NULL
    cue.rates <<- NULL
    detfn <<- NULL
    bounds <<- NULL
    sv_list <<- NULL
    fix_list <<- list(g0 = 1)
    local_input <<- FALSE
    sound.speed <<- 331
  } else if (model_type == "dist fitting"){
    #detfn = "hn"
    captures = data.frame(session = rep(1, 127 * 6), ID = rep(1:127, 6), occasion = rep(1, 127 * 6),
                          trap = rep(1:6, each = 127), stringsAsFactors = FALSE)
    captures$dist = as.vector(example.data$capt[['dist']])
    captures <<- captures[as.logical(as.vector(example.data$capt[["bincapt"]])),]
    traps <<- example.data$traps
    masks <<- example.data$mask
    par.extend <<- NULL
    ss.opts <<- NULL
    survey.length <<- NULL
    cue.rates <<- NULL
    detfn <<- NULL
    bounds <<- NULL
    sv_list <<- NULL
    fix_list <<- list(g0 = 1)
    local_input <<- FALSE
    sound.speed <<- 331
  } else if (model_type == "toa fitting"){
    #detfn = "hn"
    captures = data.frame(session = rep(1, 127 * 6), ID = rep(1:127, 6), occasion = rep(1, 127 * 6),
                          trap = rep(1:6, each = 127), stringsAsFactors = FALSE)
    captures$toa = as.vector(example.data$capt[["toa"]])
    captures <<- captures[as.logical(as.vector(example.data$capt[["bincapt"]])),]
    traps <<- example.data$traps
    masks <<- example.data$mask
    par.extend <<- NULL
    ss.opts <<- NULL
    survey.length <<- NULL
    cue.rates <<- NULL
    detfn <<- NULL
    bounds <<- NULL
    sv_list <<- NULL
    fix_list <<- list(g0 = 1)
    local_input <<- FALSE
    sound.speed <<- 331
  }

}



