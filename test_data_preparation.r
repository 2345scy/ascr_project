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
  }
}



