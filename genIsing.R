genIsing <- function(nNode, propPositive = 0.5, rewire = 0, inclusion = NULL, mean = 0.3, sd = 0.1, ...){
  
  if (is.null(inclusion)){
    
    if (rewire > 1) {
      stop("The 'rewire' argument needs to be specified as a probability between 0 and 1.")
    }
    
    # Generate graph structure:
    graph <- ifelse(bootnet::genGGM(nNode, propPositive = propPositive, p = rewire, ...) != 0, 1, 0) 
  
  } else {

    if (rewire > 0) {
      warning("The 'rewire' argument has no effect when 'inclusion' is specified. Ignoring 'rewire'.")
    } else if(inclusion > 1){
      stop("The 'inclusion' argument needs to be specified as a probability between 0 and 1.")
    }
    # Generate graph structure:
    p <- inclusion
    edges <- sample(c(0, 1), nNode * (nNode-1)/2, TRUE, prob = c(1-p, p))
    
    graph <- matrix(0, nNode, nNode)
    graph[lower.tri(graph)] <- edges
    graph <- graph + t(graph)
    
  }
  
  # Add weights to graph structure:
  graph[lower.tri(graph)] <- graph[lower.tri(graph)] * rnorm(sum(lower.tri(graph)), mean, sd)
  graph[upper.tri(graph)] <- t(graph)[upper.tri(graph)]
  
  return(graph)
}

