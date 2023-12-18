evalNet <- function(true, est, metric = NULL, directed = FALSE){
  
  # Warning if input is not a matrix: 
  if (!is.matrix(true) | !is.matrix(est)) stop("Input must be a weight matrix")
  
  if (is.null(metric)){
    metric <- c("sensitivity", "signed_sensitivity", "sensitivity_top50",
                "sensitivity_top25", "sensitivity_top10", "specificity", 
                "precision", "precision_top50", "precision_top25", 
                "precision_top10", "jaccard_index", "correlation", 
                "correlation_abs", "correlation_true", "bias", "bias_true", 
                "centrality_Pearson_cor", "centrality_Kendall_cor", 
                "centrality_top5", "centrality_top3", "centrality_top1")
  } else if(any(!(metric %in% c("sensitivity", "signed_sensitivity", "sensitivity_top50",
                                "sensitivity_top25", "sensitivity_top10", "specificity", 
                                "precision", "precision_top50", "precision_top25", 
                                "precision_top10", "jaccard_index",
                                "correlation", "correlation_abs", "correlation_true", 
                                "bias", "bias_true", "centrality_Pearson_cor", 
                                "centrality_Kendall_cor", "centrality_top5", 
                                "centrality_top3", "centrality_top1")))) {
    stop("metric invalid; needs to be 'sensitivity', 'signed_sensitivity', 'sensitivity_top50', 
    'sensitivity_top25', 'sensitivity_top10', 'specificity', 'precision', 'precision_top50', 
    'precision_top25', 'precision_top10', 'jaccard_index', 
         'correlation', 'correlation_abs', 'correlation_true', 'bias', 'bias_true', 
         'centrality_Pearson_cor', 'centrality_Kendall_cor', 'centrality_top5', 
         'centrality_top3', 'centrality_top1'")
  }
  
  # Check if centrality needs to be computed: 
  if(any(metric %in% c("centrality_Pearson_cor", "centrality_Kendall_cor", "centrality_top5", "centrality_top3", "centrality_top1"))){
    centTrue <- qgraph::centrality(true)
    centEst <- qgraph::centrality(est)
  }
  
  # Check if network is directed:
  if (directed){
    true <- c(true)
    est <- c(est)
  } else {
    true <- true[upper.tri(true, diag = FALSE)]
    est <- est[upper.tri(est, diag = FALSE)]
  }        
  
  # Output list:
  out <- list()
  
  # True positives:
  truePos <- sum(est != 0 &  true != 0)
  
  # False pos:
  falsePos <- sum(est != 0 & true == 0)
  
  # True Neg:
  trueNeg <- sum(est == 0 & true == 0)
  
  # False Neg:
  falseNeg <- sum(est == 0 & true != 0)
  
  # Sensitivity:
  if("sensitivity" %in% metric){
    out$sensitivity <- truePos / (truePos + falseNeg)
  }
  
  # Signed sensitivity:
  if("signed_sensitivity" %in% metric){
    truePosSigned <- sum(est != 0 &  true != 0 & sign(est) == sign(true))
    out$signed_sensitivity <- truePosSigned / (truePos + falseNeg)
  }
  
  # Sensitivity top 50%:
  if("sensitivity_top50" %in% metric){
    top50 <- which(abs(true) > median(abs(true[true!=0])))
    out[["sensitivity_top50"]] <- sum(est[top50]!=0 & true[top50] != 0) / sum(true[top50] != 0)
  }
  
  # Sensitivity top 25%:
  if("sensitivity_top25" %in% metric){
    top25 <- which(abs(true) > quantile(abs(true[true!=0]), 0.75))
    out[["sensitivity_top25"]] <- sum(est[top25]!=0 & true[top25] != 0) / sum(true[top25] != 0)
  }

  # Sensitivity top 10%:
  if("sensitivity_top10" %in% metric){
    top10 <- which(abs(true) > quantile(abs(true[true!=0]), 0.90))
    out[["sensitivity_top10"]] <- sum(est[top10]!=0 & true[top10] != 0) / sum(true[top10] != 0)
  }

  # Specificity:
  if("specificity" %in% metric){
    out$specificity <- trueNeg / (trueNeg + falsePos)
  }
  
  # Precision (1 - FDR):
  if("precision" %in% metric){
    out$precision <- truePos / (falsePos + truePos)
  }
  
  # precision top 50% (of estimated edges):
  if("precision_top50" %in% metric){
    top50 <- which(abs(est) > median(abs(est[est!=0])))
    out[["precision_top50"]] <- sum(est[top50]!=0 & true[top50] != 0) / sum(est[top50] != 0)
  }
  
  # precision top 25%:
  if("precision_top25" %in% metric){
    top25 <- which(abs(est) > quantile(abs(est[est!=0]), 0.75))
    out[["precision_top25"]] <- sum(est[top25]!=0 & true[top25] != 0) / sum(est[top25] != 0)
  }
  
  # precision top 10%:
  if("precision_top10" %in% metric){
    top10 <- which(abs(est) > quantile(abs(est[est!=0]), 0.90))
    out[["precision_top10"]] <- sum(est[top10]!=0 & true[top10] != 0) / sum(est[top10] != 0)
  }
  
  # Jaccard index:
  if("jaccard_index" %in% metric){
    inter <- sum(true & est)
    union <- sum(true | est)
    
    if (union == 0) { # avoid division by zero
      out$jaccard_index <- 1
    } else {
      # Calculate Jaccard Index
      out$jaccard_index <- inter / union
    }
  }

  # correlation:
  if("correlation" %in% metric){
    out$correlation <- cor(est, true)
  }
  
  # correlation between absolute edges:
  if("correlation_abs" %in% metric){
  out$abs_correlation <- cor(abs(est),abs(true))
  }
  
  # correlation between true edge weights:
  if("correlation_true" %in% metric){
    if (truePos > 0){
      trueEdges <- est != 0 & true != 0
      out$correlation_true <- cor(est[trueEdges], true[trueEdges])
    } else {
      out$correlation_true <- NA
    }
  }
  
  # average bias between edge weights:
  if("bias" %in% metric){
    out$bias <- mean(abs(est-true), na.rm = TRUE)
  }
  
  # average bias between true edge weights:
  if("bias_true" %in% metric){
    if (truePos > 0){
      trueEdges <- est != 0 & true != 0
      out$bias_true <- mean(abs(est[trueEdges] - true[trueEdges]))
    } else {
      out$bias_true <- NA
    }
  }
  
  # Pearson correlations:
  if("centrality_Pearson_cor" %in% metric){
    if(directed){
      out$in_strength_correlation <- cor(centTrue$InDegree, centEst$InDegree)
      out$out_strength_correlation <- cor(centTrue$OutDegree, centEst$OutDegree)
    } else {
      out$strength_correlation <- cor(centTrue$OutDegree, centEst$OutDegree)
    }
    out$closeness_correlation <- cor(centTrue$Closeness, centEst$Closeness)
    out$betweenness_correlation <- cor(centTrue$Betweenness, centEst$Betweenness)
  }
  
  # Kendall correlations:
  if("centrality_Kendall_cor" %in% metric){
    if(directed){
      out$in_strength_correlation_kendall <- cor(centTrue$Indegree, centEst$Indegree, method = "kendall")
      out$out_strength_correlation_kendall <- cor(centTrue$OutDegree, centEst$OutDegree, method = "kendall")
    } else{
      out$strength_correlation_kendall <- cor(centTrue$OutDegree, centEst$OutDegree, method = "kendall")
    }
    out$closeness_correlation_kendall <- cor(centTrue$Closeness, centEst$Closeness, method = "kendall")
    out$betweenness_correlation_kendall <- cor(centTrue$Betweenness, centEst$Betweenness, method = "kendall")
    
  }
  
  # Centrality top 5:
  if("centrality_top5" %in% metric){
    if(directed){
      out$in_strength_top5 <- mean(order(centEst$InDegree, decreasing = TRUE)[1:5] %in% order(centTrue$InDegree, decreasing = TRUE)[1:5])
      out$out_strength_top5 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1:5] %in% order(centTrue$OutDegree, decreasing = TRUE)[1:5])
    } else{
      out$strength_top5 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1:5] %in% order(centTrue$OutDegree, decreasing = TRUE)[1:5])
    }
    out$closeness_top5 <-  mean(order(centEst$Closeness, decreasing = TRUE)[1:5] %in% order(centTrue$Closeness, decreasing = TRUE)[1:5])
    out$betweenness_top5 <-  mean(order(centEst$Betweenness, decreasing = TRUE)[1:5] %in% order(centTrue$Betweenness, decreasing = TRUE)[1:5])
  }
  
  # Centrality top 3:
  if("centrality_top3" %in% metric){
    if(directed){
      out$in_strength_top3 <- mean(order(centEst$InDegree, decreasing = TRUE)[1:3] %in% order(centTrue$InDegree, decreasing = TRUE)[1:3])
      out$out_strength_top3 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1:3] %in% order(centTrue$OutDegree, decreasing = TRUE)[1:3])
    } else{
      out$strength_top3 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1:3] %in% order(centTrue$OutDegree, decreasing = TRUE)[1:3])
    }
    out$closeness_top3 <-  mean(order(centEst$Closeness, decreasing = TRUE)[1:3] %in% order(centTrue$Closeness, decreasing = TRUE)[1:3])
    out$betweenness_top3 <-  mean(order(centEst$Betweenness, decreasing = TRUE)[1:3] %in% order(centTrue$Betweenness, decreasing = TRUE)[1:3])
    
  }
  
  # Centrality top 1:
  if("centrality_top1" %in% metric){
    if(directed){
      out$in_strength_top1 <- mean(order(centEst$InDegree, decreasing = TRUE)[1] %in% order(centTrue$InDegree, decreasing = TRUE)[1])
      out$out_strength_top1 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1] %in% order(centTrue$OutDegree, decreasing = TRUE)[1])
    }else{
      out$strength_top1 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1] %in% order(centTrue$OutDegree, decreasing = TRUE)[1])
    }
    out$closeness_top1 <- mean(order(centEst$Closeness, decreasing = TRUE)[1] %in% order(centTrue$Closeness, decreasing = TRUE)[1])
    out$betweenness_top1 <- mean(order(centEst$Betweenness, decreasing = TRUE)[1] %in% order(centTrue$Betweenness, decreasing = TRUE)[1])
  }
  
  # Return evaluation metrics:
  return(out)
}
