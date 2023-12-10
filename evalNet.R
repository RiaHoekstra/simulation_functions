evalNet <- function(true, est, metric = NULL, directed = FALSE){
  
  # Warning if input is not a matrix: 
  if (!is.matrix(true) | !is.matrix(est)) stop("Input must be a weight matrix")
  
  # Check if network is directed:
  if (directed){
    true <- c(true)
    est <- c(est)
  } else {
    true <- true[upper.tri(true,diag=FALSE)]
    est <- est[upper.tri(est,diag=FALSE)]
  }        
  
  if (is.null(metric)){
    metric <- c("sensitivity", "signed_sensitivity", "sensitivity_top50",
                "sensitivity_top25", "sensitivity_top10", "specificity", 
                "precision", "precision_top50", "precision_top25", 
                "precision_top10", "correlation", "correlation_abs")
  } else if(any(!(metric %in% c("sensitivity", "signed_sensitivity", "sensitivity_top50",
                      "sensitivity_top25", "sensitivity_top10", "specificity", 
                      "precision", "precision_top50", "precision_top25", 
                      "precision_top10", "correlation", "correlation_abs")))) {
    stop("metric invalid; needs to be 'sensitivity', 'signed_sensitivity', 'sensitivity_top50',
                      'sensitivity_top25', 'sensitivity_top10', 'specificity', 
                      'precision', 'precision_top50', 'precision_top25', 
                      'precision_top10', 'correlation', 'correlation_abs'")
  }
  
  # Output list:
  out <- list()
  
  # True positives:
  TruePos <- sum(est != 0 &  true != 0)
  
  # False pos:
  FalsePos <- sum(est != 0 & true == 0)
  
  # True Neg:
  TrueNeg <- sum(est == 0 & true == 0)
  
  # False Neg:
  FalseNeg <- sum(est == 0 & true != 0)
  
  # Sensitivity:
  if("sensitivity" %in% metric){
    out$sensitivity <- TruePos / (TruePos + FalseNeg)
  }
  
  # Signed sensitivity:
  if("signed_sensitivity" %in% metric){
    TruePos_signed <- sum(est != 0 &  true != 0 & sign(est) == sign(true))
    out$signed_sensitivity <- TruePos_signed / (TruePos + FalseNeg)
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
    out$specificity <- TrueNeg / (TrueNeg + FalsePos)
  }
  
  # Precision (1 - FDR):
  if("precision" %in% metric){
    out$precision <- TruePos / (FalsePos + TruePos)
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

  # Correlation:
  if("correlation" %in% metric){
    out$correlation <- cor(est, true)
  }
  
  # Correlation between absolute edges:
  if("correlation_abs" %in% metric){
  out$abs_cor <- cor(abs(est),abs(true))
  }
  
  # Nonzero edge weight correlation: 
  
  return(out)
}
