#' Run RepLIS for rLIS calculation in pairwise replicability analysis
#'
#' @param pa Vector of p-values from study A
#' @param pb Vector of p-values from study B
#' @param oracle Logical. If TRUE, use oracle inputs; otherwise run EM.
#' @param pi Initial vector for pi (oracle mode)
#' @param A Initial matrix for A (oracle mode)
#' @param f1 Vector of f1 values (oracle mode)
#' @param f2 Vector of f2 values (oracle mode)
#'
#' @return A list with repLIS, fdr, and HMM parameters
#' @export

RepLIS <- function(pa, pb, oracle = FALSE, pi, A, f1, f2){
  
  if(oracle == FALSE){
    pvals.cutoff = 1e-15
    pa[pa == 0] <- min(min(pa[pa != 0]), pvals.cutoff)
    pb[pb == 0] <- min(min(pb[pb != 0]), pvals.cutoff)
    
    pi0_pa <- min(qvalue::pi0est(pa)$pi0, 0.999)
    pi0_pb <- min(qvalue::pi0est(pb)$pi0, 0.999)
    
    res <- em_hmm(pa, pb, pi0_pa, pi0_pb)
  } else{
    rLIS = replis(pi, A, f1, f2)
    FDR = fdr(rLIS)
    res <- list(replis = rLIS, fdr = FDR)
  }
  
  return(res)
}