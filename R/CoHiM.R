#' Run CoHiM for null hypotheses testing under dependence
#' 
#' @param p_matrix a n by J matrix of p-values from n studies for J SNPs
#' @param q FDR nominal
#' 
#' @return List containing
#' \describe{
#'   \item{rLIS_mat}{Matrix of pairwise replicability scores}
#'   \item{pi_mat}{Matrix of estimated stationary probabilities}
#'   \item{A_list}{List of transition matrices}
#'   \item{f1_mat}{Matrix of estimated alternative densities (study A)}
#'   \item{f2_mat}{Matrix of estimated alternative densities (study B)}
#'   \item{e_value}{Vector of e-values}
#'   \item{q_star}{Optimal FDR threshold level}
#'   \item{rej_index}{Logical vector of discovery decisions}
#' }
#' @export

CoHiM <- function(p_matrix, q = 0.05){
  n = nrow(p_matrix); J = ncol(p_matrix)
  if(n < 2){
    stop("At least two lists of p-values are required!")
  }
  n_pair = choose(n, 2)
  rLIS_mat = matrix(NA, n_pair, J)
  pi_mat = matrix(NA, n_pair, 4)
  A_list = list()
  f1_mat = matrix(NA, n_pair, J)
  f2_mat = matrix(NA, n_pair, J)
  null_prop_list = numeric(n_pair)
  i_pair = 0
  for(i1 in 1:(n-1)){
    for(i2 in (i1+1):n){
      i_pair = i_pair + 1
      res.hmm = RepLIS(p_matrix[i1, ], p_matrix[i2, ])
      rLIS_mat[i_pair, ] = res.hmm$repLIS
      pi_mat[i_pair, ] = res.hmm$pi
      A_list[[i_pair]] = res.hmm$A
      f1_mat[i_pair, ] = res.hmm$f1
      f2_mat[i_pair, ] = res.hmm$f2
      null_prop_list[i_pair] = sum(res.hmm$pi[1:3])
    }# end for i2
  } # end for i1
  min_eBH_result = min_eBH(rLIS_mat, null_prop_list)
  result = list(rLIS_mat = rLIS_mat, 
                pi_mat = pi_mat, A_list =A_list,
                f1_mat = f1_mat, f2_mat = f2_mat,
                e_value = min_eBH_result$e_value,
                q_star = min_eBH_result$q_star, 
                rej_index = min_eBH_result$rej_index)
  return(result)
} # end CoHiM