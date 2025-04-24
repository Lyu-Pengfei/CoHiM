#' Calculatet e-values and conduct eBH procedure for discoveries
#' 
#' @param rLIS_mat rLIS matrix with n_pair rows and m columns 
#' @param null_prop_list null proportion estimation for each pair
#' @param q FDR nominal
#' 
#' @return list of pairwise FDR level q_*, e-values, rejection number and decisions
#' @export

min_eBH <- function(rLIS_mat, null_prop_list, q = 0.05){
  n_pair = nrow(rLIS_mat)
  m = ncol(rLIS_mat)
  
  if(n_pair == 1){
    rLIS_sort = sort(rLIS_mat)
    rLIS_cummean = cummean(rLIS_sort)
    rej_num = sum(rLIS_cummean <= q)
    rej_threshold = rLIS_sort[rej_num]
    e_value = m*(rLIS_mat <= rej_threshold) / rLIS_cummean[rej_num]
    result = list(q_star = q, rej_num = rej_num, e_value = e_value,
                  rej_index = (e_value > 0))
    return(result)
  } else {# end n_pair
    # Generate all possible q' candidates from cumulative averages of sorted Lfdrs
    all_q_candidates <- c()
    
    sorted_lfdr <- matrix(NA, nrow = n_pair, ncol = m)
    avg_lfdr <- matrix(NA, nrow = n_pair, ncol = m)
    for (i_pair in 1:n_pair) {
      sorted_lfdr[i_pair, ] <- sort(rLIS_mat[i_pair, ])
      cumsum_lfdr <- cumsum(sorted_lfdr[i_pair, ])
      avg_lfdr[i_pair, ] <- cumsum_lfdr / (1:m)
      all_q_candidates <- c(all_q_candidates, avg_lfdr[i_pair, ])
    } # end for i_pair
    # Include q in the candidates
    all_q_candidates <- c(all_q_candidates, q)
    
    # Filter valid candidates (0 < q' <= q) and remove duplicates
    valid_candidates <- unique(all_q_candidates[all_q_candidates > 0 & all_q_candidates <= q])
    
    # Sort valid candidates in descending order
    valid_candidates <- sort(valid_candidates, decreasing = TRUE)
    # Initialize optimal q_star
    q_star <- 0
    
    # Iterate over each candidate q'
    j_start = 1; j_end = length(valid_candidates); j_chosen = 1
    while (j_start < j_end-1 | j_chosen < j_end) {
      j_chosen = ceiling((j_start + j_end)/2)
      # print(paste0("j_start = ", j_start, ", j_end = ", j_end, ", j_chosen = ", j_chosen))
      qc = valid_candidates[j_chosen]
      e_i_pair <- matrix(0, nrow = n_pair, ncol = m)
      
      for (i_pair in 1:n_pair) {
        R_i <- max(which(avg_lfdr[i_pair, ] <= qc), default = 0)
        
        if (R_i > 0) {
          threshold_i <- sorted_lfdr[i_pair, R_i]
          rejected_j <- rLIS_mat[i_pair, ] <= threshold_i
          sum_reject <- sum(rLIS_mat[i_pair, rejected_j])
          e_i_pair[i_pair, rejected_j] <- m / sum_reject
        } # end for R_i
        e_i_pair = e_i_pair * null_prop_list[i_pair]
      } # end for i_pair 
      
      # Compute e_j as the minimum e_ij across studies for each hypothesis
      e_value <- apply(e_i_pair, 2, min)
      
      # Apply e-BH procedure
      sorted_e <- sort(e_value, decreasing = TRUE)
      R <- which(sorted_e * (1:m) >= m / q)
      
      if (length(R) > 0) {
        j_end = j_chosen
      } else{
        j_start = j_chosen
      }# end if
      # print(R)
    } # end while
    if (length(R) == 0) {
      result <- list(q_star = q_star, rej_num = 0, e_value = e_value, rej_index = rep(FALSE, m))
      return(result)
    } else {
      rej_num = max(R)
      rej_threshold = sorted_e[rej_num]
      rej_index = (e_value >= rej_threshold)
      result = list(q_star = qc, rej_num = max(R), e_value = e_value, rej_index = rej_index)
      return(result)
    }
  }
} # end min_eBH