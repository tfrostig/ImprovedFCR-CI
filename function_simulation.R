
# theta (numeric) - parameters of interest
# thres (numeric)- threshold of selection
# quant_func (function) - quantile function of the distriubtion function (dist_func)
# dist_func (function) - distribution function 
# alpha - level of significance to maintain 
# alt (string) - alternative 
simulation_threshold_fcr <- function(theta, thres, quant_func, dist_func, 
                                     alpha = 0.05, alt = 'two.sided', err = 'FCR',
                                     ...) {
  m       <- length(theta)
  # sampling
  x       <- dist_func(n = m, ...)
  x       <- x + theta
  
  # direction 
  ind         <- find_selected(x, thres, alt)
  k           <- sum(ind)
  ## FCR 
  if (err == 'FCR') {
    ci <- find_confidence(x > 0, 
                          out_lambda = alpha, 
                          in_lambda  = alpha * k / m, 
                          quant_func = quant_func,
                          ...)
  }
  ## SOS 
  if(err == 'SOS') { 
    ci <- find_confidence(x > 0, 
                          out_lambda = alpha / k, 
                          in_lambda  = alpha / m,
                          quant_func = quant_func,
                          ...)
  }
  # arranging data
  ci_mat      <- cbind('x'       = x, 
                       'lower'   = x - ci$lower_bound,
                       'upper'   = x + ci$upper_bound,
                       'sign'    = sign(x), 
                       'theta'   = theta)
  selected_ci <- ci_mat[ind, , drop = FALSE] # ensure matrix 
  if (k == 0) { 
    # non were selected 
    return(c('inside' = 0, 'outside' = 0, 'k' = 0))
  } else { 
    # at least one was selected
    upper_select <- selected_ci[ ,'sign'] == 1
    lower_select <- selected_ci[ ,'sign'] == -1
    v_outside <- sum(selected_ci[upper_select, 'theta'] > selected_ci[upper_select, 'upper']) +
      sum(selected_ci[lower_select, 'theta'] < selected_ci[lower_select, 'lower'])
    
    v_inside <- sum(selected_ci[upper_select, 'theta'] < selected_ci[upper_select, 'lower']) +
      sum(selected_ci[lower_select, 'theta'] > selected_ci[lower_select, 'upper'])
    # overall coverage 
    return(c('inside' = v_inside, 'outside' = v_outside, 'k' = nrow(selected_ci)))
  }
}

# find indices of selection 
find_selected  <- function(x, thres, alt) {
  if (alt == 'two.sided') {
    ind <- abs(x) > thres
  }
  if (alt == 'right') {
    ind <- x > thres
  }
  if (alt == 'left') {
    ind <- x < thres
  }
  return(ind)
}


find_confidence <- function(is_pos, out_lambda, in_lambda, quant_func, ...) {
  cik         <- -quant_func(p = 0.5 * out_lambda, ...)
  cim         <-  quant_func(p  = 1 - 0.5 * in_lambda, ...)
  
  dik         <-  quant_func(p  = 1 - 0.5 * out_lambda, ...)
  dim         <- -quant_func(p = 0.5 * in_lambda, ...)
  ci          <- data.frame('upper_bound' = cik * is_pos + dim * !is_pos,
                            'lower_bound' = cim * is_pos + dik * !is_pos)

  return(ci)
}


