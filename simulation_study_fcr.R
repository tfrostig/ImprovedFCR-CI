# simulation study 
source('function_simulation.R')
library(tidyverse)

# cov mats 
ar_cov <- function(rho, p) {
  cov_mat <- rho^abs(outer(1:p, 1:p, FUN = "-"))
  return(cov_mat) 
}

block_cov <- function(rho, block_size, p) {
  cov_mat <-  matrix(0,ncol = p , nrow = p)
  for (i in 1:(p / block_size)) {
    start_pos <- 1 + (i - 1) * block_size
    cov_mat[start_pos:(start_pos + block_size -1), start_pos:(start_pos + block_size - 1)] <- rho 
  }
  diag(cov_mat) <- 1
  return(cov_mat)
}

# assume X is normalized 
fcr_adjustment <- function(x, p, lambda_upper, lambda_under, direction) {
  k <- length(x) 
  if (direction == 'positive') {
    upper <- x + qnorm(1 - lambda_upper) 
    lower <- x - qnorm(1 - (k / p) * lambda_under) 
  }
  if (direction == 'negative') {
    upper <- x + qnorm(1 - (k / p) * lambda_upper) 
    lower <- x - qnorm(1 - lambda_under) 
  }
  if (direction == 'two_sided') {
    upper <- lower <- rep(NA, k)
    upper[x < 0] <- x + qnorm(1 - (k / p) * lambda_upper) 
    upper[x > 0] <- x + qnorm(1 - lambda_upper)
    lower[x < 0] <- x - qnorm(1 - lambda_under) 
    lower[x > 0] <- x - qnorm(1 - (k / p) * lambda_under) 
    
  }
  
  return(cbind('lower' = lower, 'upper' = upper))
}



# multivariate normal
simulate_normal_distribution <- function(num_iter, theta, cov_mat, q, one_sided=TRUE) {
  p           <- nrow(cov_mat) 
  norm_data   <- MASS::mvrnorm(num_iter, theta, cov_mat)
  sign_mat    <- sign(norm_data)
  
  if (one_sided) {
    pval_mat    <- (1 - pnorm(norm_data))
    pval_mat    <- t(apply(pval_mat, 1, p.adjust, 'BH'))
    is_selected <- pval_mat < q 
    k_vec       <- apply(is_selected, 1, sum)
    # CI
    k_mat       <- replicate(p, k_vec)
    upper       <- norm_data - qnorm( (k_mat * (sign_mat == -1) / p + (sign_mat == 1)) * q / 2)
    lower       <- norm_data + qnorm( (k_mat * (sign_mat == 1) / p + (sign_mat == -1)) * q / 2)
  } else {
    pval_mat    <- 2 * (1 - pnorm(abs(norm_data)))
    pval_mat    <- t(apply(pval_mat, 1, p.adjust, 'BH'))
    is_selected <- pval_mat < q 
    k_vec       <- apply(is_selected, 1, sum)
    k_mat       <- replicate(p, k_vec)
    upper       <- norm_data - qnorm( (k_mat * (sign_mat == -1) / p) * q / 4 + (sign_mat == 1) * q / 2)
    lower       <- norm_data + qnorm( (k_mat * (sign_mat == 1) / p) * q / 4 + (sign_mat == -1) * q / 2)
  }
  
  # CI
  k_mat       <- replicate(p, k_vec)
  theta_mat   <- t(replicate(num_iter, theta))
  v_vec       <- apply((theta_mat > upper | theta_mat < lower) & is_selected, 1, sum)
  
  fdp <- v_vec / k_vec
  fdp[k_vec == 0] <- 0 
  
  return(c('FCR'    = mean(fdp), 
           'sd'     = sd(fdp), 
           'length' = mean(upper[k_vec > 0] - lower[k_vec > 0], na.rm = TRUE)))
}


# multivariate t
generate_multi_t <- function(num_iter, theta, p, cov_mat, df=5) {
  x <- replicate(df + 1, MASS::mvrnorm(num_iter, theta, cov_mat))
  s <- apply(x, 1:2, sd)
  x <- apply(x, 1:2, mean)
  
  #s <- matrix(rchisq(num_iter * p, df) / df, nrow = num_iter, ncol = p)
  #t_data <- t(apply(x / sqrt(s), 1, function(y) y + theta))
  return(list('x' = x, 's' = s))
}

simulate_multi_t <- function(num_iter, theta, cov_mat, df, q, one_sided=TRUE) {
  p           <- nrow(cov_mat)
  t_data_list <- generate_multi_t(num_iter, theta, p, cov_mat, df=df)
  sign_mat    <- sign(t_data_list$x)
  t_data      <- t_data_list$x / t_data_list$s
  
  if (one_sided) {
    pval_mat    <- (1 - pt(t_data, df))
    pval_mat    <- t(apply(pval_mat, 1, p.adjust, 'BH'))
    is_selected <- pval_mat < q 
    k_vec       <- apply(is_selected, 1, sum)
    k_mat       <- replicate(p, k_vec)
    upper       <- t_data_list$x - t_data_list$s * qt( (k_mat * (sign_mat == -1) / p + (sign_mat == 1)) * q / 2, df)
    lower       <- t_data_list$x + t_data_list$s * qt( (k_mat * (sign_mat == 1) / p + (sign_mat == -1)) * q / 2, df)
  } else {
    pval_mat    <- 2 * (1 - pt(abs(t_data), df))
    pval_mat    <- t(apply(pval_mat, 1, p.adjust, 'BH'))
    is_selected <- pval_mat < q 
    k_vec       <- apply(is_selected, 1, sum)
    k_mat       <- replicate(p, k_vec)
    upper       <- t_data_list$x - t_data_list$s * qt( (k_mat * (sign_mat == -1) / p) * q / 2 + (sign_mat == 1) * q / 4, df)
    lower       <- t_data_list$x + t_data_list$s * qt( (k_mat * (sign_mat == 1) / p) * q / 2 + (sign_mat == -1) * q / 4, df)
  }
  
  # results 
  theta_mat   <- t(replicate(num_iter, theta))
  v_vec       <- apply((theta_mat > upper | theta_mat < lower) & is_selected, 1, sum)
  
  fdp        <- v_vec / k_vec
  fdp[k_vec == 0] <- 0 
  return(c('FCR'    = mean(fdp), 
           'sd'     = sd(fdp),  
           'length' = mean(upper[k_vec > 0] - lower[k_vec > 0], na.rm = TRUE)))
}





# multivariate normal 
generate_mixture <- function(num_iter, theta, p, cov_mat) {
  mix_data <- MASS::mvrnorm(num_iter, rep(1, p), cov_mat) +
    MASS::mvrnorm(num_iter, rep(-1, p), cov_mat) +
    t(replicate(num_iter, theta))
}

simulate_mixture <- function(num_iter, theta, cov_mat, q) {
  p           <- nrow(cov_mat)
  mix_data    <- generate_mixture(num_iter, theta, p, cov_mat)
  pval_mat    <- 2 * (1 - EnvStats::pnormMix(abs(mix_data), mean1 = -1, mean2 = 1))
  pval_mat    <- matrix(pval_mat, ncol = p, nrow = num_iter)
  pval_mat    <- t(apply(pval_mat, 1, p.adjust, 'BH'))
  is_selected <- pval_mat < q 
  k_vec       <- apply(is_selected, 1, sum)
  sign_mat    <- sign(mix_data)
  # CI
  k_mat       <- replicate(p, k_vec)
  upper       <- mix_data - matrix(EnvStats::qnormMix( (k_mat * (sign_mat == -1) / p + (sign_mat == 1)) * q / 2,
                                                       mean1 = -1, mean2 = 1),
                                   ncol = p, nrow = num_iter)
  lower       <- mix_data + matrix(EnvStats::qnormMix( (k_mat * (sign_mat == 1) / p + (sign_mat == -1)) * q / 2,
                                                       mean1 = -1, mean2 = 1),
                                   ncol = p, nrow = num_iter)
  # results 
  theta_mat   <- t(replicate(num_iter, theta))
  v_vec       <- apply((theta_mat > upper | theta_mat < lower) & is_selected, 1, sum)
  
  fdp        <- v_vec / k_vec
  fdp[k_vec == 0] <- 0 
  return(c('FCR' = mean(fdp), 
           'sd' = sd(fdp), 
           'length' = mean(upper[k_vec > 0] - lower[k_vec > 0], na.rm = TRUE)))
}




generate_mu <- function(prop, alt_val, p) {
  alt_ind <- 
    theta   <- c(rep(-3, 10), rep(3, 10), rep(0, 80)) 
  
}


create_mus <- function(null_params = c(0, 0), discover_params = c(0, 0.5), prop_null = 0.7, n = 100) {
  null_num     <- round(n * prop_null)
  discover_num <- n - null_num 
  mus <- c(rnorm(null_num, null_params[1], null_params[2]), rnorm(discover_num, discover_params[1], discover_params[2]))
  return(sample(mus))
} 



# parameters 
q_sel    <- 0.05 
q        <- 0.05 

one_side = TRUE
prop  <- c(0.01, 0.05, seq(0.1, 1, 0.1))
set.seed(9999)

## ar matrix 
rhos       <- c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)
parameters <- expand.grid(alt_center = c(0, 3),
                          alt_sd     = c(0, 3, 5), 
                          prop_null  = c(0.99, 0.95, seq(0.9, 0.1, -0.1)),
                          rho        = rhos, 
                          one_side   = c(TRUE, FALSE),
                          p          = c(10, 50, 100),
                          dist       = c('normal', 't'))

parameters <- parameters %>% filter(!(alt_center == 0 & alt_sd == 0))

all_res <- list()

# general parameters 
n_sim_rep <- 1000 
rep_iter  <- 1

for (i in 1:nrow(parameters)) {
  cov_mat <- ar_cov(parameters$rho[i], parameters$p[i])
  tmp_res <- list() 
  for (j in 1:rep_iter) {
    theta   <- create_mus(null_params = c(0, 0), 
                          discover_params = c(parameters$alt_center[i], 
                                              parameters$alt_sd[i]), 
                          prop_null = parameters$prop_null[i],
                          n         = parameters$p[i])
    if (parameters$dist[i] == 'normal') {
      tmp_res[[j]] <- simulate_normal_distribution(n_sim_rep, theta, cov_mat, q, one_sided = parameters$one_side[i])
    }
    if (parameters$dist[i] == 't') {
      tmp_res[[j]] <- simulate_multi_t(n_sim_rep, theta, cov_mat, 5, q, one_sided = parameters$one_side[i])
    }
  }
  all_res[[i]] <- bind_rows(tmp_res) %>% colMeans()
}


all_result <- bind_cols(parameters, all_res %>% bind_rows())
all_result %>% filter(FCR > 0.05) %>% write.csv('parameters_above_threshold.csv')




# selected parameters 
sel_params <- all_result %>% filter(FCR > 0.05) %>% select(alt_center, alt_sd, prop_null, rho, one_side, p, dist)



all_res <- list()
# general parameters 
n_sim_rep <- 10000 


for (i in 1:nrow(sel_params)) {
  cov_mat <- ar_cov(sel_params$rho[i], sel_params$p[i])
  tmp_res <- list() 
  for (j in 1:rep_iter) {
    theta   <- create_mus(null_params = c(0, 0), 
                          discover_params = c(sel_params$alt_center[i], 
                                              sel_params$alt_sd[i]), 
                          prop_null = sel_params$prop_null[i],
                          n         = sel_params$p[i])
    if (sel_params$dist[i] == 'normal') {
      tmp_res[[j]] <- simulate_normal_distribution(n_sim_rep, theta, cov_mat, q, one_sided = sel_params$one_side[i])
    }
    if (sel_params$dist[i] == 't') {
      tmp_res[[j]] <- simulate_multi_t(n_sim_rep, theta, cov_mat, 5, q, one_sided = sel_params$one_side[i])
    }
  }
  all_res[[i]] <- bind_rows(tmp_res) %>% colMeans()
}


all_result <- bind_cols(sel_params, all_res %>% bind_rows())
all_result %>% filter(FCR > 0.05) %>% View()

all_result %>% filter(FCR > 0.05) %>% write.csv('parameters_above_threshold_after_selection.csv')


