set.seed(998)
# simulation of normal distribution 
library(tidyverse)

source('manuscript_scripts/polhyedral_cis.R')
source('manuscript_scripts/WFB_cis.R')


# Positive EFRON example 
one_sided_sim <- function(mu, q_sel, q) {
  x  <- rnorm(length(mu), mu, 1)
  ct <- qnorm(q_sel * rank(x) / m)
  
  pv           <- pnorm(x)
  selected_ind <- p.adjust(pv, method = 'BH') < q_sel
  
  k            <- sum(selected_ind)

  ci_df        <- data.frame('x'      = x, 
                             'mu'     = mu, 
                             'is_selected' = selected_ind, 
                             'threshold' = ct, 
                             'by__lower' = x + qnorm(q / 2 * k / m), 
                             'by__upper' = x - qnorm(q / 2 * k / m))
  
  
  ci_df['bf__lower'] <- x +( (x > 0) * qnorm(q / 2 * k / m) + (x < 0) * qnorm(q / 2) )
  ci_df['bf__upper'] <- x -( (x < 0) * qnorm(q / 2 * k / m) + (x > 0) * qnorm(q / 2) )
  
  
  
  reid_threshold    <- min(ci_df$x[!ci_df$is_selected])
  
  
  reid_cis <- list()
  #ci_df_filter <- ci_df[ci_df$is_selected, ]
  ci_df_filter  <- ci_df %>% filter(is_selected)
  short_res <- qc_res <- prat_res <- reid_res <- NULL 
  for (i in 1:nrow(ci_df_filter)) {
    # short_res <- rbind(short_res, 
    #                    wrap_shortest(ci_df_filter$x[i], ci_df_filter$threshold[i], q))
    qc_res    <- rbind(qc_res, 
                       wrap_qc(ci_df_filter$x[i], 0.8, -ci_df_filter$threshold[i], q))
    # prat_res  <- rbind(prat_res,
    #                    wrap_pratt(ci_df_filter$x[i], 1.2, ci_df_filter$threshold[i], q))
    reid_res  <- rbind(reid_res,
                       quantile_finder(y=ci_df_filter$x[i], 
                                       sigma=1, 
                                       p=q, 
                                       a=reid_threshold, 
                                       b=-reid_threshold))
  }
  
  
  ci_df_filter <- ci_df_filter %>%
    mutate(wfb_qc__lower = qc_res[ ,1], 
           wfb_qc__upper = qc_res[ ,2],
           reid__lower   = reid_res[ ,1],
           reid__upper   = reid_res[ ,2]) 
  
  plot_ci_df     <- ci_df_filter %>%
    pivot_longer(cols = contains('__'), 
                 names_to = c('type', 'bound'
                 ), 
                 names_pattern = "(.*)__(.*)" )
  
  
  plot_ci_df$value[!plot_ci_df$is_selected] <- NA 
  plot_ci_df <- plot_ci_df %>% filter(!(type %in% 'wfb_short')) %>% 
    mutate(type = recode(type, 
                         'bf' = 'New', 'by' = 'BY', 
                         'reid' = 'TN', 'wfb_qc' = 'WFB-QC'))
  
  
  
  plot_ci_df <- plot_ci_df %>% 
    mutate(type = factor(type, levels = c('New', 'BY', 'WFB-QC', 'TN'),
                         ordered=TRUE))
  
  return(plot_ci_df)
}

m  <- 10000 
r  <- 1000

q_sel  <- 0.1 
q      <- 0.1

mu <- c(rep(0, m - r), rnorm(r, -3, 1))

plot_ci_df <- one_sided_sim(mu, q_sel, q)


# creating plot 
colors <- RColorBrewer::brewer.pal(4, 'Dark2')
names(colors) <- c('New', 'BY', 'TN', 'WFB-QC')
linetypes     <- c('New' = 'solid', 'BY' = '64', 
                   'TN' = '56', 'WFB-QC' = '38')
pointtypes    <- c('New' = 3, 'BY' = 20,
                   'TN'  = 8, 'WFB-QC' = 6)


# Lines 
plot_ci_df %>%
  ggplot(aes(x = x, y = mu, size = is_selected))+
  geom_abline(slope=1, intercept=0, color = 'grey') +
  geom_point(alpha = 0.025) +
  geom_line(plot_ci_df %>%
              filter(bound == 'lower', is_selected),
            mapping = aes(color = type, x = x, y = value, linetype = type), size = 1.1) +
  geom_line(plot_ci_df %>%
              filter(bound == 'upper', is_selected),
            mapping = aes(color = type, x = x, y = value, linetype = type), size = 1.1) +
  ylab('True Effect Size') +
  xlab('Estimate') +
  theme_bw() +
  xlim(c(-7.5 ,0)) +
  ylim(c(-10, 2.5)) +
  scale_size_manual(values = c(0.001, 1.2), guide = "none") +
  scale_color_manual(values = colors,
                     name = '')  +
  scale_linetype_manual(values = linetypes, guide = "none") +
  theme(legend.position = 'bottom', text = element_text(size=14))


ggsave('fcr_paper/tex/images/efron_negative_example.pdf', height = 6, width = 10, dpi = 300)


result_list_one_sided <- list()
for (i in 1:1000) {
  result_list_one_sided[[i]] <- one_sided_sim(mu, q_sel, q) %>% mutate(iteration = i)
}

summary_sim <- bind_rows(result_list_one_sided) %>% 
  filter(is_selected) %>% 
  group_by(type, x, mu, iteration) %>%
  pivot_wider(id_cols = c(x, mu, type, iteration),
              names_from = bound,
              values_from = value) %>% 
  summarise(center = (upper + lower) / 2,
            mu     = mean(mu),
            mu_var = var(mu),
            length = upper - lower,
            is_over_est  = mu < lower,
            is_under_est = mu > upper) %>%
  group_by(type, iteration) %>% 
  summarise(dist_center = mean(abs(center - mu), na.rm = TRUE),
            length      = mean(length, na.rm = TRUE), 
            fdr_upper   = mean(is_over_est, na.rm = TRUE),
            fdr_under   = mean(is_under_est, na.rm = TRUE)) %>% 
  group_by(type) %>% 
  summarise_if(is.numeric, function(x) glue::glue('{round(mean(x), 3)} ({round(sd(x) / sqrt(1000), 4)})'))
  
write.csv(summary_sim, 'fcr_paper/Results/efron_negative.csv')






## two sided-comparision 



# Positive EFRON example 
two_sided_sim <- function(mu, q) {
  
  x  <- rnorm(length(mu), mu, 1)
  #ct <- qnorm(q_sel * 2 *  rank(x) / m)
  
  pv           <- 2 * (1 - pnorm(abs(x)))
  selected_ind <- abs(x) > 4 
  
  k            <- sum(selected_ind)
  sum(abs(x) > ct)
  
  ci_df        <- data.frame('x' = x, 
                             'mu' = mu, 
                             'is_selected' = selected_ind, 
                             'threshold' = 4, 
                             'by__lower' = x + qnorm(q / 2 * k / m), 
                             'by__upper' = x - qnorm(q / 2 * k / m))
  
  
  ci_df['bf__lower'] <- x +( (x > 0) * qnorm(q / 4 * k / m) + (x < 0) * qnorm(q / 2) )
  ci_df['bf__upper'] <- x -( (x < 0) * qnorm(q / 4 * k / m) + (x > 0) * qnorm(q / 2) )
  
  
  
  reid_threshold    <- 4
  
  
  reid_cis <- list()
  #ci_df_filter <- ci_df[ci_df$is_selected, ]
  ci_df_filter  <- ci_df %>% filter(is_selected)
  short_res <- qc_res <- prat_res <- reid_res <- NULL 
  for (i in 1:nrow(ci_df_filter)) {
    # short_res <- rbind(short_res, 
    #                    wrap_shortest(ci_df_filter$x[i], ci_df_filter$threshold[i], q))
    qc_res    <- rbind(qc_res, 
                       wrap_qc(ci_df_filter$x[i], 0.8, ci_df_filter$threshold[i], q))
    # prat_res  <- rbind(prat_res,
    #                    wrap_pratt(ci_df_filter$x[i], 1.2, ci_df_filter$threshold[i], q))
    reid_res  <- rbind(reid_res,
                       quantile_finder(y=ci_df_filter$x[i], 
                                       sigma=1, 
                                       p=q, 
                                       a=-reid_threshold, 
                                       b=reid_threshold))
  }
  
  
  ci_df_filter <- ci_df_filter %>%
    mutate(#wfb_short__lower = short_res[ ,1], 
           #wfb_short__upper = short_res[,2],
           wfb_qc__lower = qc_res[ ,1], 
           wfb_qc__upper = qc_res[ ,2],
           reid__lower   = reid_res[ ,1],
           reid__upper   = reid_res[ ,2]) 
  
  plot_ci_df     <- ci_df_filter %>%
    pivot_longer(cols = contains('__'), 
                 names_to = c('type', 'bound'), 
                 names_pattern = "(.*)__(.*)" )
  
  
  
  plot_ci_df$value[!plot_ci_df$is_selected] <- NA 
  plot_ci_df <- plot_ci_df %>% filter(!(type %in% 'wfb_short')) %>% 
    mutate(type = recode(type, 
                         'bf' = 'New', 'by' = 'BY', 
                         'reid' = 'TN', 'wfb_qc' = 'WFB-QC'))
  return(plot_ci_df)
}


m  <- 10000 
r  <- 1000

q      <- 0.1

mu <- c(rep(0, m - 2 * r), rnorm(r, 3, 1), rnorm(r, -3, 1))
plot_ci_df <- two_sided_sim(mu, q)

# creating plot 
colors <- RColorBrewer::brewer.pal(4, 'Dark2')
names(colors) <- c('New', 'BY', 'TN', 'WFB-QC')


plot_ci_df %>% 
  ggplot(aes(x = x, y = mu, size = is_selected))+ 
  geom_point(alpha = 0.025) + 
  geom_line(plot_ci_df %>%
              filter(bound == 'lower', is_selected, x > 0),
            mapping = aes(color = type, x = x, y = value, linetype = type), size = 1.2) + 

  geom_line(plot_ci_df %>%
              filter(bound == 'upper', is_selected, x > 0),
            mapping = aes(color = type, x = x, y = value, linetype = type), size = 1.2) + 

  
  geom_line(plot_ci_df %>%
              filter(bound == 'lower', is_selected, x < 0),
            mapping = aes(color = type, x = x, y = value, linetype = type), size = 1.2) +

  
  geom_line(plot_ci_df %>%
              filter(bound == 'upper', is_selected, x < 0),
            mapping = aes(color = type, x = x, y = value, linetype = type), size = 1.2) + 

  
  
  scale_linetype_manual(values = linetypes, guide = "none") +
  
  ylab('True Effect Size') + 
  xlab('Estimate') +  
  theme_bw() +
  #xlim(c(-7.5, 0)) +
  # ylim(c(-2.5, 10)) + 
  scale_size_manual(values = c(0.001, 1.2), guide = 'none') + 
  scale_color_manual(values = colors,
                     name = '')  +
  theme(legend.position = 'bottom', text = element_text(size=14))


ggsave('fcr_paper/tex/images/efron_two_side_example.pdf', height = 6, width = 10, dpi = 300)

result_list_two_sided <- list()
for (i in 1:1000) {
  result_list_two_sided[[i]] <- two_sided_sim(mu, q) %>% mutate(iteration = i)
}


summary_sim <- bind_rows(result_list_two_sided) %>% 
  filter(is_selected) %>% 
  group_by(type, x, mu, iteration) %>%
  pivot_wider(id_cols = c(x, mu, type, iteration),
              names_from = bound,
              values_from = value) %>% 
  summarise(center = (upper + lower) / 2,
            mu     = mean(mu),
            mu_var = var(mu),
            length = upper - lower,
            is_over_est  = mu < lower,
            is_under_est = mu > upper) %>%
  group_by(type, iteration) %>% 
  summarise(dist_center = mean(abs(center - mu), na.rm = TRUE),
            length      = mean(length, na.rm = TRUE), 
            fdr_upper   = mean(is_over_est, na.rm = TRUE),
            fdr_under   = mean(is_under_est, na.rm = TRUE)) %>% 
  group_by(type) %>% 
  summarise_if(is.numeric, function(x) glue::glue('{round(mean(x), 3)} ({round(sd(x) / sqrt(5000), 4)})'))


write.csv(summary_sim, 'fcr_paper/Results/efron_two_side.csv')
