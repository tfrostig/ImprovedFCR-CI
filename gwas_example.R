library(tidyverse)

source('WFB_cis.R')


dat <- read.csv('gwas_diabetes_example.csv')


dat_all_combined <- dat %>% 
  select(rs, chr, position, Region, contains('AllCombined')) %>%
  filter(complete.cases(.))

k  <- 10 
m  <- 393453 
a  <- 0.05 
ct <- qnorm(1 - 10^-5 / 2) 



dat_all_combined %>% 
  mutate_if(is.numeric, log) %>% 
  mutate(sd = (AllCombined_Upper - AllCombined_Lower) / (2 * qnorm(1 - a / 2)) ) %>% 
  mutate(AllCombined_BY_Lower = AllCombined_OR + sd * qnorm(k * a / (2 * m)),
         AllCombined_BY_Upper = AllCombined_OR - sd * qnorm(k * a / (2 * m))) %>% 
  mutate_if(is.numeric, exp) %>% 
  mutate(increase = AllCombined_BY_Upper / AllCombined_Upper) %>% 
  select(rs, contains('AllCombined'), increase) %>% 
  write.csv('results_gwas_diabetes.csv')


secondary <- dat_all_combined %>% 
  select(rs, AllCombined_OR, AllCombined_Lower, 
         AllCombined_Upper)

colnames(secondary) <- c('name', 'mean', 'lower_marginal', 'upper_marginal')


# find ME 
secondary_log <- secondary %>% 
  mutate_if(is.numeric, log) %>% 
  mutate(length = upper_marginal - lower_marginal,
         ME     = length / 2, 
         se     = ME / qnorm(0.975),
         stat   = mean / se, 
         is.selected = TRUE,
  ) 


shortest_ci <- secondary_log %>% 
  filter(is.selected) %>%
  select(stat) %>%
  unlist() %>% 
  sapply(., Shortest.CI, ct=ct, alpha=0.05) %>% 
  t()

qc_ci <- secondary_log %>% 
  filter(is.selected) %>%
  select(stat) %>%
  unlist() %>% 
  sapply(., QC.CI, ct=ct, alpha=0.05) %>% 
  t() %>%
  matrix(., ncol = 2)

secondary_log <- secondary_log %>% 
  mutate(lower_bf = mean - qnorm(1 - 0.025 * k / m) * se,
         upper_bf = mean + qnorm(0.975) * se) %>% 
  mutate(lower_by = mean - qnorm(1 - 0.025 * k / m) * se,
         upper_by = mean + qnorm(1 - 0.025 * k / m) * se) %>% 
  mutate(lower_bonferroni = mean - qnorm(1 - 0.025 / m) * se,
         upper_bonferroni = mean + qnorm(1 - 0.025 / m) * se) %>% 
  mutate(lower_shortest   = unlist(shortest_ci[ ,1]) * se,
         upper_shortest   = unlist(shortest_ci[ ,2]) * se) %>% 
  mutate(lower_qc         = unlist(qc_ci[ ,1]) * se, 
         upper_qc         = unlist(qc_ci[ ,2]) * se)

full_long <- bind_rows(secondary_log %>% 
                         pivot_longer(cols = contains('_'), 
                                      names_to = c('ci_type', 'bound'), 
                                      names_pattern = "(.*)_(.*)", values_to = "value", 
                                      values_transform = exp) %>% 
                         mutate(summary_type = 'Odds Ratio'),
                       secondary_log %>% 
                         pivot_longer(cols = contains('_'), 
                                      names_to = c('ci_type', 'bound'), 
                                      names_pattern = "(.*)_(.*)", values_to = "value") %>% 
                         mutate(summary_type = 'log(Odds Ratio)'))

# check difference in estimate of coef 
rel_names <- full_long %>% 
  filter(bound == 'by', ci_type == 'upper' & value < 1, is.selected) %>% 
  select(name) %>% 
  unlist() %>% 
  unique()

full_long %>% 
  filter(name %in% rel_names, ci_type == 'lower') %>%
  group_by(bound, summary_type) %>%
  summarise(min(value)) 



full_long %>% filter(name %in% rel_names, bound %in% c('by', 'marginal'))


full_wide <- full_long %>% 
  pivot_wider(names_from = ci_type, values_from = value) %>% 
  mutate(lower    = if_else(is.selected | bound == 'marginal', lower, NA_real_),
         upper    = if_else(is.selected | bound == 'marginal', upper, NA_real_),
         mean     = if_else(bound == 'marginal', mean, NA_real_)) %>% 
  mutate(crucial_val = 1 * (summary_type == 'Odds Ratio')) %>%
  mutate(mean = (summary_type == 'Odds Ratio') * exp(mean) +
           (summary_type != 'Odds Ratio') * mean)


full_wide %>% 
  filter(summary_type == 'log(Odds Ratio)')%>%
  ggplot(aes(x = name, 
             y = mean, color = bound), lwd = 0.8) +
  geom_hline(aes(yintercept = crucial_val)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    linetype = bound), 
                position = position_dodge2(w = 1),
                lwd = 0.8) +
  facet_grid(. ~ summary_type, scales = 'free') + 
  geom_point(position = position_dodge2(0.9)) +
  scale_linetype_manual(values = c(1 , 1, 1, 2, 1, 1),
                        labels = c('New', 'Simultaneous', 'BY', 'Marginal', 'WFB-QC', 'TN')) +
  scale_color_brewer(labels = c('New', 'Simultaneous', 'BY', 'Marginal', 'WFB-QC', 'TN'),
                     palette = 'Dark2') +
  theme_bw() + 
  theme(legend.position = 'bottom', legend.title = element_blank(),
        #axis.text.x = element_text( hjust = 1),
        legend.key.size = unit(2, "cm")) + 
  xlab('SNPs') + 
  ylab('') + 
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_flip()


full_long %>% pivot_wider(id_cols = c(name, bound, summary_type),
                          values_from = value, names_from = ci_type) %>% 
  group_by(summary_type, bound) %>%
  summarise(len = mean(upper - lower))


full_long %>% 
  filter(ci_type == 'lower', summary_type == 'Odds Ratio') %>% 
  group_by(bound) %>% 
  summarise(sum(value >= 1))


