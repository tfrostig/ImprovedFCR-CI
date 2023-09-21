library(tidyverse)

source('manuscript_scripts/WFB_cis.R')

secondary <- data.frame('name'  = c('Pelvic Infection', 'Additional Antibiotic Use', 'Additional Analgesia', 
                                    'Unplanned Ward Admission', 'Unplanned Consultation', 'Abdominal or Pelvic Pain', 
                                    'Presenece of Vaginal Bleeding', 'Delay of Usual Activities'),
                        'mean'  = c(0.6, 0.81, 0.72, 0.72, 0.6, 0.8, 1, 0.7),
                        'lower_marginal' = c(0.37, 0.65, 0.57, 0.37, 0.37, 0.68, 0.93, 0.27),
                        'upper_marginal' = c(0.96, 1.01, 0.92, 1.39, 0.97, 0.95, 1.07, 1.84))


# find ME 
secondary_log <- secondary %>% 
  mutate_if(is.numeric, log) %>% 
  mutate(length = upper_marginal - lower_marginal,
         ME     = length / 2, 
         se     = ME / qnorm(0.975),
         stat   = mean / se, 
         pval   = pnorm(stat),
         is.selected = pval < 0.025,
         rank   = rank(stat)
  ) 



shortest_ci <- secondary_log %>% 
  filter(is.selected) %>%
  select(stat) %>%
  unlist() %>% 
  sapply(., Shortest.CI, ct=1.96, alpha=0.05) %>% 
  t()


qc_ci <- secondary_log %>% 
  filter(is.selected) %>%
  select(stat) %>%
  unlist() %>% 
  sapply(., QC.CI, ct=1.96, alpha=0.05) %>% 
  t() %>%
  matrix(., ncol = 2)


secondary_log <- secondary_log %>% 
  mutate(lower_bf = mean - qnorm(0.975) * se,
         upper_bf = mean + qnorm(1 - 0.025 * sum(is.selected) / n()) * se) %>% 
  mutate(lower_by = mean - qnorm(1 - 0.025 * sum(is.selected) / n()) * se,
         upper_by = upper_bf) %>% 
  mutate(lower_bonferroni = mean - qnorm(1 - 0.025 / n()) * se,
         upper_bonferroni = mean + qnorm(1 - 0.025 / n()) * se) %>% 
  mutate(lower_shortest   = unlist(shortest_ci[ ,1]) * se,
         upper_shortest   = unlist(shortest_ci[ ,2]) * se) %>% 
  mutate(lower_qc         = unlist(qc_ci[ ,1]) * se, 
         upper_qc         = unlist(qc_ci[ ,2]) * se)

full_long <- bind_rows(secondary_log %>% 
                         pivot_longer(cols = contains('_'), 
                                      names_to = c('ci_type', 'bound'), 
                                      names_pattern = "(.*)_(.*)", values_to = "value", 
                                      values_transform = exp) %>% 
                         mutate(summary_type = 'Relative Risk'),
                       secondary_log %>% 
                         pivot_longer(cols = contains('_'), 
                                      names_to = c('ci_type', 'bound'), 
                                      names_pattern = "(.*)_(.*)", values_to = "value") %>% 
                         mutate(summary_type = 'log(Relative Risk)'))

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
  mutate(crucial_val = 1 * (summary_type == 'Relative Risk')) %>%
  mutate(mean = (summary_type == 'Relative Risk') * exp(mean) +
           (summary_type != 'Relative Risk') * mean)


full_wide %>% 
  arrange(rank) %>%
  filter(is.selected, summary_type == 'log(Relative Risk)')%>%
  #filter(type != "wfb_short") %>%
  ggplot(aes(x = factor(rank, levels = rank, labels = name, ordered = TRUE), 
             y = mean, color = bound), lwd = 0.8) +
  geom_hline(aes(yintercept = crucial_val)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    linetype = bound), 
                position = position_dodge2(w = 1),
                lwd = 0.8) +
  #facet_grid(. ~ summary_type, scales = 'free') + 
  geom_point(position = position_dodge2(0.9)) +
  scale_linetype_manual(values = c(1 , 1, 1, 2, 1, 1),
                        labels = c('New', 'Simultaneous', 'BY', 'Marginal', 'WFB-QC', 'TN')) +
  scale_color_brewer(labels = c('New', 'Simultaneous', 'BY', 'Marginal', 'WFB-QC', 'TN'),
                     palette = 'Dark2') +
  theme_bw() + 
  theme(legend.position = 'bottom', legend.title = element_blank(),
        #axis.text.x = element_text( hjust = 1),
        legend.key.size = unit(2, "cm")) + 
  guides(color = guide_legend(override.aes = list(shape = NULL))) +
  xlab('Secondary Outcomes') + 
  ylab('log(Relative Risk)') + 
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  coord_flip()

