# Author:   Matthew Wiens, Metrum Research Group
# Purpose:  This script performs replicates of the 
#           logistic regression analysis in tte-demo.Rmd to 
#           confirm that the results are not due to random
#           variability


library(tidyverse)
library(mrgsolve)
library(survival)
library(survminer)

set.seed(123)

theme_set(theme_bw())

# patient size
n <- 2000

# Simulation replicates
replicates <- 1000

mod <- mrgsolve::mread(here::here("model", "nonmem", "simmod", "pk2cmt.cpp"))

data <- bbr::nm_join(bbr::read_model(here::here("model/nonmem/106")))

dose_rec <- filter(data, EVID==1)

pars <- distinct(dose_rec,ID,CL,V2,Q,V3,KA,AMT,RF,ACTARM)

dose <- tibble( AMT = 100, ID = 1, TIME=0:5 * 21, EVID=1, CMT=1)   %>%
  inner_join(pars %>% slice_head(n = 1) %>% select(CL, V2, Q, V3, KA), 
             by = character()) %>%
  mutate(CL = CL*4)



run_replicate <- function(i) {

  
  indiv_dose <- dose %>%
    select(-ID) %>% 
    inner_join(tibble(ID = 1:n,
                      ETA1 = rnorm(n, 0, 0.3)),
               by = character()) %>%
    mutate(CL = CL * exp(ETA1)) %>%
    arrange(ID, TIME)
  
  out2 <- mrgsim_df(mod, 
                    indiv_dose,
                    recover = "ACTARM,RF,CL,DOSE",  
                    carry.out = "EVID",
                    recsort = 3, 
                    tgrid = 0:(6 * 21)) %>% 
    as_tibble() %>% 
    distinct() %>%
    filter(EVID != 1) %>%
    filter(TIME > 0) %>%
    group_by(ID) %>%
    mutate(CAVG = AUC/TIME) %>%
    ungroup()
  
  CAVGC1_df = out2 %>% 
    filter(TIME == 21) %>% 
    select(ID, CAVGC1 = CAVG) 
  
  exposures = out2 %>%
    inner_join(CAVGC1_df, by = "ID") %>%
    select(ID, TIME, CAVG, CAVGC1) %>% 
    arrange(ID, TIME)
  
  last_time <- max(exposures$TIME)
  
  
  TTEs <- tibble(ID = 1:n) %>%
    mutate(can_have_event = rbernoulli(n, 0.75)) %>%
    mutate(actual_event_time = 1 + ceiling(rweibull(n, shape = 0.45, scale = 100))) %>%
    mutate(EV = can_have_event * as.numeric(actual_event_time <= last_time) ,
           EVTIME = if_else(can_have_event,
                            pmin(actual_event_time, last_time), 
                            last_time))
  analysis_dat <- exposures %>%
    inner_join(TTEs, by = c("ID" = "ID", "TIME" = "EVTIME")) %>%
    mutate(CAVGTTE_quartile = as.factor(glue::glue("Q{ntile(CAVG, 4)}")),
           CAVGC1_quartile = as.factor(glue::glue("Q{ntile(CAVGC1, 4)}")))
  
  
  cavg_results <- glm(EV ~ CAVG,
      family = "binomial",
    data = analysis_dat) %>%
    broom::tidy(conf.int = T) %>%
    filter(term == "CAVG") 
      
  c1_results <- glm(EV ~ CAVGC1,
                      family = "binomial",
                      data = analysis_dat) %>%
    broom::tidy(conf.int = T) %>%
    filter(term == "CAVGC1") 
  
  bind_rows(cavg_results, c1_results)

}  


simulation_results <- map_dfr(1:1000, run_replicate)

# Number of times 0 is in the confidence interval

simulation_results %>%
  group_by(term) %>%
  summarise(conf_in_contain_0 = mean((conf.low <= 0) & (conf.high >= 0)))

# Create figure

simulation_figure <- simulation_results %>%
  select(-std.error, -statistic,  -p.value) %>%
  pivot_longer(cols = c(estimate, conf.low, conf.high), names_to = "Statistic", values_to = "value") %>%
  mutate(Statistic = factor(Statistic, levels = c("estimate", "conf.low", "conf.high")),
         term = factor(term)) %>%
  mutate(Statistic = forcats::fct_recode(Statistic,
                                         "Point Estimate\n\n" = "estimate", 
                                         "Lower Bound of\n95% CI\n" = "conf.low", 
                                         "Upper Bound of\n95% CI\n" = "conf.high"),
         term = forcats::fct_recode(term, 
                                    "plain('Average Concentration to Event (Cavg')[TE]~plain(')')" = "CAVG",
                                    "plain('Average Concentration in Cycle 1 (Cavg')[C1]~plain(')')" = "CAVGC1")) %>%
  ggplot(aes(x = value, group = Statistic, linetype = Statistic)) +
  geom_vline(size = 2, linetype = 4, xintercept = 0) +
  geom_density() +
  facet_wrap(~term, ncol = 1, labeller = label_parsed) +
  labs(y = "", x = "Coefficent of Exposure in the Logistic Regression Analysis") +
  theme(axis.text.y = element_blank())

ggsave(filename = "supp_fig_1.png", 
       plot = simulation_figure,
       path = here::here("deliv", "figure", "publication"),
       device = "png",
       width = 178,
       height = 150,
       units = "mm",
       dpi = 400)
