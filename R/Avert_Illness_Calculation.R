####################################################################################################
## WHO: assessment of the influenza vaccine impact
##
## Program to estimate the confidence intervals for estimates of influenza prevented by vaccination
##
####################################################################################################


## load packages ##
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, openxlsx, tidyverse)

## load functions ##
source("HelperFunction.R")
source("AbsenceVaxFunction.R")
source("PresenceVaxFunction.R")

#program starts
start.time <- Sys.time()



##### Simulation Prep #####

# Setting the number of simulations
nsim <- 5000

# data with inputs 
data_import <- read_excel("Phase 2 data input tool_final.xlsm",  
                          sheet = "Export") %>% 
  filter(!is.na(mnth_hosp)) %>% 
  mutate(month = row_number()) %>% 
  select(-c(deployment, seasonality))

# scenario
scenario <- unique(data_import$scenario)

# length of time period (in months)
time_period <- nrow(data_import)

# Replicate the dataset nsim times
temp_sim <- do.call("rbind", replicate(nsim, data_import, simplify = FALSE)) %>% 
  mutate(sim_index = rep(1:nsim, each = time_period),
         year = (month - 1) %/% 12)

# Ratio parameters
set.seed(123)
ratio_param <- temp_sim %>% 
  group_by(sim_index, year, 
           mult_nonhosp, mult_ma, mult_deathhosp, death_overall) %>%  # just to keep multipliers 
  dplyr::summarize(
    total_hosp = sum(mnth_hosp),
    .groups = "drop"
  ) %>% 
  mutate(
    sim_total_hosp = rpois(n(), lambda = total_hosp),
    sim_total_nohosp = rpois(n(), lambda = total_hosp * mult_nonhosp),
    sim_total_macase = rpois(n(), lambda = (total_hosp + total_hosp * mult_nonhosp) * mult_ma),
    sim_total_death = ifelse(!is.na(death_overall),
                             rpois(n(), lambda = total_hosp * mult_deathhosp * 1 / death_overall),
                             rpois(n(), lambda = total_hosp * mult_deathhosp) )
  ) %>% 
  mutate(
    sim_hnhratio = case_when(
      sim_total_hosp < 1  ~ sim_total_nohosp,
      sim_total_hosp >= 1 ~ sim_total_nohosp / sim_total_hosp
    ),
    sim_maratio = sim_total_macase / (sim_total_hosp + sim_total_nohosp),
    sim_dhratio = sim_total_death / sim_total_hosp
  )

# Vaccine coverage parameter
set.seed(123)
vaccine_param <- temp_sim %>% 
  # by each season, total vaccine coverage 
  group_by(sim_index, year, 
           target_pop_size) %>% 
  dplyr::summarize(
    total_coverage = sum(mnth_coverage),
    .groups = "drop"
  ) %>% 
  mutate(coverage_se = total_coverage / target_pop_size ) %>%
  mutate(senorm_cov = rnorm(n(), mean = 0, sd = .75),
         sim_cov_se = coverage_se * senorm_cov,
         sim_cov    = total_coverage + sim_cov_se,
         adj_sim_vc = case_when(sim_cov >= 0 ~ sim_cov,
                                sim_cov < 0 ~ 0)) 

# VE sampling
# assumption: VE hosp must be larger than VE ill and VE must not be negative
set.seed(123)
adj_sim_ve_ill <- c()
adj_sim_ve_hosp <- c()

for (i in 0:max(temp_sim$year)) {
  
  # Initialize temporary vectors for this parameter set
  temp_adj_sim_ve_ill <- c()
  temp_adj_sim_ve_hosp <- c()
  
  # filter to corresponding season data for input parameter
  temp_season <- temp_sim %>% 
    filter(year == i)
  
  if(all(temp_season$adjve_ill == temp_season$adjve_hosp)){
    
    # Start sampling
    # VE against illness estimates
    ve_illness_batch <- beta_pert(unique(temp_season$adjve_ill), lower = unique(temp_season$adjve_ill_lcl), upper = unique(temp_season$adjve_ill_ucl), 
                                  samples = nsim)
    
    # assign same VE to both illness and hospitalization
    temp_adj_sim_ve_ill <- c(temp_adj_sim_ve_ill, ve_illness_batch)
    temp_adj_sim_ve_hosp <- c(temp_adj_sim_ve_hosp, ve_illness_batch)
    
  } else {
    # Start sampling
    while (length(temp_adj_sim_ve_ill) < nsim) {
      
      # VE against illness estimates
      ve_illness_batch <- beta_pert(unique(temp_season$adjve_ill), lower = unique(temp_season$adjve_ill_lcl), upper = unique(temp_season$adjve_ill_ucl), 
                                    samples = nsim) 
      # VE against hospitalization estimates
      ve_hosp_batch <- beta_pert(unique(temp_season$adjve_hosp), lower = unique(temp_season$adjve_hosp_lcl), upper = unique(temp_season$adjve_hosp_ucl), 
                                 samples = nsim)
      
      df_batch <- tibble(ve_ill_sample = ve_illness_batch,
                         ve_hosp_sample = ve_hosp_batch)
      
      #discard if assumption not met
      df_acceptable <- df_batch %>% 
        filter(ve_hosp_sample > ve_ill_sample) %>% 
        mutate(ve_ill_sample = ifelse(ve_ill_sample < 0, 0, ve_ill_sample),
               ve_hosp_sample = ifelse(ve_hosp_sample < 0, 0, ve_hosp_sample))
      
      temp_adj_sim_ve_ill <- c(temp_adj_sim_ve_ill, df_acceptable$ve_ill_sample)
      temp_adj_sim_ve_hosp <- c(temp_adj_sim_ve_hosp, df_acceptable$ve_hosp_sample)
    }
  }
  
  # Ensure we have exactly nsim samples for this parameter set
  temp_adj_sim_ve_ill <- temp_adj_sim_ve_ill[1:nsim]
  temp_adj_sim_ve_hosp <- temp_adj_sim_ve_hosp[1:nsim]
  
  # Append the samples to the overall vectors
  adj_sim_ve_ill <- c(adj_sim_ve_ill, temp_adj_sim_ve_ill)
  adj_sim_ve_hosp <- c(adj_sim_ve_hosp, temp_adj_sim_ve_hosp)
  
}
ve_param <- bind_cols(
  adj_sim_ve_ill = adj_sim_ve_ill,
  adj_sim_ve_hosp = adj_sim_ve_hosp
  ) %>%
  mutate(sim_index = rep(1:nsim, max(temp_sim$year)+1),
         year = rep(unique(temp_sim$year), each = nsim))


# final sim data set
set.seed(123)
dat_sim <- temp_sim %>% 
  left_join(ratio_param %>% 
              select(sim_index, year, sim_hnhratio, sim_maratio, sim_dhratio),
            by = c("sim_index", "year")) %>% 
  left_join(vaccine_param %>% 
              select(sim_index, year, total_coverage, adj_sim_vc),
            by = c("sim_index", "year")) %>% 
  left_join(ve_param,
            by = c("sim_index", "year")) %>% 
  # simulated hospiatlization
  mutate(sim_hosp = rpois(n(), lambda = ceiling(mnth_hosp))) %>% 
  # monthly coverage
  mutate(
    propvax = mnth_coverage / total_coverage,
    propvax = ifelse(is.nan(propvax), 0, propvax),
    sim_mnthvc = adj_sim_vc * propvax
  ) %>% 
  group_by(sim_index, year) %>%
  mutate(starting_ve_ill = waning_ve(adj_sim_ve_ill[1])[[1]]) %>%
  ungroup() %>% 
  mutate(
    effvax_pop = target_pop_size * sim_mnthvc * starting_ve_ill/100,  
    noteff_pop = target_pop_size * sim_mnthvc * (1-starting_ve_ill/100)  
  ) 

##### averted illness #####

# current influenza vaccination
if (scenario == 1) { 
  
  avert <- dat_sim %>%
    group_by(sim_index) %>%
    group_modify(~ calculate_presence(.x)) %>%
    bind_rows() %>% 
    #hypothetical absence of vaccination program
    mutate(absence_hosp = hypo_hosp, 
           absence_nohosp = hypo_nohosp,
           absence_maill = hypo_maill,
           absence_death = hypo_death) %>% 
    #observed presence of vaccination program
    mutate(presence_hosp = unv_hosp + nev_hosp, 
           presence_nohosp = unv_nohosp + nev_nohosp,
           presence_maill = unv_maill + nev_maill,
           presence_death = unv_death + nev_death) %>% 
    #averted burden
    mutate(avert_hosp = absence_hosp - presence_hosp,
           avert_nohosp = absence_nohosp - presence_nohosp,
           avert_maill = absence_maill - presence_maill,
           avert_death = absence_death - presence_death) 
  
} else { #hypothetical influenza vaccination 
  
  avert <- dat_sim %>%
    group_by(sim_index) %>%
    group_modify(~ calculate_absence(.x)) %>%
    bind_rows() %>% 
    #observed absence of vaccination program
    mutate(absence_hosp = sim_hosp, 
           absence_nohosp = sim_nohosp,
           absence_maill = sim_maill,
           absence_death = sim_death) %>% 
    #hypothetical presence of vaccination program
    mutate(presence_hosp = hypo_unv_hosp + hypo_nev_hosp, 
           presence_nohosp = hypo_unv_nohosp + hypo_nev_nohosp,
           presence_maill = hypo_unv_maill + hypo_nev_maill,
           presence_death = hypo_unv_death + hypo_nev_death) %>% 
    #averted burden
    mutate(avert_hosp = absence_hosp - presence_hosp,
           avert_nohosp = absence_nohosp - presence_nohosp,
           avert_maill = absence_maill - presence_maill,
           avert_death = absence_death - presence_death) 
  
}

# program end 
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken


##### Inputs as it is #####
#get parameter data ready for calculation
run_param <- data_import %>% 
    mutate(sim_dhratio = ifelse(is.na(death_overall), mult_deathhosp, mult_deathhosp * 1/death_overall)) %>% 
    rename(sim_hosp = mnth_hosp,
           sim_mnthvc = mnth_coverage,
           sim_hnhratio = mult_nonhosp,
           sim_maratio = mult_ma,
           adj_sim_ve_ill = adjve_ill,
           adj_sim_ve_hosp = adjve_hosp)
  
# Apply the waning_ve function to dat_sim$adj_sim_ve_ill
run_param$starting_ve_ill <- sapply(run_param$adj_sim_ve_ill, function(report_ve) {
  starting_ve_ill <- waning_ve(report_ve)[[1]]  
  return(starting_ve_ill)
})
  
# effectively/not effectively population at each month
run_param <- run_param %>% 
  mutate(
    effvax_pop = target_pop_size * sim_mnthvc * starting_ve_ill/100,  
    noteff_pop = target_pop_size * sim_mnthvc * (1-starting_ve_ill/100)  
  ) 
  
# current influenza vaccination
if (scenario == 1) { 
    
  avert_param <- run_param %>%
      calculate_presence() %>%
      #hypothetical absence of vaccination program
      mutate(absence_hosp = hypo_hosp, 
             absence_nohosp = hypo_nohosp,
             absence_maill = hypo_maill,
             absence_death = hypo_death) %>% 
      #observed presence of vaccination program
      mutate(presence_hosp = unv_hosp + nev_hosp, 
             presence_nohosp = unv_nohosp + nev_nohosp,
             presence_maill = unv_maill + nev_maill,
             presence_death = unv_death + nev_death) %>% 
      #averted burden
      mutate(avert_hosp = absence_hosp - presence_hosp,
             avert_nohosp = absence_nohosp - presence_nohosp,
             avert_maill = absence_maill - presence_maill,
             avert_death = absence_death - presence_death) 
    
} else { #hypothetical influenza vaccination 
    
  avert_param <- run_param %>%
      calculate_absence() %>%
      #observed absence of vaccination program
      mutate(absence_hosp = sim_hosp, 
             absence_nohosp = sim_nohosp,
             absence_maill = sim_maill,
             absence_death = sim_death) %>% 
      #hypothetical presence of vaccination program
      mutate(presence_hosp = hypo_unv_hosp + hypo_nev_hosp, 
             presence_nohosp = hypo_unv_nohosp + hypo_nev_nohosp,
             presence_maill = hypo_unv_maill + hypo_nev_maill,
             presence_death = hypo_unv_death + hypo_nev_death) %>% 
      #averted burden
      mutate(avert_hosp = absence_hosp - presence_hosp,
             avert_nohosp = absence_nohosp - presence_nohosp,
             avert_maill = absence_maill - presence_maill,
             avert_death = absence_death - presence_death) 
    
}




##### output #####


# for each simulated dataset, calculate output
output <- avert %>% 
  group_by(sim_index, year) %>% 
  summarise(absence_hosp_sum = sum(absence_hosp, na.rm = TRUE),
            presence_hosp_sum = sum(presence_hosp, na.rm = TRUE),
            avert_hosp_sum = sum(avert_hosp, na.rm = TRUE),
            
            absence_nohosp_sum = sum(absence_nohosp, na.rm = TRUE),
            presence_nohosp_sum = sum(presence_nohosp, na.rm = TRUE),
            avert_nohosp_sum = sum(avert_nohosp, na.rm = TRUE),
            
            absence_maill_sum = sum(absence_maill, na.rm = TRUE),
            presence_maill_sum = sum(presence_maill, na.rm = TRUE),
            avert_maill_sum = sum(avert_maill, na.rm = TRUE),
            
            absence_death_sum = sum(absence_death, na.rm = TRUE),
            presence_death_sum = sum(presence_death, na.rm = TRUE),
            avert_death_sum = sum(avert_death, na.rm = TRUE),
            
            vax_sum = sum(effvax_pop, na.rm = TRUE) + sum(noteff_pop, na.rm = TRUE)) %>% 
  mutate(pf_hosp = avert_hosp_sum/absence_hosp_sum,
         pf_nohosp = avert_nohosp_sum/absence_nohosp_sum,
         pf_maill = avert_maill_sum/absence_maill_sum,
         pf_death = avert_death_sum/absence_death_sum,
         pf_overall = (avert_hosp_sum + avert_nohosp_sum)/(absence_hosp_sum + absence_nohosp_sum) ) %>% 
  mutate(nnv_hosp = vax_sum/avert_hosp_sum,
         nnv_nohosp = vax_sum/avert_nohosp_sum,
         nnv_maill = vax_sum/avert_maill_sum,
         nnv_death = vax_sum/avert_death_sum,
         nnv_overall = vax_sum/(avert_hosp_sum + avert_nohosp_sum) )

# point estimate with uncertainty 
result <- output %>% 
  group_by(year) %>% 
  summarise(
            #burden in presence of vaccine program
            absence_hosp_mean = mean(absence_hosp_sum, na.rm = TRUE),
            absence_hosp_lcl = quantile(absence_hosp_sum, 0.025, na.rm = TRUE),
            absence_hosp_ucl = quantile(absence_hosp_sum, 0.975, na.rm = TRUE),
            absence_nohosp_mean = mean(absence_nohosp_sum, na.rm = TRUE),
            absence_nohosp_lcl = quantile(absence_nohosp_sum, 0.025, na.rm = TRUE),
            absence_nohosp_ucl = quantile(absence_nohosp_sum, 0.975, na.rm = TRUE),
            absence_maill_mean = mean(absence_maill_sum, na.rm = TRUE),
            absence_maill_lcl = quantile(absence_maill_sum, 0.025, na.rm = TRUE),
            absence_maill_ucl = quantile(absence_maill_sum, 0.975, na.rm = TRUE),
            absence_death_mean = mean(absence_death_sum, na.rm = TRUE),
            absence_death_lcl = quantile(absence_death_sum, 0.025, na.rm = TRUE),
            absence_death_ucl = quantile(absence_death_sum, 0.975, na.rm = TRUE),
            #burden in presence of vaccine program
            presence_hosp_mean = mean(presence_hosp_sum, na.rm = TRUE),
            presence_hosp_lcl = quantile(presence_hosp_sum, 0.025, na.rm = TRUE),
            presence_hosp_ucl = quantile(presence_hosp_sum, 0.975, na.rm = TRUE),
            presence_nohosp_mean = mean(presence_nohosp_sum, na.rm = TRUE),
            presence_nohosp_lcl = quantile(presence_nohosp_sum, 0.025, na.rm = TRUE),
            presence_nohosp_ucl = quantile(presence_nohosp_sum, 0.975, na.rm = TRUE),
            presence_maill_mean = mean(presence_maill_sum, na.rm = TRUE),
            presence_maill_lcl = quantile(presence_maill_sum, 0.025, na.rm = TRUE),
            presence_maill_ucl = quantile(presence_maill_sum, 0.975, na.rm = TRUE),
            presence_death_mean = mean(presence_death_sum, na.rm = TRUE),
            presence_death_lcl = quantile(presence_death_sum, 0.025, na.rm = TRUE),
            presence_death_ucl = quantile(presence_death_sum, 0.975, na.rm = TRUE),
            #averted burden
            avert_hosp_mean = mean(avert_hosp_sum, na.rm = TRUE),
            avert_hosp_lcl = quantile(avert_hosp_sum, 0.025, na.rm = TRUE),
            avert_hosp_ucl = quantile(avert_hosp_sum, 0.975, na.rm = TRUE),
            avert_nohosp_mean = mean(avert_nohosp_sum, na.rm = TRUE),
            avert_nohosp_lcl = quantile(avert_nohosp_sum, 0.025, na.rm = TRUE),
            avert_nohosp_ucl = quantile(avert_nohosp_sum, 0.975, na.rm = TRUE),
            avert_maill_mean = mean(avert_maill_sum, na.rm = TRUE),
            avert_maill_lcl = quantile(avert_maill_sum, 0.025, na.rm = TRUE),
            avert_maill_ucl = quantile(avert_maill_sum, 0.975, na.rm = TRUE),
            avert_death_mean = mean(avert_death_sum, na.rm = TRUE),
            avert_death_lcl = quantile(avert_death_sum, 0.025, na.rm = TRUE),
            avert_death_ucl = quantile(avert_death_sum, 0.975, na.rm = TRUE),
            #prevented fraction (avert/burden in absence of vaccination)
            pf_hosp_mean = mean(pf_hosp, na.rm = TRUE),
            pf_hosp_lcl = quantile(pf_hosp, 0.025, na.rm = TRUE),
            pf_hosp_ucl = quantile(pf_hosp, 0.975, na.rm = TRUE),
            pf_nohosp_mean = mean(pf_nohosp, na.rm = TRUE),
            pf_nohosp_lcl = quantile(pf_nohosp, 0.025, na.rm = TRUE),
            pf_nohosp_ucl = quantile(pf_nohosp, 0.975, na.rm = TRUE),
            pf_maill_mean = mean(pf_maill, na.rm = TRUE),
            pf_maill_lcl = quantile(pf_maill, 0.025, na.rm = TRUE),
            pf_maill_ucl = quantile(pf_maill, 0.975, na.rm = TRUE),
            pf_death_mean = mean(pf_death, na.rm = TRUE),
            pf_death_lcl = quantile(pf_death, 0.025, na.rm = TRUE),
            pf_death_ucl = quantile(pf_death, 0.975, na.rm = TRUE),
            pf_overall_mean = mean(pf_overall, na.rm = TRUE),
            pf_overall_lcl = quantile(pf_overall, 0.025, na.rm = TRUE),
            pf_overall_ucl = quantile(pf_overall, 0.975, na.rm = TRUE),
            
            avert_overall_mean = mean(avert_hosp_sum + avert_nohosp_sum, na.rm = TRUE),
            avert_overall_lcl = quantile(avert_hosp_sum + avert_nohosp_sum, 0.025, na.rm = TRUE),
            avert_overall_ucl = quantile(avert_hosp_sum + avert_nohosp_sum, 0.975, na.rm = TRUE),
            
            #NNV (number vaccinated/number of cases averted)
            nnv_hosp_mean = mean(nnv_hosp, na.rm = TRUE),
            nnv_hosp_lcl = quantile(nnv_hosp, 0.025, na.rm = TRUE),
            nnv_hosp_ucl = quantile(nnv_hosp, 0.975, na.rm = TRUE),
            nnv_nohosp_mean = mean(nnv_nohosp, na.rm = TRUE),
            nnv_nohosp_lcl = quantile(nnv_nohosp, 0.025, na.rm = TRUE),
            nnv_nohosp_ucl = quantile(nnv_nohosp, 0.975, na.rm = TRUE),
            nnv_maill_mean = mean(nnv_maill, na.rm = TRUE),
            nnv_maill_lcl = quantile(nnv_maill, 0.025, na.rm = TRUE),
            nnv_maill_ucl = quantile(nnv_maill, 0.975, na.rm = TRUE),
            nnv_death_mean = mean(nnv_death, na.rm = TRUE),
            nnv_death_lcl = quantile(nnv_death, 0.025, na.rm = TRUE),
            nnv_death_ucl = quantile(nnv_death, 0.975, na.rm = TRUE),
            nnv_overall_mean = mean(nnv_overall, na.rm = TRUE),
            nnv_overall_lcl = quantile(nnv_overall, 0.025, na.rm = TRUE),
            nnv_overall_ucl = quantile(nnv_overall, 0.975, na.rm = TRUE)
            )

result_param <- avert_param %>% 
  mutate(year = case_when(month < 13 ~ 0,
                          month < 25 ~ 1,
                          month < 37 ~ 2,
                          month < 49 ~ 3, 
                          TRUE ~ 4)) %>%
  group_by(year) %>% 
  summarise(absence_hosp_param = sum(absence_hosp, na.rm = TRUE),
            presence_hosp_param = sum(presence_hosp, na.rm = TRUE),
            avert_hosp_param = sum(avert_hosp, na.rm = TRUE),
            
            absence_nohosp_param = sum(absence_nohosp, na.rm = TRUE),
            presence_nohosp_param = sum(presence_nohosp, na.rm = TRUE),
            avert_nohosp_param = sum(avert_nohosp, na.rm = TRUE),
            
            absence_maill_param = sum(absence_maill, na.rm = TRUE),
            presence_maill_param = sum(presence_maill, na.rm = TRUE),
            avert_maill_param = sum(avert_maill, na.rm = TRUE),
            
            absence_death_param = sum(absence_death, na.rm = TRUE),
            presence_death_param = sum(presence_death, na.rm = TRUE),
            avert_death_param = sum(avert_death, na.rm = TRUE),
            
            vax_param = sum(effvax_pop, na.rm = TRUE) + sum(noteff_pop, na.rm = TRUE)) %>%  
  mutate(pf_hosp_param = avert_hosp_param/absence_hosp_param,
         pf_nohosp_param = avert_nohosp_param/absence_nohosp_param,
         pf_maill_param = avert_maill_param/absence_maill_param,
         pf_death_param = avert_death_param/absence_death_param,
         pf_overall_param = (avert_hosp_param + avert_nohosp_param) / (absence_hosp_param + absence_nohosp_param) ) %>%  
  mutate(nnv_hosp_param = vax_param/avert_hosp_param,
         nnv_nohosp_param = vax_param/avert_nohosp_param,
         nnv_maill_param = vax_param/avert_maill_param,
         nnv_death_param = vax_param/avert_death_param,
         nnv_overall_param = vax_param / (avert_hosp_param + avert_nohosp_param) ) %>% 
  select(-vax_param)

result_final <- result %>% 
  left_join(result_param, 
            by = "year") %>% 
  filter(year != 0)

# export results and scenarios to input tool
wb <- loadWorkbook("Phase 2 data input tool_final.xlsm")
if (!"R Output" %in% names(wb)) {
  addWorksheet(wb, sheetName = "R Output")
}
writeData(wb, sheet = "R Output", x = result_final, withFilter = TRUE)
setRowHeights(wb, sheet = "Inputs", rows = 20:33, heights = 0)
saveWorkbook(wb, file = "Phase 2 data input tool_final.xlsm", overwrite = TRUE)

