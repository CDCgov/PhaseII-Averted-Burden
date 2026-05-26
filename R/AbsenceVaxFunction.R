
# Function to calculate averted burden in the absence of a vaccine program 
calculate_absence <- function(dat) {
  
  ### Observed
  # list of matrix to fill in: (unvaccinated)
  remaining_S <- numeric(time_period)
  sim_nohosp <- numeric(time_period)
  sim_maill <- numeric(time_period)
  sim_death <- numeric(time_period)
  
  # list of matrix to fill in: observed scenario (hospitalization rate)
  hosp_rate <- numeric(time_period)
  
  # list of matrix to fill in: hypothetical scenario
  hypo_remaining_unvS <- numeric(time_period)
  hypo_unv_hosp <- numeric(time_period)
  hypo_unv_nohosp <- numeric(time_period)
  hypo_unv_maill <- numeric(time_period)
  hypo_unv_death <- numeric(time_period)
  hypo_remaining_nevS <- matrix(NA, nrow = time_period, ncol = time_period)
  hypo_nev_hosp_matrix <- matrix(NA, nrow = time_period, ncol = time_period)
  hypo_nev_nohosp_matrix <- matrix(NA, nrow = time_period, ncol = time_period)
  hypo_nev_maill_matrix <- matrix(NA, nrow = time_period, ncol = time_period)
  hypo_nev_death_matrix <- matrix(NA, nrow = time_period, ncol = time_period)
  
  # set the initial value for observed
  remaining_S[1] <- dat$target_pop_size[1] - dat$sim_hosp[1] - (dat$sim_hosp[1] * dat$sim_hnhratio[1])
  sim_nohosp[1] <- dat$sim_hosp[1] * dat$sim_hnhratio[1]
  sim_maill[1] <- (dat$sim_hosp[1] + sim_nohosp[1]) * dat$sim_maratio[1]
  sim_death[1] <- dat$sim_hosp[1] * dat$sim_dhratio[1]
  hosp_rate[1] <- dat$sim_hosp[1] / dat$target_pop_size[1]
  
  # set the initial value for hypothetical scenario
  hypo_remaining_unvS[1] <- dat$target_pop_size[1] - dat$sim_hosp[1] - sim_nohosp[1] - dat$effvax_pop[1] - dat$noteff_pop[1]
  hypo_unv_hosp[1] <- dat$sim_hosp[1]
  hypo_unv_nohosp[1] <- sim_nohosp[1]
  hypo_unv_maill[1] <- (hypo_unv_hosp[1] + hypo_unv_nohosp[1]) * dat$sim_maratio[1]
  hypo_unv_death[1] <- dat$sim_hosp[1] * dat$sim_dhratio[1]
  diag(hypo_remaining_nevS) <- dat$noteff_pop
  diag(hypo_nev_hosp_matrix) <- 0
  diag(hypo_nev_nohosp_matrix) <- 0
  diag(hypo_nev_maill_matrix) <- 0
  diag(hypo_nev_death_matrix) <- 0
  
  # effectively vaccinated 
  hypo_ev_return <- matrix(NA, nrow = time_period, ncol = time_period)
  ve_value <- sapply(1:time_period, function(d) waning_ve(dat$adj_sim_ve_ill[d])[[3]] / 100)
  for (m in 1:time_period) {
    
    for (i in 1:time_period){
      if (i >= m) {
        ve_index <- i - m + 1
        hypo_ev_return[i, m] <- dat$target_pop_size[i] * dat$sim_mnthvc[m] * ve_value[ve_index, m] 
      }
    }
  }
  
  # Calculate initial infections and track the remaining population for each month
  waning_ve_ill <- sapply(1:time_period, function(d) waning_ve(dat$adj_sim_ve_ill[d])[[2]] / 100)
  waning_ve_hosp <- sapply(1:time_period, function(d) waning_ve(dat$adj_sim_ve_hosp[d])[[2]] / 100)
  residual_protection <- (1 - waning_ve_hosp) / (1 - waning_ve_ill)
  for (m in 2:time_period) {
    
    ### Observed
    # hospitalization (case divided by non-vaccinated Susceptible)
    hosp_rate[m] <- dat$sim_hosp[m] / remaining_S[m-1]
    sim_nohosp[m] <- (hosp_rate[m] * dat$sim_hnhratio[m]) * remaining_S[m-1]
    sim_maill[m] <- (dat$sim_hosp[m] + sim_nohosp[m]) * dat$sim_maratio[m]
    sim_death[m] <- (dat$sim_hosp[m]) * dat$sim_dhratio[m]
    remaining_S[m] <- remaining_S[m-1] - dat$sim_hosp[m] - sim_nohosp[m] +
      (ifelse(m > 8, dat$sim_hosp[m-8], 0) + ifelse(m > 8, sim_nohosp[m-8], 0))   # 8-months protection (natural immunity)
      # (ifelse(m > 8, dat$sim_hosp[m-8], 0) + ifelse(m > 8, sim_nohosp[m-8], 0) - ifelse(m > 8, sim_death[m-8], 0))   # 8-months protection (natural immunity)
      
    
    ### Hypothetical scenario
    # not effectively vaccinated
    for (i in 1:(m-1)) {
      
      hypo_nev_hosp_matrix[m, i] <- hypo_remaining_nevS[m-1, i] * hosp_rate[m] * residual_protection[m-i, i]  
      hypo_nev_nohosp_matrix[m, i] <- hypo_remaining_nevS[m-1, i] * (hosp_rate[m] * dat$sim_hnhratio[m])
      hypo_nev_maill_matrix[m, i] <- (hypo_nev_hosp_matrix[m, i] + hypo_nev_nohosp_matrix[m, i]) * dat$sim_maratio[m]
      hypo_nev_death_matrix[m, i] <- hypo_nev_hosp_matrix[m, i] * dat$sim_dhratio[m]
      hypo_remaining_nevS[m, i] <- hypo_remaining_nevS[m-1, i] - hypo_nev_hosp_matrix[m, i] - hypo_nev_nohosp_matrix[m, i] 
    }
    
    # Unvaccinated
    hypo_unv_hosp[m] <- hypo_remaining_unvS[m-1] * hosp_rate[m]  
    hypo_unv_nohosp[m] <- hypo_remaining_unvS[m-1] * (hosp_rate[m] * dat$sim_hnhratio[m])
    hypo_unv_maill[m] <- (hypo_unv_hosp[m] + hypo_unv_nohosp[m]) * dat$sim_maratio[m]
    hypo_unv_death[m] <- hypo_unv_hosp[m] * dat$sim_dhratio[m]
    hypo_remaining_unvS[m] <- hypo_remaining_unvS[m-1] - 
      # Cases at month m
      hypo_unv_hosp[m] - hypo_unv_nohosp[m] -
      # vaccinated people at month m
      dat$effvax_pop[m] - dat$noteff_pop[m] +
      # return to Unvaccinated susceptible due to waning vaccine-induced immunity
      sum(hypo_ev_return[m, ], na.rm = TRUE) + 
      # return to Unvaccinated susceptible as their natural immunity expires after 8 months
      (ifelse(m > 8, hypo_unv_hosp[m-8], 0) + ifelse(m > 8, hypo_unv_nohosp[m-8], 0) ) +
      # (ifelse(m > 8, hypo_unv_hosp[m-8], 0) + ifelse(m > 8, hypo_unv_nohosp[m-8], 0) - ifelse(m > 8, hypo_unv_death[m-8], 0)) +
      (ifelse(m > 8, sum(hypo_nev_hosp_matrix[m-8, ], na.rm = TRUE), 0) + ifelse(m > 8, sum(hypo_nev_nohosp_matrix[m-8, ], na.rm = TRUE), 0) )
      # (ifelse(m > 8, sum(hypo_nev_hosp_matrix[m-8, ], na.rm = TRUE), 0) + ifelse(m > 8, sum(hypo_nev_nohosp_matrix[m-8, ], na.rm = TRUE), 0) - ifelse(m > 8, sum(hypo_nev_death_matrix[m-8, ], na.rm = TRUE), 0))
      
    
  }
  
  # Convert the matrix to a data frame and assign column names (observed)
  sim_nohosp_df <- as.data.frame(sim_nohosp)
  colnames(sim_nohosp_df) <- "sim_nohosp"
  sim_maill_df <- as.data.frame(sim_maill)
  colnames(sim_maill_df) <- "sim_maill"
  sim_death_df <- as.data.frame(sim_death)
  colnames(sim_death_df) <- "sim_death"
  
  # Convert the matrix to a data frame and assign column names (hypothetical nonvaccinated)
  hypo_unv_hosp_df <- as.data.frame(hypo_unv_hosp)
  colnames(hypo_unv_hosp_df) <- "hypo_unv_hosp"
  hypo_unv_nohosp_df <- as.data.frame(hypo_unv_nohosp)
  colnames(hypo_unv_nohosp_df) <- "hypo_unv_nohosp"
  hypo_unv_maill_df <- as.data.frame(hypo_unv_maill)
  colnames(hypo_unv_maill_df) <- "hypo_unv_maill"
  hypo_unv_death_df <- as.data.frame(hypo_unv_death)
  colnames(hypo_unv_death_df) <- "hypo_unv_death"
  
  # Convert the matrix to a data frame and assign column names (hypothetical not effectivley vaccinated)
  hypo_nev_hosp_matrix_df <- as.data.frame(hypo_nev_hosp_matrix)
  colnames(hypo_nev_hosp_matrix_df) <- paste0("hypo_nev_hosp_m", 1:time_period)
  hypo_nev_nohosp_matrix_df <- as.data.frame(hypo_nev_nohosp_matrix)
  colnames(hypo_nev_nohosp_matrix_df) <- paste0("hypo_nev_nohosp_m", 1:time_period)
  hypo_nev_maill_matrix_df <- as.data.frame(hypo_nev_maill_matrix)
  colnames(hypo_nev_maill_matrix_df) <- paste0("hypo_nev_maill_m", 1:time_period)
  hypo_nev_death_matrix_df <- as.data.frame(hypo_nev_death_matrix)
  colnames(hypo_nev_death_matrix_df) <- paste0("hypo_nev_death_m", 1:time_period)
  
  # Bind columns to the original data
  result <- bind_cols(dat, sim_nohosp_df, sim_maill_df, sim_death_df, 
                      hypo_unv_hosp_df, hypo_unv_nohosp_df, hypo_unv_maill_df, hypo_unv_death_df, 
                      hypo_nev_hosp_matrix_df, hypo_nev_nohosp_matrix_df, hypo_nev_maill_matrix_df, hypo_nev_death_matrix_df) %>% 
    mutate(hypo_nev_hosp = rowSums(across(starts_with("hypo_nev_hosp_m")), na.rm = TRUE),
           hypo_nev_nohosp = rowSums(across(starts_with("hypo_nev_nohosp_m")), na.rm = TRUE),
           hypo_nev_maill = rowSums(across(starts_with("hypo_nev_maill_m")), na.rm = TRUE),
           hypo_nev_death = rowSums(across(starts_with("hypo_nev_death_m")), na.rm = TRUE)) 
  
  return(result)
  
}