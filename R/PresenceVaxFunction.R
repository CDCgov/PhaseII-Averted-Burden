
# function to calculate hypothetical hospitalizations
calculate_presence <- function(dat) {
  
  # list of matrix to fill in: observed scenario (unvaccinated)
  remaining_unvS <- numeric(time_period)
  unv_hosp <- numeric(time_period)
  unv_nohosp <- numeric(time_period)
  unv_maill <- numeric(time_period)
  unv_death <- numeric(time_period)
  
  # list of matrix to fill in: observed scenario (not effectively vaccinated)
  remaining_nevS <- matrix(NA, nrow = time_period, ncol = time_period)
  nev_hosp_matrix <- matrix(NA, nrow = time_period, ncol = time_period)
  nev_nohosp_matrix <- matrix(NA, nrow = time_period, ncol = time_period)
  nev_maill_matrix <- matrix(NA, nrow = time_period, ncol = time_period)
  nev_death_matrix <- matrix(NA, nrow = time_period, ncol = time_period)
  
  # list of matrix to fill in: observed scenario (hospitalization rate)
  hosp_rate <- numeric(time_period)
  
  # list of matrix to fill in: hypothetical scenario
  hypo_remaining_S <- numeric(time_period)
  hypo_hosp <- numeric(time_period)
  hypo_nohosp <- numeric(time_period)
  hypo_maill <- numeric(time_period)
  hypo_death <- numeric(time_period)
  
  # Set the initial value for observed scenario
  remaining_unvS[1] <- dat$target_pop_size[1] - dat$sim_hosp[1] - (dat$sim_hosp[1] * dat$sim_hnhratio[1]) - dat$effvax_pop[1] - dat$noteff_pop[1]
  unv_hosp[1] <- dat$sim_hosp[1]
  unv_nohosp[1] <- dat$sim_hosp[1] * dat$sim_hnhratio[1]
  unv_maill[1] <- (unv_hosp[1] + unv_nohosp[1]) * dat$sim_maratio[1]
  unv_death[1] <- unv_hosp[1] * dat$sim_dhratio[1]
  diag(remaining_nevS) <- dat$noteff_pop
  diag(nev_hosp_matrix) <- 0
  diag(nev_nohosp_matrix) <- 0
  diag(nev_maill_matrix) <- 0
  diag(nev_death_matrix) <- 0
  hosp_rate[1] <- dat$sim_hosp[1] / dat$target_pop_size[1]
  
  # Set the initial value for hypothetical scenario
  hypo_remaining_S[1] <- dat$target_pop_size[1] - dat$sim_hosp[1] - (dat$sim_hosp[1] * dat$sim_hnhratio[1])
  hypo_hosp[1] <- dat$sim_hosp[1]
  hypo_nohosp[1] <- dat$sim_hosp[1] * dat$sim_hnhratio[1]
  hypo_maill[1] <- (hypo_hosp[1] + hypo_nohosp[1]) * dat$sim_maratio[1]
  hypo_death[1] <- hypo_hosp[1] * dat$sim_dhratio[1]
  
  # effectively vaccinated 
  ev_return <- matrix(NA, nrow = time_period, ncol = time_period)
  ve_value <- sapply(1:time_period, function(d) waning_ve(dat$adj_sim_ve_ill[d])[[3]] / 100)
  for (m in 1:time_period) {
    
    for (i in 1:time_period){
      if (i >= m) {
        ve_index <- i - m + 1
        ev_return[i, m] <- dat$target_pop_size[i] * dat$sim_mnthvc[m] * ve_value[ve_index, m] 
      }
    }
    
  }
  
  # Fill the rest
  waning_ve_ill <- sapply(1:time_period, function(d) waning_ve(dat$adj_sim_ve_ill[d])[[2]] / 100)
  waning_ve_hosp <- sapply(1:time_period, function(d) waning_ve(dat$adj_sim_ve_hosp[d])[[2]] / 100)
  residual_protection <- (1 - waning_ve_hosp) / (1 - waning_ve_ill)
  for (m in 2:time_period) {
    
    # derive hospitalization rate lambda
    denom <- remaining_unvS[m-1]
    for (i in 1:(m-1)) {  
      # denominator
      denom <- denom + (remaining_nevS[m-1, i] * residual_protection[m-i, i])
    }
    hosp_rate[m] <- dat$sim_hosp[m] / denom
    
    # Not effectively vaccinated                          
    for (i in 1:(m-1)) {
      nev_hosp_matrix[m, i] <- remaining_nevS[m-1, i] * hosp_rate[m] * residual_protection[m-i, i] 
      nev_nohosp_matrix[m, i] <- remaining_nevS[m-1, i] * (hosp_rate[m] * dat$sim_hnhratio[m])
      nev_maill_matrix[m, i] <- (nev_hosp_matrix[m, i] + nev_nohosp_matrix[m, i]) * dat$sim_maratio[m]
      nev_death_matrix[m, i] <- nev_hosp_matrix[m, i] * dat$sim_dhratio[m]
      remaining_nevS[m, i] <- remaining_nevS[m-1, i] - nev_hosp_matrix[m, i] - nev_nohosp_matrix[m, i]
    }
    
    # Unvaccinated
    unv_hosp[m] <- remaining_unvS[m-1] * hosp_rate[m]  
    unv_nohosp[m] <- remaining_unvS[m-1] * (hosp_rate[m] * dat$sim_hnhratio[m])
    unv_maill[m] <- (unv_hosp[m] + unv_nohosp[m]) * dat$sim_maratio[m]
    unv_death[m] <- unv_hosp[m] * dat$sim_dhratio[m]
    remaining_unvS[m] <- remaining_unvS[m-1] - 
      # Cases at month m
      unv_hosp[m] - unv_nohosp[m] -
      # vaccinated people at month m
      dat$effvax_pop[m] - dat$noteff_pop[m] +
      # return to Unvaccinated susceptible due to waning vaccine-induced immunity
      sum(ev_return[m, ], na.rm = TRUE) +
      # return to Unvaccinated susceptible as their natural immunity expires after 8 months
      (ifelse(m > 8, unv_hosp[m-8], 0) + ifelse(m > 8, unv_nohosp[m-8], 0)) + 
      # (ifelse(m > 8, unv_hosp[m-8], 0) + ifelse(m > 8, unv_nohosp[m-8], 0)- ifelse(m > 8, unv_death[m-8], 0)) + 
      (ifelse(m > 8, sum(nev_hosp_matrix[m-8, ], na.rm = TRUE), 0) + ifelse(m > 8, sum(nev_nohosp_matrix[m-8, ], na.rm = TRUE), 0))
      # (ifelse(m > 8, sum(nev_hosp_matrix[m-8, ], na.rm = TRUE), 0) + ifelse(m > 8, sum(nev_nohosp_matrix[m-8, ], na.rm = TRUE), 0) - ifelse(m > 8, sum(nev_death_matrix[m-8, ], na.rm = TRUE), 0))
      
                                  
    # Hypothetical scenario: No vaccine
    hypo_hosp[m] <- hosp_rate[m] * hypo_remaining_S[m-1]
    hypo_nohosp[m] <- (hosp_rate[m] * dat$sim_hnhratio[m]) * hypo_remaining_S[m-1]
    hypo_maill[m] <- (hypo_hosp[m] + hypo_nohosp[m]) * dat$sim_maratio[m]
    hypo_death[m] <- hypo_hosp[m] * dat$sim_dhratio[m]
    hypo_remaining_S[m] <- hypo_remaining_S[m-1] - hypo_hosp[m] - hypo_nohosp[m] +
      (ifelse(m > 8, hypo_hosp[m-8], 0) + ifelse(m > 8, hypo_nohosp[m-8], 0))  # 8-month protection (natural immunity)
      # (ifelse(m > 8, hypo_hosp[m-8], 0) + ifelse(m > 8, hypo_nohosp[m-8], 0) - ifelse(m > 8, hypo_death[m-8], 0))  # 8-month protection (natural immunity)
      
  }
  
    
  # Convert the matrix to a data frame and assign column names
  unv_hosp_df <- as.data.frame(unv_hosp)
  colnames(unv_hosp_df) <- "unv_hosp"
  unv_nohosp_df <- as.data.frame(unv_nohosp)
  colnames(unv_nohosp_df) <- "unv_nohosp"
  unv_maill_df <- as.data.frame(unv_maill)
  colnames(unv_maill_df) <- "unv_maill"
  unv_death_df <- as.data.frame(unv_death)
  colnames(unv_death_df) <- "unv_death"
  
  nev_hosp_df <- as.data.frame(nev_hosp_matrix)
  colnames(nev_hosp_df) <- paste0("nev_hosp_m", 1:time_period)
  nev_nohosp_df <- as.data.frame(nev_nohosp_matrix)
  colnames(nev_nohosp_df) <- paste0("nev_nohosp_m", 1:time_period)
  nev_maill_df <- as.data.frame(nev_maill_matrix)
  colnames(nev_maill_df) <- paste0("nev_maill_m", 1:time_period)
  nev_death_df <- as.data.frame(nev_death_matrix)
  colnames(nev_death_df) <- paste0("nev_death_m", 1:time_period)
  
  hypo_hosp_df <- as.data.frame(hypo_hosp)
  colnames(hypo_hosp_df) <- "hypo_hosp"
  hypo_nohosp_df <- as.data.frame(hypo_nohosp)
  colnames(hypo_nohosp_df) <- "hypo_nohosp"
  hypo_maill_df <- as.data.frame(hypo_maill)
  colnames(hypo_maill_df) <- "hypo_maill"
  hypo_death_df <- as.data.frame(hypo_death)
  colnames(hypo_death_df) <- "hypo_death"
  
  # Bind columns to the original data
  result <- bind_cols(dat, 
                      unv_hosp_df, unv_nohosp_df, unv_maill_df, unv_death_df,
                      nev_hosp_df, nev_nohosp_df, nev_maill_df, nev_death_df,
                      hypo_hosp_df, hypo_nohosp_df, hypo_maill_df, hypo_death_df) %>% 
    mutate(nev_hosp = rowSums(across(starts_with("nev_hosp_m")), na.rm = TRUE),
           nev_nohosp = rowSums(across(starts_with("nev_nohosp_m")), na.rm = TRUE),
           nev_maill = rowSums(across(starts_with("nev_maill_m")), na.rm = TRUE),
           nev_death = rowSums(across(starts_with("nev_death_m")), na.rm = TRUE)) 
  
  return(result)
}

  
  