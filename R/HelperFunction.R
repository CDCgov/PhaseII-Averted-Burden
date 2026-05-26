# Helper function for waning VE
waning_ve <- function(report_ve, num_months = time_period) {
  
  # convert to % if not already
  report_ve <- ifelse(report_ve < 1, report_ve * 100, report_ve)
  
  # Define the cubic decline VE function
  cubic_ve <- function(starting_ve, months) {
    ve <- starting_ve + (-1.37 * 2 * months) + (0.18 * 2^2 * months^2) - (0.03 * 2^3 * months^3)
    pmax(ve, 0)  # Set negative VE values to 0
  }
  
  # Calculate the average VE for a given starting VE
  calculate_average_ve <- function(starting_ve) {
    months <- 0:7  # 8 months (0 to 7)
    ve_values <- cubic_ve(starting_ve, months)
    mean(ve_values)  # Average VE over 8 months
  }
  
  # Optimize starting VE to minimize the difference between the average VE and target VE
  optimal_starting_ve <- optimize(
    f = function(starting_ve) abs(calculate_average_ve(starting_ve) - report_ve),
    interval = c(0, 100)  # Search range for starting VE
  )
  
  # VE for each month 
  months <- 0:(num_months-1)   # set to the number of months user specified
  ve_months <- cubic_ve(optimal_starting_ve$minimum, months)
  
  # VE difference between current month and previous month (how much did VE wane)
  ve_diff <- c(0, abs(diff(ve_months)))
  
  # Return the adjusted starting VE
  return(list(optimal_starting_ve$minimum, ve_months, ve_diff))
  
}

# Helper function to draw samples from beta PERT distribution
beta_pert <- function(most_likely, lower, upper, samples) {
  # Calculate mean 
  mean_pert <- (lower + 4 * most_likely + upper) / 6
  
  # alpha and beta 
  alpha <- 1 + 4 * (mean_pert - lower) / (upper - lower)
  beta <- 1 + 4 * (upper - mean_pert) / (upper - lower)
  
  # Draw sample from the beta distribution
  sample <- rbeta(samples, alpha, beta)

  # Rescale the beta distribution to [lower, upper]
  pert_sample <- lower + sample * (upper - lower)
  
  return(pert_sample)
}
