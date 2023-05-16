setwd('path/to/data')

data <- read.csv('flower_tampon_calibrated_initial_slopes.csv')

# calculate estimated concentration from initial slopes based on formula from
# 'flower_cystatin_a_calibration_curve.xlsx'
calculate_estimated_concentration <- function(initial_slope) {
  conc <- (0.1619 * initial_slope) + 0.7716
  return(10**conc)
}

estimated_concentration <- c()
estimated_concentration_lower_ci <- c()
estimated_concentration_upper_ci <- c()
for (i in 1:nrow(data)) {
  estimated_concentration <- c(estimated_concentration, calculate_estimated_concentration(data[i,]$initial_slopes))
  estimated_concentration_lower_ci <- c(estimated_concentration_lower_ci, calculate_estimated_concentration(data[i,]$lower_ci))
  estimated_concentration_upper_ci <- c(estimated_concentration_upper_ci, calculate_estimated_concentration(data[i,]$upper_ci))
}
data$estimated_concentration <- estimated_concentration
data$estimated_concentration_lower_ci <- estimated_concentration_lower_ci
data$estimated_concentration_upper_ci <- estimated_concentration_upper_ci

write.csv(data, file='flower_tampon_calibrated_initial_slopes_estimated_concentrations.csv', row.names=FALSE)
