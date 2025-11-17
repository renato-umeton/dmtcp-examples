# Monte Carlo Pi Estimation in R
# This script estimates Pi using random sampling

estimate_pi <- function(num_points) {
  cat("Starting Monte Carlo Pi estimation with", format(num_points, big.mark=","), "points...\n")
  
  inside_circle <- 0
  
  for (i in 1:num_points) {
    # Generate random point in [0,1] x [0,1]
    x <- runif(1, 0, 1)
    y <- runif(1, 0, 1)
    
    # Check if point is inside quarter circle
    if (x^2 + y^2 <= 1) {
      inside_circle <- inside_circle + 1
    }
    
    # Progress updates every 50 million points
    if (i %% 50000000 == 0 && i > 0) {
      progress <- (i / num_points) * 100
      cat(sprintf("Progress: %s points processed (%.1f%%)\n", 
                  format(i, big.mark=","), progress))
    }
  }
  
  # Calculate Pi estimate
  pi_estimate <- 4 * inside_circle / num_points
  
  return(list(
    pi_estimate = pi_estimate,
    inside_circle = inside_circle,
    num_points = num_points
  ))
}

# Main execution
set.seed(42)  # For reproducibility
num_points <- 500000000

# Record start time
start_time <- Sys.time()

# Run simulation
result <- estimate_pi(num_points)

# Record end time
end_time <- Sys.time()
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Print results
cat("\n=== Results ===\n")
cat(sprintf("Estimated Pi: %.8f\n", result$pi_estimate))
cat(sprintf("Actual Pi:    %.8f\n", pi))
cat(sprintf("Error:        %.8f\n", abs(result$pi_estimate - pi)))
cat(sprintf("Used %s points\n", format(result$num_points, big.mark=",")))
cat(sprintf("Points inside circle: %s\n", format(result$inside_circle, big.mark=",")))
cat(sprintf("Elapsed time: %.2f seconds\n", elapsed_time))
