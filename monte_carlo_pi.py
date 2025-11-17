import random
import time
import math

def estimate_pi(num_points):
    """
    Estimate Pi using Monte Carlo simulation.
    """
    print(f"Starting Monte Carlo Pi estimation with {num_points:,} points...")
    
    inside_circle = 0
    
    for i in range(num_points):
        # Generate random point in [0,1] x [0,1]
        x = random.random()
        y = random.random()
        
        # Check if point is inside quarter circle
        if x**2 + y**2 <= 1.0:
            inside_circle += 1
        
        # Progress updates every 50 million points
        if i % 50000000 == 0 and i > 0:
            progress = (i / num_points) * 100
            print(f"Progress: {i:,} points processed ({progress:.1f}%)")
    
    # Calculate Pi estimate
    pi_estimate = 4.0 * inside_circle / num_points
    
    return pi_estimate, inside_circle

if __name__ == "__main__":
    # Set random seed for reproducibility
    random.seed(42)
    
    # Number of random points
    num_points = 500000000
    
    # Run simulation
    start_time = time.time()
    pi_estimate, inside_circle = estimate_pi(num_points)
    elapsed_time = time.time() - start_time
    
    # Print results
    print("\n=== Results ===")
    print(f"Estimated Pi: {pi_estimate:.8f}")
    print(f"Actual Pi:    {math.pi:.8f}")
    print(f"Error:        {abs(pi_estimate - math.pi):.8f}")
    print(f"Used {num_points:,} points")
    print(f"Points inside circle: {inside_circle:,}")
    print(f"Elapsed time: {elapsed_time:.2f} seconds")
