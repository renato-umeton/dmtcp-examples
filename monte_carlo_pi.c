#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main() {
    long long num_points = 500000000;
    long long inside_circle = 0;
    double x, y;
    
    printf("Starting Monte Carlo Pi estimation with %lld points...\n", num_points);
    
    // Seed the random number generator
    srand(time(NULL));
    
    // Monte Carlo simulation
    for (long long i = 0; i < num_points; i++) {
        // Generate random point in [0,1] x [0,1]
        x = (double)rand() / RAND_MAX;
        y = (double)rand() / RAND_MAX;
        
        // Check if point is inside quarter circle
        if (x*x + y*y <= 1.0) {
            inside_circle++;
        }
        
        // Progress updates every 50 million points
        if (i % 50000000 == 0 && i > 0) {
            printf("Progress: %lld points processed (%.1f%%)\n", 
                   i, (double)i / num_points * 100);
        }
    }
    
    // Calculate Pi estimate
    double pi_estimate = 4.0 * inside_circle / num_points;
    
    printf("\n=== Results ===\n");
    printf("Estimated Pi: %.8f\n", pi_estimate);
    printf("Actual Pi:    %.8f\n", M_PI);
    printf("Error:        %.8f\n", fabs(pi_estimate - M_PI));
    printf("Used %lld points\n", num_points);
    
    return 0;
}
