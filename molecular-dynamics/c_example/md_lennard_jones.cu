/*
 * Molecular Dynamics Simulation: Lennard-Jones Fluid with DMTCP Checkpointing
 * Uses CUDA for GPU acceleration
 * Compile: nvcc -o md_lj md_lennard_jones.cu -lm
 * Run with DMTCP: dmtcp_launch ./md_lj
 * Checkpoint: dmtcp_command --checkpoint
 * Restart: dmtcp_restart ckpt_*.dmtcp
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <time.h>

#define N_PARTICLES 8192        // Number of particles
#define BOX_SIZE 50.0          // Simulation box size (Angstroms)
#define DT 0.001               // Time step (ps)
#define N_STEPS 100000         // Total simulation steps
#define CHECKPOINT_FREQ 1000   // Steps between checkpoints
#define CUTOFF 2.5             // LJ cutoff distance
#define EPSILON 0.238          // LJ epsilon (kcal/mol)
#define SIGMA 3.4              // LJ sigma (Angstroms)

typedef struct {
    float x, y, z;
} Vec3;

// CUDA kernel for force calculation using Lennard-Jones potential
__global__ void calculate_forces(Vec3* positions, Vec3* forces, int n, float box_size) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n) return;
    
    float fx = 0.0f, fy = 0.0f, fz = 0.0f;
    Vec3 pi = positions[i];
    
    for (int j = 0; j < n; j++) {
        if (i == j) continue;
        
        Vec3 pj = positions[j];
        
        // Minimum image convention
        float dx = pi.x - pj.x;
        float dy = pi.y - pj.y;
        float dz = pi.z - pj.z;
        
        dx -= box_size * roundf(dx / box_size);
        dy -= box_size * roundf(dy / box_size);
        dz -= box_size * roundf(dz / box_size);
        
        float r2 = dx*dx + dy*dy + dz*dz;
        
        if (r2 < CUTOFF * CUTOFF && r2 > 0.0f) {
            float r2inv = 1.0f / r2;
            float r6inv = r2inv * r2inv * r2inv;
            float sigma_r6 = SIGMA * SIGMA * SIGMA * SIGMA * SIGMA * SIGMA * r6inv;
            float sigma_r12 = sigma_r6 * sigma_r6;
            
            // Force magnitude: 24*epsilon * (2*sigma^12/r^14 - sigma^6/r^8)
            float force_mag = 24.0f * EPSILON * r2inv * (2.0f * sigma_r12 - sigma_r6);
            
            fx += force_mag * dx;
            fy += force_mag * dy;
            fz += force_mag * dz;
        }
    }
    
    forces[i].x = fx;
    forces[i].y = fy;
    forces[i].z = fz;
}

// CUDA kernel for velocity Verlet integration
__global__ void integrate_positions(Vec3* positions, Vec3* velocities, Vec3* forces, 
                                   int n, float dt, float box_size, float mass) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n) return;
    
    // Update positions: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
    positions[i].x += velocities[i].x * dt + 0.5f * forces[i].x / mass * dt * dt;
    positions[i].y += velocities[i].y * dt + 0.5f * forces[i].y / mass * dt * dt;
    positions[i].z += velocities[i].z * dt + 0.5f * forces[i].z / mass * dt * dt;
    
    // Apply periodic boundary conditions
    positions[i].x -= box_size * floorf(positions[i].x / box_size);
    positions[i].y -= box_size * floorf(positions[i].y / box_size);
    positions[i].z -= box_size * floorf(positions[i].z / box_size);
    
    // Half-step velocity update: v(t+dt/2) = v(t) + 0.5*a(t)*dt
    velocities[i].x += 0.5f * forces[i].x / mass * dt;
    velocities[i].y += 0.5f * forces[i].y / mass * dt;
    velocities[i].z += 0.5f * forces[i].z / mass * dt;
}

__global__ void update_velocities(Vec3* velocities, Vec3* forces, 
                                 int n, float dt, float mass) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n) return;
    
    // Complete velocity update: v(t+dt) = v(t+dt/2) + 0.5*a(t+dt)*dt
    velocities[i].x += 0.5f * forces[i].x / mass * dt;
    velocities[i].y += 0.5f * forces[i].y / mass * dt;
    velocities[i].z += 0.5f * forces[i].z / mass * dt;
}

void initialize_system(Vec3* h_positions, Vec3* h_velocities) {
    srand(42);
    
    // Place particles on a lattice
    int n_side = (int)ceil(pow(N_PARTICLES, 1.0/3.0));
    float spacing = BOX_SIZE / n_side;
    
    for (int i = 0; i < N_PARTICLES; i++) {
        int ix = i % n_side;
        int iy = (i / n_side) % n_side;
        int iz = i / (n_side * n_side);
        
        h_positions[i].x = ix * spacing + spacing * 0.5f;
        h_positions[i].y = iy * spacing + spacing * 0.5f;
        h_positions[i].z = iz * spacing + spacing * 0.5f;
        
        // Random velocities (Maxwell-Boltzmann at T=300K would be better)
        h_velocities[i].x = ((float)rand() / RAND_MAX - 0.5f) * 0.1f;
        h_velocities[i].y = ((float)rand() / RAND_MAX - 0.5f) * 0.1f;
        h_velocities[i].z = ((float)rand() / RAND_MAX - 0.5f) * 0.1f;
    }
}

float calculate_kinetic_energy(Vec3* h_velocities, float mass) {
    float ke = 0.0f;
    for (int i = 0; i < N_PARTICLES; i++) {
        float v2 = h_velocities[i].x * h_velocities[i].x +
                   h_velocities[i].y * h_velocities[i].y +
                   h_velocities[i].z * h_velocities[i].z;
        ke += 0.5f * mass * v2;
    }
    return ke;
}

int main() {
    printf("=== Molecular Dynamics with DMTCP Checkpointing ===\n");
    printf("Particles: %d\n", N_PARTICLES);
    printf("Box size: %.2f A\n", BOX_SIZE);
    printf("Time step: %.4f ps\n", DT);
    printf("Total steps: %d\n", N_STEPS);
    printf("Checkpoint frequency: %d steps\n\n", CHECKPOINT_FREQ);
    
    // Allocate host memory
    Vec3 *h_positions = (Vec3*)malloc(N_PARTICLES * sizeof(Vec3));
    Vec3 *h_velocities = (Vec3*)malloc(N_PARTICLES * sizeof(Vec3));
    Vec3 *h_forces = (Vec3*)malloc(N_PARTICLES * sizeof(Vec3));
    
    // Initialize system
    initialize_system(h_positions, h_velocities);
    
    // Allocate device memory
    Vec3 *d_positions, *d_velocities, *d_forces;
    cudaMalloc(&d_positions, N_PARTICLES * sizeof(Vec3));
    cudaMalloc(&d_velocities, N_PARTICLES * sizeof(Vec3));
    cudaMalloc(&d_forces, N_PARTICLES * sizeof(Vec3));
    
    // Copy to device
    cudaMemcpy(d_positions, h_positions, N_PARTICLES * sizeof(Vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_velocities, h_velocities, N_PARTICLES * sizeof(Vec3), cudaMemcpyHostToDevice);
    
    // Setup kernel launch parameters
    int threadsPerBlock = 256;
    int blocksPerGrid = (N_PARTICLES + threadsPerBlock - 1) / threadsPerBlock;
    
    float mass = 39.948f;  // Argon mass (amu)
    
    FILE *energy_file = fopen("energy.dat", "w");
    fprintf(energy_file, "# Step Time(ps) KineticEnergy(kcal/mol)\n");
    
    // Main MD loop
    clock_t start_time = clock();
    
    for (int step = 0; step < N_STEPS; step++) {
        // Calculate forces
        calculate_forces<<<blocksPerGrid, threadsPerBlock>>>(d_positions, d_forces, N_PARTICLES, BOX_SIZE);
        
        // Integrate positions (half-step velocities)
        integrate_positions<<<blocksPerGrid, threadsPerBlock>>>(d_positions, d_velocities, d_forces, 
                                                                N_PARTICLES, DT, BOX_SIZE, mass);
        
        // Calculate new forces
        calculate_forces<<<blocksPerGrid, threadsPerBlock>>>(d_positions, d_forces, N_PARTICLES, BOX_SIZE);
        
        // Complete velocity update
        update_velocities<<<blocksPerGrid, threadsPerBlock>>>(d_velocities, d_forces, 
                                                              N_PARTICLES, DT, mass);
        
        // Periodic output and checkpointing
        if (step % CHECKPOINT_FREQ == 0) {
            cudaMemcpy(h_velocities, d_velocities, N_PARTICLES * sizeof(Vec3), cudaMemcpyDeviceToHost);
            
            float ke = calculate_kinetic_energy(h_velocities, mass);
            float time = step * DT;
            
            fprintf(energy_file, "%d %.6f %.6f\n", step, time, ke);
            fflush(energy_file);
            
            printf("Step %d / %d (%.1f%%) - Time: %.3f ps - KE: %.2f kcal/mol\n", 
                   step, N_STEPS, 100.0*step/N_STEPS, time, ke);
            
            // DMTCP will handle checkpointing externally via dmtcp_command
            // The application just needs to be checkpoint-safe (no open network connections, etc.)
        }
    }
    
    clock_t end_time = clock();
    double elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    
    printf("\nSimulation completed in %.2f seconds\n", elapsed);
    printf("Performance: %.2f steps/second\n", N_STEPS / elapsed);
    
    fclose(energy_file);
    
    // Cleanup
    cudaFree(d_positions);
    cudaFree(d_velocities);
    cudaFree(d_forces);
    free(h_positions);
    free(h_velocities);
    free(h_forces);
    
    printf("\nResults saved to energy.dat\n");
    printf("To checkpoint: dmtcp_command --checkpoint\n");
    printf("To restart: dmtcp_restart ckpt_*.dmtcp\n");
    
    return 0;
}
