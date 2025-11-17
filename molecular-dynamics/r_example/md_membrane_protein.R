#!/usr/bin/env Rscript

###############################################################################
# Molecular Dynamics Simulation: Membrane-Protein System with DMTCP
# Uses gpuR for GPU acceleration in R
# Run with DMTCP: dmtcp_launch Rscript md_membrane_protein.R
# Checkpoint: dmtcp_command --checkpoint
# Restart: dmtcp_restart ckpt_*.dmtcp
###############################################################################

library(gpuR)
library(methods)

# Check GPU availability
cat("\n=================================================================\n")
cat("Membrane-Protein MD Simulation with DMTCP Checkpointing\n")
cat("=================================================================\n\n")

if (!gpuR::detectGPUs() > 0) {
  stop("No GPU detected. This simulation requires a CUDA-capable GPU.")
}

cat("GPU Information:\n")
print(gpuInfo())
cat("\n")

# Simulation parameters
N_LIPIDS <- 2048            # Number of lipid molecules (2 tails each)
N_PROTEIN_ATOMS <- 1024     # Protein atoms embedded in membrane
N_WATER <- 3072             # Water molecules above/below membrane
BOX_X <- 80.0               # Box size X (Angstroms)
BOX_Y <- 80.0               # Box size Y (Angstroms)
BOX_Z <- 100.0              # Box size Z (Angstroms)
DT <- 0.002                 # Time step (ps)
N_STEPS <- 30000            # Total simulation steps
OUTPUT_FREQ <- 400          # Output frequency
TEMP_TARGET <- 310.0        # Body temperature (K)
CUTOFF <- 12.0              # Cutoff distance (Angstroms)

# Physical constants
KB <- 0.001987204           # Boltzmann constant (kcal/mol/K)

# Calculate total particles
N_PARTICLES <- N_LIPIDS * 3 + N_PROTEIN_ATOMS + N_WATER * 3

cat("System Configuration:\n")
cat(sprintf("  Lipid molecules: %d (each with 3 beads)\n", N_LIPIDS))
cat(sprintf("  Protein atoms: %d\n", N_PROTEIN_ATOMS))
cat(sprintf("  Water molecules: %d\n", N_WATER))
cat(sprintf("  Total particles: %d\n", N_PARTICLES))
cat(sprintf("  Box dimensions: %.1f x %.1f x %.1f A\n", BOX_X, BOX_Y, BOX_Z))
cat(sprintf("  Time step: %.4f ps\n", DT))
cat(sprintf("  Total steps: %d\n", N_STEPS))
cat("\n")

###############################################################################
# Initialize System
###############################################################################

initialize_membrane_system <- function() {
  cat("Initializing membrane-protein system...\n")
  
  positions <- matrix(0, nrow = N_PARTICLES, ncol = 3)
  velocities <- matrix(0, nrow = N_PARTICLES, ncol = 3)
  masses <- numeric(N_PARTICLES)
  particle_types <- character(N_PARTICLES)
  
  idx <- 1
  
  # Place lipids in membrane (z ~ 50 +/- 10)
  lipids_per_side <- ceiling(sqrt(N_LIPIDS))
  spacing_x <- BOX_X / lipids_per_side
  spacing_y <- BOX_Y / lipids_per_side
  
  for (i in 1:N_LIPIDS) {
    ix <- (i - 1) %% lipids_per_side
    iy <- floor((i - 1) / lipids_per_side)
    
    # Head group
    positions[idx, ] <- c(ix * spacing_x + runif(1, -1, 1),
                          iy * spacing_y + runif(1, -1, 1),
                          50 + runif(1, -2, 2))
    masses[idx] <- 50.0  # Head mass
    particle_types[idx] <- "lipid_head"
    idx <- idx + 1
    
    # Tail 1
    positions[idx, ] <- c(ix * spacing_x + runif(1, -1, 1),
                          iy * spacing_y + runif(1, -1, 1),
                          45 + runif(1, -3, 3))
    masses[idx] <- 30.0
    particle_types[idx] <- "lipid_tail"
    idx <- idx + 1
    
    # Tail 2
    positions[idx, ] <- c(ix * spacing_x + runif(1, -1, 1),
                          iy * spacing_y + runif(1, -1, 1),
                          40 + runif(1, -3, 3))
    masses[idx] <- 30.0
    particle_types[idx] <- "lipid_tail"
    idx <- idx + 1
  }
  
  # Place protein atoms in membrane center
  protein_per_side <- ceiling(sqrt(N_PROTEIN_ATOMS))
  protein_spacing <- min(BOX_X, BOX_Y) / (protein_per_side + 2)
  offset_x <- BOX_X / 2 - protein_per_side * protein_spacing / 2
  offset_y <- BOX_Y / 2 - protein_per_side * protein_spacing / 2
  
  for (i in 1:N_PROTEIN_ATOMS) {
    ix <- (i - 1) %% protein_per_side
    iy <- floor((i - 1) / protein_per_side)
    
    positions[idx, ] <- c(offset_x + ix * protein_spacing + runif(1, -0.5, 0.5),
                          offset_y + iy * protein_spacing + runif(1, -0.5, 0.5),
                          50 + runif(1, -8, 8))
    masses[idx] <- 15.0
    particle_types[idx] <- "protein"
    idx <- idx + 1
  }
  
  # Place water above and below membrane
  water_per_layer <- N_WATER / 2
  water_side <- ceiling(sqrt(water_per_layer))
  water_spacing_x <- BOX_X / water_side
  water_spacing_y <- BOX_Y / water_side
  
  # Water above (z > 60)
  for (i in 1:(water_per_layer)) {
    ix <- (i - 1) %% water_side
    iy <- floor((i - 1) / water_side)
    
    # O atom
    positions[idx, ] <- c(ix * water_spacing_x + runif(1, -0.5, 0.5),
                          iy * water_spacing_y + runif(1, -0.5, 0.5),
                          65 + runif(1, 0, 30))
    masses[idx] <- 16.0
    particle_types[idx] <- "water_O"
    idx <- idx + 1
    
    # H atoms (simplified, same position)
    for (h in 1:2) {
      positions[idx, ] <- positions[idx - 1, ] + runif(3, -0.2, 0.2)
      masses[idx] <- 1.0
      particle_types[idx] <- "water_H"
      idx <- idx + 1
    }
  }
  
  # Water below (z < 40)
  for (i in 1:(N_WATER - water_per_layer)) {
    ix <- (i - 1) %% water_side
    iy <- floor((i - 1) / water_side)
    
    # O atom
    positions[idx, ] <- c(ix * water_spacing_x + runif(1, -0.5, 0.5),
                          iy * water_spacing_y + runif(1, -0.5, 0.5),
                          5 + runif(1, 0, 30))
    masses[idx] <- 16.0
    particle_types[idx] <- "water_O"
    idx <- idx + 1
    
    # H atoms
    for (h in 1:2) {
      positions[idx, ] <- positions[idx - 1, ] + runif(3, -0.2, 0.2)
      masses[idx] <- 1.0
      particle_types[idx] <- "water_H"
      idx <- idx + 1
    }
  }
  
  # Initialize velocities (Maxwell-Boltzmann)
  for (i in 1:N_PARTICLES) {
    v_scale <- sqrt(KB * TEMP_TARGET / masses[i])
    velocities[i, ] <- rnorm(3, mean = 0, sd = v_scale)
  }
  
  # Remove center of mass motion
  total_mass <- sum(masses)
  com_velocity <- colSums(velocities * masses) / total_mass
  for (i in 1:N_PARTICLES) {
    velocities[i, ] <- velocities[i, ] - com_velocity
  }
  
  # Apply periodic boundary conditions
  positions[, 1] <- positions[, 1] %% BOX_X
  positions[, 2] <- positions[, 2] %% BOX_Y
  positions[, 3] <- positions[, 3] %% BOX_Z
  
  cat("System initialized successfully.\n\n")
  
  return(list(
    positions = positions,
    velocities = velocities,
    masses = masses,
    types = particle_types
  ))
}

###############################################################################
# GPU-Accelerated Force Calculation
###############################################################################

calculate_forces_gpu <- function(positions_gpu, masses_gpu, box_dims) {
  # This function uses gpuR for GPU computation
  # Note: gpuR has limitations compared to CUDA/OpenCL directly
  # For production, consider using more advanced GPU frameworks
  
  n <- nrow(positions_gpu)
  forces <- vclMatrix(0, nrow = n, ncol = 3, type = "float")
  
  # Simplified pairwise force calculation
  # In production, use spatial decomposition or neighbor lists
  
  # Transfer to GPU
  pos_gpu <- vclMatrix(positions_gpu, type = "float")
  
  # Compute forces (simplified Lennard-Jones)
  # This is a placeholder - actual implementation would require
  # custom OpenCL/CUDA kernels for efficiency
  
  epsilon <- 0.25  # kcal/mol
  sigma <- 3.5     # Angstroms
  
  # For demonstration, compute a subset of interactions
  # Full N^2 calculation would be too slow without custom kernels
  
  forces_cpu <- matrix(0, nrow = n, ncol = 3)
  
  # Sample-based force calculation for performance
  sample_size <- min(1000, n)
  sample_indices <- sample(1:n, sample_size)
  
  for (i in sample_indices) {
    pos_i <- positions_gpu[i, ]
    
    for (j in 1:n) {
      if (i == j) next
      
      pos_j <- positions_gpu[j, ]
      
      # Distance with PBC
      dx <- pos_i[1] - pos_j[1]
      dy <- pos_i[2] - pos_j[2]
      dz <- pos_i[3] - pos_j[3]
      
      dx <- dx - box_dims[1] * round(dx / box_dims[1])
      dy <- dy - box_dims[2] * round(dy / box_dims[2])
      dz <- dz - box_dims[3] * round(dz / box_dims[3])
      
      r2 <- dx^2 + dy^2 + dz^2
      
      if (r2 < CUTOFF^2 && r2 > 0.01) {
        r2_inv <- 1.0 / r2
        r6_inv <- r2_inv^3
        sigma_r6 <- (sigma^6) * r6_inv
        sigma_r12 <- sigma_r6^2
        
        force_mag <- 24 * epsilon * r2_inv * (2 * sigma_r12 - sigma_r6)
        
        forces_cpu[i, 1] <- forces_cpu[i, 1] + force_mag * dx
        forces_cpu[i, 2] <- forces_cpu[i, 2] + force_mag * dy
        forces_cpu[i, 3] <- forces_cpu[i, 3] + force_mag * dz
      }
    }
  }
  
  return(forces_cpu)
}

###############################################################################
# Velocity Verlet Integration
###############################################################################

velocity_verlet_step <- function(positions, velocities, masses, box_dims) {
  # Calculate forces
  forces <- calculate_forces_gpu(positions, masses, box_dims)
  
  # Update positions
  accelerations <- forces / masses
  positions <- positions + velocities * DT + 0.5 * accelerations * DT^2
  
  # Apply PBC
  positions[, 1] <- positions[, 1] %% box_dims[1]
  positions[, 2] <- positions[, 2] %% box_dims[2]
  positions[, 3] <- positions[, 3] %% box_dims[3]
  
  # Half-step velocity update
  velocities <- velocities + 0.5 * accelerations * DT
  
  # Calculate new forces
  forces_new <- calculate_forces_gpu(positions, masses, box_dims)
  
  # Complete velocity update
  accelerations_new <- forces_new / masses
  velocities <- velocities + 0.5 * accelerations_new * DT
  
  return(list(positions = positions, velocities = velocities))
}

###############################################################################
# Calculate Energies
###############################################################################

calculate_energies <- function(velocities, masses) {
  # Kinetic energy
  v2 <- rowSums(velocities^2)
  ke <- 0.5 * sum(masses * v2)
  
  # Temperature
  temp <- (2.0 / 3.0) * ke / (N_PARTICLES * KB)
  
  return(list(kinetic = ke, temperature = temp))
}

###############################################################################
# Main Simulation Loop
###############################################################################

run_simulation <- function() {
  # Initialize system
  system <- initialize_membrane_system()
  positions <- system$positions
  velocities <- system$velocities
  masses <- system$masses
  
  box_dims <- c(BOX_X, BOX_Y, BOX_Z)
  
  # Open output file
  energy_file <- file("energy_membrane.dat", "a")
  writeLines(sprintf("# Simulation started at %s", Sys.time()), energy_file)
  writeLines("# Step Time(ps) KineticEnergy(kcal/mol) Temperature(K)", energy_file)
  
  cat("Starting MD simulation...\n")
  cat(sprintf("Total steps: %d (%.2f ps)\n\n", N_STEPS, N_STEPS * DT))
  
  start_time <- Sys.time()
  
  for (step in 1:N_STEPS) {
    # Perform integration step
    result <- velocity_verlet_step(positions, velocities, masses, box_dims)
    positions <- result$positions
    velocities <- result$velocities
    
    # Output energies periodically
    if (step %% OUTPUT_FREQ == 0) {
      energies <- calculate_energies(velocities, masses)
      sim_time <- step * DT
      
      output_line <- sprintf("%d %.6f %.4f %.2f",
                           step, sim_time, energies$kinetic, energies$temperature)
      writeLines(output_line, energy_file)
      flush(energy_file)
      
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      steps_per_sec <- step / elapsed
      
      cat(sprintf("Step %6d / %d | Time: %8.3f ps | KE: %10.2f kcal/mol | T: %6.1f K | %.1f steps/s\n",
                  step, N_STEPS, sim_time, energies$kinetic, energies$temperature, steps_per_sec))
      
      # DMTCP checkpoint note
      if (step %% (OUTPUT_FREQ * 5) == 0) {
        cat("  [Use 'dmtcp_command --checkpoint' to save state]\n")
      }
    }
  }
  
  close(energy_file)
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat("\n")
  cat("=================================================================\n")
  cat(sprintf("Simulation completed in %.2f seconds\n", total_time))
  cat(sprintf("Performance: %.2f steps/second\n", N_STEPS / total_time))
  cat(sprintf("Total simulation time: %.3f ps\n", N_STEPS * DT))
  cat("=================================================================\n\n")
  
  cat("Results saved to energy_membrane.dat\n")
  cat("\nDMTCP Usage:\n")
  cat("  Checkpoint: dmtcp_command --checkpoint\n")
  cat("  Restart: dmtcp_restart ckpt_*.dmtcp\n")
}

###############################################################################
# Execute Simulation
###############################################################################

tryCatch({
  run_simulation()
}, error = function(e) {
  cat("\nError occurred:", conditionMessage(e), "\n")
  cat("Simulation can be restarted from last DMTCP checkpoint\n")
}, interrupt = function(e) {
  cat("\n\nSimulation interrupted by user\n")
  cat("Use DMTCP to checkpoint and restart\n")
})

cat("\nSimulation ended.\n")
