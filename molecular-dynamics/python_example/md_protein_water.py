"""
Molecular Dynamics Simulation: Protein in Water with DMTCP Checkpointing
Uses CuPy for GPU acceleration (CUDA-based NumPy)
Run with DMTCP: dmtcp_launch python md_protein_water.py
Checkpoint: dmtcp_command --checkpoint
Restart: dmtcp_restart ckpt_*.dmtcp
"""

import numpy as np
import cupy as cp
import time
import pickle
import os
from datetime import datetime

# Simulation parameters
N_WATER = 4096          # Number of water molecules
N_PROTEIN_ATOMS = 512   # Simplified protein atoms
BOX_SIZE = 60.0         # Simulation box (Angstroms)
DT = 0.002              # Time step (ps)
N_STEPS = 50000         # Total steps
CHECKPOINT_FREQ = 500   # Steps between energy outputs
TEMP_TARGET = 300.0     # Target temperature (K)
CUTOFF = 10.0           # Cutoff distance (Angstroms)

# Physical constants
KB = 0.001987204  # Boltzmann constant (kcal/mol/K)

class WaterProteinMD:
    """Molecular dynamics simulation of protein in water using GPU"""
    
    def __init__(self, n_water, n_protein, box_size):
        self.n_water = n_water
        self.n_protein = n_protein
        self.n_particles = n_water * 3 + n_protein  # Water = 3 atoms each
        self.box_size = box_size
        self.step = 0
        
        print(f"Initializing MD system:")
        print(f"  Water molecules: {n_water}")
        print(f"  Protein atoms: {n_protein}")
        print(f"  Total particles: {self.n_particles}")
        print(f"  Box size: {box_size} A")
        
        # Initialize on GPU
        self._initialize_system()
        
    def _initialize_system(self):
        """Initialize particle positions, velocities, and masses on GPU"""
        
        # Generate initial positions (lattice + random perturbation)
        n_side = int(np.ceil(self.n_particles ** (1/3)))
        spacing = self.box_size / n_side
        
        positions = []
        for i in range(self.n_particles):
            ix = i % n_side
            iy = (i // n_side) % n_side
            iz = i // (n_side * n_side)
            
            x = ix * spacing + np.random.uniform(-0.1, 0.1)
            y = iy * spacing + np.random.uniform(-0.1, 0.1)
            z = iz * spacing + np.random.uniform(-0.1, 0.1)
            
            positions.append([x % self.box_size, 
                            y % self.box_size, 
                            z % self.box_size])
        
        self.positions = cp.array(positions, dtype=cp.float32)
        
        # Initialize velocities (Maxwell-Boltzmann distribution)
        velocities = np.random.randn(self.n_particles, 3).astype(np.float32)
        velocities *= np.sqrt(KB * TEMP_TARGET / 18.0)  # Approximate mass
        
        # Remove center of mass motion
        velocities -= np.mean(velocities, axis=0)
        
        self.velocities = cp.array(velocities, dtype=cp.float32)
        
        # Assign masses (simplified: water O=16, H=1, protein=12 amu average)
        masses = []
        for i in range(self.n_water):
            masses.extend([16.0, 1.0, 1.0])  # O, H, H
        masses.extend([12.0] * self.n_protein)  # Protein atoms
        
        self.masses = cp.array(masses, dtype=cp.float32).reshape(-1, 1)
        
        # Initialize forces array
        self.forces = cp.zeros((self.n_particles, 3), dtype=cp.float32)
        
    def _apply_pbc(self, dr):
        """Apply periodic boundary conditions to distance vector"""
        return dr - self.box_size * cp.round(dr / self.box_size)
    
    def _calculate_forces_gpu(self):
        """Calculate forces using GPU-accelerated pairwise interactions"""
        
        self.forces.fill(0)
        
        # Simplified Lennard-Jones potential for all particles
        # In reality, water would use specialized potentials (TIP3P, SPC/E)
        
        # Compute all pairwise distances (memory-intensive but GPU-efficient)
        # For production, use cell lists or neighbor lists
        
        # Expand positions for broadcasting
        pos_i = self.positions[:, cp.newaxis, :]  # (N, 1, 3)
        pos_j = self.positions[cp.newaxis, :, :]  # (1, N, 3)
        
        # Distance vectors with PBC
        dr = pos_i - pos_j
        dr = self._apply_pbc(dr)
        
        # Distance squared
        r2 = cp.sum(dr**2, axis=2)
        
        # Avoid self-interaction and apply cutoff
        mask = (r2 > 0.01) & (r2 < CUTOFF**2)
        
        # Lennard-Jones parameters (simplified)
        epsilon = 0.2  # kcal/mol
        sigma = 3.0    # Angstroms
        
        # Calculate LJ force
        r2_inv = cp.where(mask, 1.0 / r2, 0)
        r6_inv = r2_inv**3
        sigma_r6 = (sigma**6) * r6_inv
        sigma_r12 = sigma_r6**2
        
        # Force magnitude: 24*epsilon/r^2 * (2*sigma^12/r^12 - sigma^6/r^6)
        force_mag = cp.where(mask, 
                            24 * epsilon * r2_inv * (2 * sigma_r12 - sigma_r6),
                            0)
        
        # Force vectors
        force_vectors = dr * force_mag[:, :, cp.newaxis]
        
        # Sum forces on each particle
        self.forces = cp.sum(force_vectors, axis=1)
        
    def _velocity_verlet_step(self):
        """Perform one velocity Verlet integration step"""
        
        # Update positions: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
        acceleration = self.forces / self.masses
        self.positions += self.velocities * DT + 0.5 * acceleration * DT**2
        
        # Apply periodic boundary conditions
        self.positions = self.positions % self.box_size
        
        # Half-step velocity update
        self.velocities += 0.5 * acceleration * DT
        
        # Calculate new forces
        self._calculate_forces_gpu()
        
        # Complete velocity update
        acceleration = self.forces / self.masses
        self.velocities += 0.5 * acceleration * DT
        
        self.step += 1
        
    def _calculate_energies(self):
        """Calculate kinetic and potential energies"""
        
        # Kinetic energy
        v2 = cp.sum(self.velocities**2, axis=1)
        ke = 0.5 * cp.sum(self.masses.flatten() * v2)
        
        # Simplified potential energy (just for monitoring)
        pos_i = self.positions[:, cp.newaxis, :]
        pos_j = self.positions[cp.newaxis, :, :]
        dr = self._apply_pbc(pos_i - pos_j)
        r2 = cp.sum(dr**2, axis=2)
        
        mask = (r2 > 0.01) & (r2 < CUTOFF**2)
        r2_inv = cp.where(mask, 1.0 / r2, 0)
        r6_inv = r2_inv**3
        
        epsilon = 0.2
        sigma = 3.0
        sigma_r6 = (sigma**6) * r6_inv
        sigma_r12 = sigma_r6**2
        
        pe = 4 * epsilon * cp.sum(sigma_r12 - sigma_r6) / 2  # Divide by 2 for double counting
        
        # Temperature
        temp = (2.0 / 3.0) * ke / (self.n_particles * KB)
        
        return float(ke.get()), float(pe.get()), float(temp.get())
    
    def run_simulation(self, n_steps):
        """Run the MD simulation for n_steps"""
        
        energy_file = open('energy_protein_water.dat', 'a')
        energy_file.write(f"# Restarted at step {self.step}\n")
        energy_file.write("# Step Time(ps) KE(kcal/mol) PE(kcal/mol) TotalE Temperature(K)\n")
        
        print(f"\nStarting simulation from step {self.step}")
        print("=" * 70)
        
        start_time = time.time()
        last_time = start_time
        
        for i in range(n_steps):
            self._velocity_verlet_step()
            
            if self.step % CHECKPOINT_FREQ == 0:
                ke, pe, temp = self._calculate_energies()
                total_e = ke + pe
                sim_time = self.step * DT
                
                energy_file.write(f"{self.step} {sim_time:.6f} {ke:.4f} {pe:.4f} "
                                f"{total_e:.4f} {temp:.2f}\n")
                energy_file.flush()
                
                current_time = time.time()
                elapsed = current_time - last_time
                steps_per_sec = CHECKPOINT_FREQ / elapsed if elapsed > 0 else 0
                
                print(f"Step {self.step:6d} | Time: {sim_time:8.3f} ps | "
                      f"KE: {ke:10.2f} | PE: {pe:10.2f} | "
                      f"T: {temp:6.1f} K | {steps_per_sec:.1f} steps/s")
                
                last_time = current_time
                
                # DMTCP handles checkpointing externally
                # Application state in GPU memory will be preserved
        
        energy_file.close()
        
        total_time = time.time() - start_time
        print("\n" + "=" * 70)
        print(f"Simulation completed {n_steps} steps in {total_time:.2f} seconds")
        print(f"Performance: {n_steps/total_time:.2f} steps/second")
        print(f"Total simulation time: {self.step * DT:.3f} ps")

def main():
    print("=" * 70)
    print("Protein in Water MD Simulation with DMTCP Checkpointing")
    print("=" * 70)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"GPU Device: {cp.cuda.Device()}")
    print(f"CUDA Available: {cp.cuda.is_available()}")
    print()
    
    # Check if we're restarting (state file would exist from manual save)
    # Note: DMTCP handles the actual process/memory checkpoint
    # This is just for demonstration of application-level state
    
    md = WaterProteinMD(N_WATER, N_PROTEIN_ATOMS, BOX_SIZE)
    
    try:
        md.run_simulation(N_STEPS)
    except KeyboardInterrupt:
        print("\n\nSimulation interrupted by user")
        print("Use DMTCP to checkpoint and restart")
    
    print("\nTo checkpoint this simulation:")
    print("  In another terminal: dmtcp_command --checkpoint")
    print("\nTo restart from checkpoint:")
    print("  dmtcp_restart ckpt_*.dmtcp")

if __name__ == "__main__":
    main()
