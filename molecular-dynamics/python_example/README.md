# Python/CuPy Molecular Dynamics with DMTCP

## Description
This example simulates a simplified protein in water system using CuPy for GPU acceleration. The simulation demonstrates DMTCP checkpointing with Python and CUDA via CuPy.

## Scientific Details
- **System**: 4096 water molecules + 512 protein atoms
- **Potential**: Simplified Lennard-Jones (for demonstration)
- **Integration**: Velocity Verlet algorithm
- **GPU Acceleration**: CuPy (CUDA-accelerated NumPy)
- **Ensemble**: NVE (constant number, volume, energy)

## Requirements
- Python 3.7+
- CuPy (GPU-accelerated NumPy)
- NumPy
- DMTCP
- CUDA-capable GPU

## Installation

### Install CuPy
```bash
# For CUDA 11.x
pip install cupy-cuda11x

# For CUDA 12.x
pip install cupy-cuda12x

# Check CuPy installation
python -c "import cupy as cp; print(cp.__version__); print(cp.cuda.is_available())"
```

### Install DMTCP
```bash
# Ubuntu/Debian
sudo apt-get install dmtcp

# Or from source
git clone https://github.com/dmtcp/dmtcp.git
cd dmtcp
./configure
make
sudo make install
```

### Install other dependencies
```bash
pip install numpy
```

## Usage

### Running with DMTCP
```bash
# Start the DMTCP coordinator (in a separate terminal)
dmtcp_coordinator

# Launch the simulation under DMTCP
dmtcp_launch python md_protein_water.py
```

### Checkpointing
While the simulation is running, in another terminal:
```bash
# Create a checkpoint
dmtcp_command --checkpoint

# Checkpoint files will be created
ls ckpt_*.dmtcp
```

### Restarting from Checkpoint
```bash
# Ensure coordinator is running
dmtcp_coordinator

# Restart from checkpoint
dmtcp_restart ckpt_*.dmtcp
```

## Output
- `energy_protein_water.dat`: Time series of kinetic, potential energy and temperature
- Console output showing simulation progress and performance

## DMTCP Features Demonstrated
1. Python process checkpointing
2. CuPy GPU array state preservation
3. CUDA context checkpointing
4. Automatic restart with full state restoration
5. No modification to core simulation code

## Code Structure
- `WaterProteinMD` class: Main simulation engine
  - `_initialize_system()`: Set up initial conditions
  - `_calculate_forces_gpu()`: GPU-accelerated force calculation
  - `_velocity_verlet_step()`: Integration step
  - `run_simulation()`: Main loop

## Performance Notes
- Uses GPU for all array operations via CuPy
- Pairwise distance calculations are vectorized on GPU
- For production simulations, implement neighbor lists or cell lists
- Current implementation is O(N^2) but GPU-accelerated

## Verification Steps
1. Simulation should run continuously
2. Energy output every 500 steps
3. Temperature should stabilize around target (300K with proper thermostat)
4. Checkpoint and restart should preserve:
   - Particle positions and velocities
   - GPU array states
   - Simulation step counter
   - Energy trajectories

## Advanced Usage

### Custom checkpoint intervals
Modify DMTCP to checkpoint at specific intervals:
```bash
dmtcp_launch --interval 300 python md_protein_water.py  # Every 5 minutes
```

### Multiple restarts
```bash
# List available checkpoints
ls -lt ckpt_*.dmtcp

# Restart from specific checkpoint
dmtcp_restart ckpt_<timestamp>.dmtcp
```

## Notes
- DMTCP automatically handles Python interpreter state
- CuPy GPU arrays are checkpointed transparently
- The simulation continues seamlessly from checkpoint
- For very large systems, consider using GPU-resident neighbor lists
- Real protein-water simulations would use proper force fields (AMBER, CHARMM, etc.)
