# C/CUDA Molecular Dynamics with DMTCP

## Description
This example simulates a Lennard-Jones fluid (argon-like particles) using CUDA for GPU acceleration. The simulation can be checkpointed and restarted using DMTCP.

## Scientific Details
- **System**: 8192 argon atoms in a periodic box
- **Potential**: Lennard-Jones 12-6
- **Integration**: Velocity Verlet algorithm
- **GPU Acceleration**: CUDA kernels for force calculation and integration

## Requirements
- NVIDIA CUDA Toolkit (nvcc compiler)
- DMTCP (Distributed MultiThreaded CheckPointing)
- CUDA-capable GPU

## Installation of DMTCP
```bash
# Ubuntu/Debian
sudo apt-get install dmtcp

# From source
git clone https://github.com/dmtcp/dmtcp.git
cd dmtcp
./configure
make
sudo make install
```

## Compilation
```bash
nvcc -o md_lj md_lennard_jones.cu -lm
```

## Usage

### Running with DMTCP
```bash
# Start the DMTCP coordinator (in a separate terminal)
dmtcp_coordinator

# Launch the simulation under DMTCP
dmtcp_launch ./md_lj
```

### Checkpointing
While the simulation is running, in another terminal:
```bash
# Create a checkpoint
dmtcp_command --checkpoint

# Check checkpoint files
ls ckpt_*.dmtcp
```

### Restarting from Checkpoint
```bash
# Make sure coordinator is running
dmtcp_coordinator

# Restart from checkpoint
dmtcp_restart ckpt_*.dmtcp
```

## Output
- `energy.dat`: Time series of kinetic energy
- Console output showing simulation progress

## DMTCP Features Demonstrated
1. Transparent checkpointing of GPU-accelerated code
2. CUDA memory state preservation
3. Automatic restart capability
4. No code modification needed for basic checkpointing

## Verification
The simulation should:
1. Run continuously for 100,000 steps
2. Output energy data every 1,000 steps
3. Be interruptible and restartable at any checkpoint
4. Continue seamlessly from the last checkpoint

## Notes
- DMTCP handles CUDA context and GPU memory automatically
- The application must be checkpoint-safe (no unclosed files, network connections, etc.)
- Checkpoints include full GPU state
