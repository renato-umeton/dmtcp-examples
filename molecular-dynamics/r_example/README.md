# R/gpuR Molecular Dynamics with DMTCP

## Description
This example simulates a membrane-protein system with surrounding water using R and gpuR for GPU acceleration. The simulation demonstrates DMTCP checkpointing with R scripts running GPU computations.

## Scientific Details
- **System**: Lipid bilayer membrane with embedded protein and water layers
  - 2048 lipid molecules (simplified 3-bead model)
  - 1024 protein atoms
  - 3072 water molecules
- **Potential**: Coarse-grained Lennard-Jones
- **Integration**: Velocity Verlet algorithm
- **GPU Acceleration**: gpuR (OpenCL-based GPU computing in R)
- **Ensemble**: NVE (microcanonical)

## Requirements
- R (version 4.0+)
- gpuR package
- OpenCL runtime
- DMTCP
- GPU with OpenCL support

## Installation

### Install R packages
```R
# Install from CRAN
install.packages("gpuR")

# Or install from source for latest version
install.packages("devtools")
devtools::install_github("cdeterman/gpuR")

# Verify installation
library(gpuR)
gpuInfo()
detectGPUs()
```

### Install OpenCL
```bash
# Ubuntu/Debian (NVIDIA)
sudo apt-get install nvidia-opencl-dev

# Ubuntu/Debian (AMD)
sudo apt-get install amdgpu-pro-opencl

# Check OpenCL
clinfo
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

## Usage

### Running with DMTCP
```bash
# Start the DMTCP coordinator (in separate terminal)
dmtcp_coordinator

# Launch the R simulation under DMTCP
dmtcp_launch Rscript md_membrane_protein.R

# Alternative: make script executable and run directly
chmod +x md_membrane_protein.R
dmtcp_launch ./md_membrane_protein.R
```

### Checkpointing
While the simulation is running, in another terminal:
```bash
# Create a checkpoint
dmtcp_command --checkpoint

# List checkpoint files
ls -lh ckpt_*.dmtcp
```

### Restarting from Checkpoint
```bash
# Ensure coordinator is running
dmtcp_coordinator

# Restart from checkpoint
dmtcp_restart ckpt_*.dmtcp

# The simulation will continue from where it was checkpointed
```

## Output
- `energy_membrane.dat`: Time series of kinetic energy and temperature
- Console output with simulation progress
- Checkpoint hints every 2000 steps

## Code Structure

### Main Components
1. **System Initialization** (`initialize_membrane_system`)
   - Places lipids in bilayer configuration
   - Embeds protein in membrane center
   - Adds water layers above and below
   - Initializes Maxwell-Boltzmann velocities

2. **Force Calculation** (`calculate_forces_gpu`)
   - GPU-accelerated pairwise interactions
   - Periodic boundary conditions
   - Simplified Lennard-Jones potential

3. **Integration** (`velocity_verlet_step`)
   - Velocity Verlet algorithm
   - Position and velocity updates

4. **Energy Analysis** (`calculate_energies`)
   - Kinetic energy calculation
   - Temperature estimation

## DMTCP Features Demonstrated
1. R process checkpointing
2. gpuR/OpenCL state preservation
3. GPU memory checkpointing
4. Seamless restart capability
5. No code modifications needed for checkpointing

## Performance Notes
- gpuR provides GPU acceleration through OpenCL
- Current implementation uses simplified force calculation for demonstration
- For production simulations:
  - Implement neighbor lists
  - Use optimized OpenCL kernels
  - Consider spatial decomposition
- Performance depends on GPU model and OpenCL implementation

## System Visualization
```
     +----------------------------------+
     |         Water Layer (top)        |  z > 60A
     +----------------------------------+
     |      Lipid Heads (upper)         |
     |  ~~~  Membrane  ~~~              |  z ~ 50A
     |    [Protein Embedded]            |
     |  ~~~  Membrane  ~~~              |
     |      Lipid Heads (lower)         |
     +----------------------------------+
     |        Water Layer (bottom)      |  z < 40A
     +----------------------------------+
```

## Verification Steps
1. Check GPU detection at startup
2. Verify particle counts match expected
3. Monitor energy conservation (NVE ensemble)
4. Temperature should fluctuate around 310K
5. Checkpoint and restart should preserve:
   - All particle positions and velocities
   - Simulation time counter
   - GPU state
   - Energy trajectories

## Advanced Usage

### Periodic checkpointing
```bash
# Checkpoint every 5 minutes automatically
dmtcp_launch --interval 300 Rscript md_membrane_protein.R
```

### Multiple checkpoints
```bash
# Keep last 3 checkpoints
dmtcp_launch --checkpoint-interval 600 --num-checkpoints 3 Rscript md_membrane_protein.R
```

### Restart from specific checkpoint
```bash
# List checkpoints by timestamp
ls -lt ckpt_*.dmtcp | head -5

# Restart from specific checkpoint
dmtcp_restart ckpt_<specific_timestamp>.dmtcp
```

## Limitations and Extensions

### Current Limitations
- Simplified force field (not production-quality)
- Sample-based force calculation for performance
- No proper thermostat/barostat
- Basic coarse-grained model

### Possible Extensions
1. Implement Martini force field for realistic membrane
2. Add Nosé-Hoover thermostat
3. Implement proper neighbor lists
4. Add trajectory output (XYZ, DCD formats)
5. Analyze membrane properties (thickness, area per lipid)
6. Study protein-lipid interactions

## Troubleshooting

### GPU not detected
```R
# Check GPU availability
library(gpuR)
gpuInfo()
detectGPUs()

# If no GPUs found, check OpenCL installation
system("clinfo")
```

### DMTCP issues
```bash
# Check DMTCP is running
dmtcp_command --status

# View coordinator log
cat /tmp/dmtcp-*.log
```

### Performance issues
- Reduce system size (N_LIPIDS, N_PROTEIN_ATOMS, N_WATER)
- Increase OUTPUT_FREQ to reduce I/O
- Check GPU memory usage
- Ensure OpenCL drivers are up to date

## References
- gpuR documentation: https://cran.r-project.org/package=gpuR
- DMTCP: http://dmtcp.sourceforge.net/
- Martini force field: http://cgmartini.nl/
- Membrane MD best practices: doi:10.1021/acs.jctc.5b00935
