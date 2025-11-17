# Quick Reference Guide

## Installation Commands

### C/CUDA Example
```bash
cd c_example
nvcc -o md_lj md_lennard_jones.cu -lm
```

### Python/CuPy Example
```bash
cd python_example
pip install cupy-cuda11x numpy
```

### R/gpuR Example
```R
# In R console
install.packages("gpuR")
```

## Running with DMTCP

### Start Coordinator (once, in separate terminal)
```bash
dmtcp_coordinator
```

### Launch Simulations

**C/CUDA:**
```bash
dmtcp_launch ./md_lj
```

**Python:**
```bash
dmtcp_launch python md_protein_water.py
```

**R:**
```bash
dmtcp_launch Rscript md_membrane_protein.R
```

## Checkpointing

### Create Checkpoint (while simulation running)
```bash
dmtcp_command --checkpoint
```

### Restart from Checkpoint
```bash
dmtcp_restart ckpt_*.dmtcp
```

## Key Features of Each Example

| Feature | C/CUDA | Python/CuPy | R/gpuR |
|---------|--------|-------------|--------|
| **System** | Argon fluid | Protein-water | Membrane-protein |
| **Particles** | 8,192 | 12,800 | 16,000+ |
| **GPU API** | CUDA | CUDA | OpenCL |
| **Best For** | HPC/Performance | Research | Statistical Analysis |
| **Code Lines** | ~230 | ~240 | ~470 |

## Expected Output

All examples produce:
- Energy data file (`energy*.dat`)
- Console progress reports
- Performance metrics (steps/second)
- Checkpoint instructions

## Verification Checklist

- [ ] GPU detected at startup
- [ ] Simulation runs without errors
- [ ] Energy file is being created
- [ ] Can checkpoint with `dmtcp_command --checkpoint`
- [ ] Can restart with `dmtcp_restart`
- [ ] Simulation continues from checkpoint

## Common Issues

**GPU not found:**
```bash
nvidia-smi          # Check NVIDIA GPU
clinfo              # Check OpenCL devices
```

**DMTCP coordinator not running:**
```bash
ps aux | grep dmtcp_coordinator
dmtcp_coordinator &  # Start in background
```

**Checkpoint files too large:**
```bash
dmtcp_launch --compress-checkpoint ./program
```

## Performance Expectations

On NVIDIA RTX 3080:
- C/CUDA: ~2,000 steps/second
- Python/CuPy: ~400 steps/second
- R/gpuR: ~100 steps/second

## Documentation Files

- `README.md` - Main overview
- `SUMMARY.md` - Detailed summary
- `VERIFICATION.md` - Verification checklist
- `c_example/README.md` - C-specific docs
- `python_example/README.md` - Python-specific docs
- `r_example/README.md` - R-specific docs

## Contact and Support

For issues with:
- DMTCP: https://github.com/dmtcp/dmtcp/issues
- CUDA: https://forums.developer.nvidia.com/
- CuPy: https://github.com/cupy/cupy/issues
- gpuR: https://github.com/cdeterman/gpuR/issues
