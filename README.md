# PyAMTB

A Python package for tight-binding model calculations for altermagnets.

## Introduction

PyAMTB (Python Altermagnet Tight Binding) is built on top of the PythTB package, providing specialized tight-binding model calculations for altermagnets. It extends PythTB's capabilities by adding direct support for POSCAR structure files and altermagnet-specific features.

## Features

- Tight-binding model calculations for altermagnets
- Support for various lattice structures
- Band structure calculations
- Easy configuration through TOML files
- Command-line interface for quick calculations

## Installation

### From PyPI

```bash
pip install pyamtb --upgrade
```

### From source

```bash
git clone https://github.com/ooteki-teo/pyamtb.git
cd pyamtb
pip install -e .
```

## Usage

### Command Line Interface

The package provides a command-line interface for easy calculations:

```bash
# Show help and available commands
pyamtb --help

# Calculate distances between atoms
pyamtb distance --poscar POSCAR --element1 Mn --element2 N
# check the ditance, determin max_distance

# create a template.toml file
pyamtb template 

# Calculate band structure using configuration file
pyamtb calculate --config config.toml --poscar POSCAR

```

### Configuration

The package uses TOML files for configuration. Here's an example configuration file (`tbparas.toml`):

```toml
# Basic parameters
poscar="Mn2N.vasp"
dimk = 2                    # **important** Dimension of k-space (1, 2, or 3) 
dimr = 3                    # Dimension of real space, do not change this
nspin = 2                   # Number of spin components (1 or 2)
a0 = 1.0                    # Lattice constant scaling factor

# Band structure calculation
k_path = [
[0.0,0.0],
[0.5,0.0],
[0.5,0.5],
[0.0,0.5],
[0.0,0.0]
]  # k-point path
nkpt = 100             # Number of k-points
k_labels = ["Γ", "X", "M", "Y", "Γ"]  # k-point labels

# Hopping parameters
t0 = 1.0                      # Reference hopping strength
t0_distance = 2.5             # **important** Reference distance 
hopping_decay = 1.0                 # Decay parameter 
max_neighbors = 2             # Maximum number of neighbors to consider 
max_distance = 3.5           # **important** Maximum hopping distance 
min_distance = 0.1                 # Minimum hopping distance

# Onsite energy and magnetism
onsite_energy = [0.0, 0.0, 0.0]    # Onsite energy for each atom
magnetic_moment = 0.1         # Magnetic moment
magnetic_order = "+-0"       # **important** Magnetic order pattern

# Output settings
output_filename = "band_structure"   # **important**
output_format = "png"
savedir = "."

# some other options
is_print_tb_model = false
is_print_tb_model_hop = false
is_check_flat_bands = true
is_black_degenerate_bands = true  # plot the degenerate band in black, otherwise in blue/red for spin polarized
energy_threshold = 0.00001
```

### Python API

You can also use the package in your Python code:

```python
from pyamtb import Parameters, calculate_band_structure, create_pythtb_model

# Load parameters from TOML file
params = Parameters("config.toml")

# Create tight-binding model
model = create_pythtb_model("POSCAR", params)

# Calculate band structure
calculate_band_structure(model, params)
```

## Development

### Running Tests

```bash
pip install -e ".[dev]"
pytest
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```bibtex
@software{pyamtb,
  author = {Dinghui Wang, Junting Zhang, Yu Xie},
  title = {PyAMTB: A Python package for tight-binding model calculations},
  year = {2025},
  url = {https://github.com/ooteki-teo/pyamtb.git}
}
``` 
