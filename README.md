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
pip install pyamtb
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

# Calculate band structure using configuration file
pyamtb calculate --config config.toml --poscar POSCAR

# Calculate distances between atoms
pyamtb distance --poscar POSCAR --element1 Mn --element2 N
```

### Configuration

The package uses TOML files for configuration. Here's an example configuration file (`config.toml`):

```toml
# Basic parameters
dimk = 3                    # Dimension of k-space (1, 2, or 3)
dimr = 3                    # Dimension of real space
nspin = 2                   # Number of spin components (1 or 2)
a0 = 1.0                    # Lattice constant scaling factor

# Band structure calculation
k_path = ["G", "X", "M", "G"]  # k-point path
num_k_points = 100             # Number of k-points
k_labels = ["Γ", "X", "M", "Γ"]  # k-point labels

# Hopping parameters
t0 = 1.0                      # Reference hopping strength
t0_distance = 2.0             # Reference distance
lambda_ = 1.0                 # Decay parameter
max_neighbors = 2             # Maximum number of neighbors to consider
max_distance = 10.0           # Maximum hopping distance
mindist = 0.1                 # Minimum hopping distance

# Onsite energy and magnetism
onsite_energy = [0.0, 0.0]    # Onsite energy for each atom
magnetic_moment = 1.0         # Magnetic moment
magnetic_order = "++--"       # Magnetic order pattern

# Output settings
output_filename = "band_structure"
output_format = "png"
savedir = "."

# Debug options
is_print_tb_model = false
is_print_tb_model_hop = false
is_check_flat_bands = true
adjust_degenerate_bands = true
energy_threshold = 0.001
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
  author = {Wang Dinghui},
  title = {PyAMTB: A Python package for tight-binding model calculations},
  year = {2024},
  url = {https://github.com/ooteki-teo/pyamtb.git}
}
``` 