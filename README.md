# 💥 DCF — Descriptor-Based Analysis of Atomic Structures

This repository contains Python scripts to simulate particle collisions in atomic structures and extract structural descriptors, following the methodology described in the associated scientific publication.

---

## ⚠️ CIF File Compatibility

These scripts require **CIF (Crystallographic Information Files)** to define atomic structures.

> **Recommended:** Use CIF files **exported from VESTA**, which are typically compatible with ASE.  
> ⚠️ CIFs from other tools may cause errors during the replication step.

A helper script, `primitive.py`, is included to preprocess these files:
- Input: `input.cif`
- Output: `1.cif` (simplified and ready for simulation)

---

## 📜 Scripts Overview

### `run.py` — Single Structure Simulation

Simulates particle trajectories in a periodic atomic structure defined by `1.cif`. Computes:

- Mean Free Path (λ)
- Relaxation Time (τ)
- Diffusivity (D)
- Angular Entropy & Dominant Rotational Symmetry

#### 🔧 Required Inputs
- `1.cif` (structure file)
- `param.dat` (simulation parameters)

#### ✅ Example: `param.dat`

```ini
# Simulation parameters for run.py

# Name of the CIF file to load
cif_file = 1.cif

# Replication factors for the unit cell (Nx, Ny, Nz)
replication = 2,2,1

# Number of steps per simulation run (Nstep in paper)
num_steps = 1000 #(10000 in the paper)

# Number of launches per configuration
num_launches = 100

# Diameter of the ball (in angstroms)
ball_diameter = 0.5

# Radius of the pins (in angstroms), scaled relative to carbon
pin_radius = 0.5

# Time step for the simulation (in picoseconds)
time_step = 0.1

# Temperature of the system (in Kelvin)
temperature = 300

# Mass of the ball (in eV·ps²/Å²)
ball_mass = 1.0

# Range of initial horizontal position variation (in angstroms)
x0_variation = 6.0

# Fraction of pins to randomly remove (disorder level)
disorder = 0.0

# Number of independent simulation runs
num_runs = 1

```
#### **DCF** — Descriptor-Based Analysis of Atomic Structures

This project simulates the behavior of a "ball" launched into a 2D atomic structure obtained from .cif files.
It calculates properties such as the mean free path, relaxation time, and diffusivity. The simulation also
performs statistical analysis using Gaussian Mixture Models (GMM) and Fourier transforms of collision angles.

---

## 📁 Files used to obtain descriptors

- `run_descriptor.py`: Main simulation script..
- `param_desc.dat`: Configuration file with simulation parameters. 
- `structures/`: Directory where you place your `.cif` structures.  
  > **Note:** The script will collect all `.cif` files in this directory and, if they are named with leading numbers (e.g. `1.cif`, `2.cif`, `3.cif`), it will process them in ascending numeric order. Otherwise, it will fall back to the default lexicographical order of the file system. 
- `descritor.csv`: Output file containing calculated descriptors.
- `fit_<structure>.png`: Gaussian fit plots of free paths (output for each structure).



#### ✅ Example: `param_desc.dat`

```ini
# Simulation parameters for run_descriptor.py
# Directory containing .cif structures
structures_path = structures

# Output CSV for all descriptors
output_descritor_csv = descritor.csv

# Output PNG for Gaussian mixture fit (base name, structure name will be appended)
output_fit_figure = fit

# Unit cell replication in x, y, z
replication = 2,2,1

# Ball diameter in angstroms
ball_diameter = 0.5

# Pin radius in angstroms
pin_radius = 0.5

# Number of time steps (Nsteps in paper)
num_steps = 1000 #(10000 in the paper)

# Number of simulations (launches)
num_launches = 100

# Time step (ps)
time_step = 0.1

# Temperature (K)
temperature = 300

# Ball mass (eV·ps²/Å²)
ball_mass = 1.0

# Initial x position variation (Å) #variation for initial condition
x0_variation = 6.0

# Disorder level (0.0 to 1.0)
disorder = 0.0

```

 # ====== Atomic radii ======
atomic_radii = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
}

This dictionary defines the atomic radii for selected elements, relative 
to carbon (C = 0.76). These values are used to estimate the size of atoms in the structure
and determine collision detection during the simulation. Before running the simulation,
make sure that all atomic species in your .cif files are listed in the atomic_radii dictionary.
If a given element is not included, the simulation will default to using the carbon radius.


## ⚙️ Requirements

Install the required packages using:

```bash
pip install -r requirements.txt

Contents of requirements.txt:

numpy
pandas
scipy
matplotlib
scikit-learn
ase


