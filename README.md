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
num_steps = 10000

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
x0_variation = 8.0

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
- `structures/`: Folder where you place yours structures `.cif`.
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

# Number of time steps (Nsteps in paper)
num_steps = 10000

# Number of simulations (launches)
num_launches = 100

# Time step (ps)
time_step = 0.1

# Temperature (K)
temperature = 300

# Ball mass (eV·ps²/Å²)
ball_mass = 1.0

# Initial x position variation (Å) #variation for initial condition
x0_variation = 8.0

# Disorder level (0.0 to 1.0)
disorder = 0.0


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


