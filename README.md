# DCF  
# Descriptor-Based Analysis of Atomic Structures

This repository contains two Python scripts that simulate particle collisions in atomic structures and extract structural descriptors, based on the methodology described in the associated scientific article.

---

## ⚠️ Important Note on CIF Files

These scripts work with **CIF (Crystallographic Information Files)** to define atomic structures.

> ⚠️ We strongly recommend using CIF files **exported from VESTA**, as they are typically compatible with ASE.  
> CIF files from other tools may cause issues during the replication step.

To improve compatibility, the repository includes a helper script called **`primitive.py`**, which reads a file named `input.cif` and exports a simplified version as `1.cif` — suitable for use by both `run.py` and `run_descriptor.py`.

---

## 📄 Scripts

### `run.py` — Single Structure Simulation

This script simulates particle trajectories through a periodic atomic structure defined in `1.cif`. It computes:

- Mean free path (λ)
- Relaxation time (τ)
- Diffusivity (D)
- Angular entropy and dominant rotational symmetry

#### 🔧 Inputs
- A CIF file named `1.cif` (placed in the root directory)
- A parameter file named `param.dat`

---

### ✅ `param.dat` (for `run.py`)

```ini
# Simulation parameters for run.py
cif_file = 1.cif
ball_diameter = 0.5
pin_radius = 0.5
time_step = 0.1
num_steps = 10000
temperature = 300
ball_mass = 1.0
x0_variation = 6.0
num_launches = 100
disorder = 0.0
replication = 2,2,1
num_runs = 1

output.txt: summary of transport and angular statistics

PNG and PDF figures:

histogram_paths_config_1.png

histogram_times_config_1.png

histogram_angles_config_1.png

fourier_spectrum_angular_config_1.png and .pdf

collisions_2D_config_1.png and .pdf

run_descriptor.py — Descriptor Generator for Multiple Structures
This script runs simulations on multiple CIF files and generates high-level descriptors for each structure.
It is suitable for dataset generation and use in machine learning or materials informatics.

🔧 Inputs
A folder containing .cif files (e.g., structures/)

A parameter file named param_desc.dat

✅ param_desc.dat (for run_descriptor.py)

# Simulation parameters for run_descriptor.py
structures_path = structures
output_descritor_csv = descritor.csv
output_fit_figure = fit
ball_diameter = 0.5
time_step = 0.1
num_steps = 10000
temperature = 300
ball_mass = 1.0
x0_variation = 6.0
num_launches = 100
disorder = 0.0
replication = 2,2,1
All .cif files in the folder specified by structures_path will be processed in numerical/alphabetical order.

📤 Outputs
descritor.csv: one row per structure with the following descriptors:

Mean, std, median, skewness, and kurtosis of free paths

Gaussian Mixture Model parameters (means and weights)

Angular entropy and symmetry intensities (1-fold to 9-fold)

Structure file name

GMM fit plots:

One figure per structure, named fit_<structure>.png

⚙️ Requirements
Install the required Python libraries with:

pip install -r requirements.txt
✅ requirements.txt
numpy
pandas
scipy
matplotlib
scikit-learn
ase

