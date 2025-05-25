# DCF  
# Descriptor-Based Analysis of Atomic Structures

This repository contains two Python scripts that simulate particle collisions in atomic structures and extract structural descriptors, based on the methodology described in the associated scientific article.

---

## ⚠️ Important Note on CIF Files

These scripts expect a CIF file named **`1.cif`** as input. We strongly recommend using CIF files **exported from VESTA**, as they are generally compatible.

> ⚠️ CIF files generated from other tools may cause issues during the replication step.

To improve compatibility, we provide a helper script called **`primitive.py`** (included in this repository).  
This script reads a file named **`input.cif`** and generates a simplified version saved as **`1.cif`**, which can then be used directly by both `run.py` and `run_descriptor.py`.

---

## 📄 Scripts

### `run.py`

This script simulates particles colliding with a replicated atomic structure loaded from `1.cif`. It calculates physical properties such as:

- Mean free path (λ)
- Relaxation time (τ)
- Diffusivity (D)
- Angular entropy and rotational order (via FFT)

**Input:**
- A CIF file named `1.cif` (preprocessed if necessary)

**Outputs:**
- `output.txt`: summary of the simulation
- PNG images:
  - Histograms of free paths, relaxation times, and collision angles
  - Angular order Fourier spectrum
  - 2D map of collision events

**Configuration:**
Simulation parameters (temperature, time steps, disorder, etc.) are located at the **top of the script** for easy modification.

---

### `run_descriptor.py`

This script performs similar simulations but focuses on generating **a numerical descriptor** that can be used for data analysis or machine learning.

**Outputs:**
- `descritor.csv`: containing statistical and physical metrics:
  - Mean, std, median, skewness, kurtosis of free paths
  - GMM fit parameters (means, weights)
  - Angular entropy and symmetry intensities
- `fit.png`: Gaussian Mixture Model fit over the histogram

---

## ⚙️ Requirements

Install dependencies via:

```bash
pip install -r requirements.txt
