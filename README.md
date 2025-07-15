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

# Simulador de Difusividade com Lançamentos de Bola

Este projeto simula o comportamento de uma "bola" lançada em uma estrutura atômica 2D obtida de arquivos `.cif`, calculando propriedades como caminho livre médio, tempo de relaxação, difusividade, e realizando análise estatística com Gaussian Mixture Models (GMM) e transformadas de Fourier dos ângulos de colisão.

---

## 📁 Estrutura do Projeto

- `main.py`: Script principal da simulação.
- `param_desc.dat`: Arquivo com os parâmetros da simulação.
- `structures/`: Pasta onde você coloca arquivos `.cif`.
- `descritor.csv`: Saída com os descritores calculados.
- `fit_<estrutura>.png`: Ajustes Gaussianos para os caminhos livres.
- `examples/`: Exemplos de saída.

---

## ⚙️ Requisitos

Instale os pacotes com:

```bash
pip install -r requirements.txt

