import numpy as np
import pandas as pd
from scipy.stats import kurtosis, skew
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
from ase.io import read

# ====== CONFIGURAÇÕES ======
cif_file = "1.cif"
output_descritor_csv = "descritor.csv"
output_fit_figure = "fit.png"

ball_diameter = 0.3
ball_radius = ball_diameter / 2
time_step = 0.1
num_steps = 1000
temperature = 300
ball_mass = 1.0
x0_variation = 6.0
num_launches = 600
disorder = 0.0
replication = (2, 2, 1)

kB = 8.617333262145e-5
thermal_velocity = np.sqrt(2 * kB * temperature / ball_mass)
gravity_acceleration = -0.0

atomic_radii = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
}
carbon_radius = atomic_radii['C']

# ====== PREPARAÇÃO ======
atoms = read(cif_file)
cell = atoms.get_cell()[:2, :2]
inv_cell = np.linalg.inv(cell.T)
atoms_replicated = atoms.repeat(replication)
pin_positions = atoms_replicated.get_positions()[:, :2]
symbols = atoms_replicated.get_chemical_symbols()
pin_radii = [atomic_radii.get(s, carbon_radius) / carbon_radius * 0.6 for s in symbols]
pin_radii = np.array(pin_radii)

Nx, Ny = replication[:2]
total_cell = cell.T @ np.array([Nx, Ny])
x0_base = (cell.T @ np.array([Nx / 2, Ny / 2]))[0]
y_top = np.max(pin_positions[:, 1])
y_bottom = np.min(pin_positions[:, 1])

def apply_periodic_boundary(pos):
    frac = inv_cell @ pos
    frac = frac % [Nx, Ny]
    return (cell.T @ frac)

def simulate_collisions(x0):
    y_initial = y_top + 2 * np.random.uniform(-x0_variation, x0_variation)
    pos = np.array([x0, y_initial])
    vx0 = np.random.uniform(-0.2, 0.2) * thermal_velocity
    vy0 = -np.sqrt(thermal_velocity**2 - vx0**2)
    vel = np.array([vx0, vy0])

    traj = [pos.copy()]
    collisions = []
    collision_angles = []

    for _ in range(num_steps):
        vel[1] += gravity_acceleration * time_step
        pos += vel * time_step
        pos = apply_periodic_boundary(pos)

        if pos[1] >= y_initial:
            return None

        for i, pin in enumerate(pin_positions):
            r_vec = pos - pin
            dist = np.linalg.norm(r_vec)
            total_radius = pin_radii[i] + ball_radius
            if dist < total_radius:
                normal = r_vec / dist
                angle = np.degrees(np.arctan2(normal[1], normal[0]))
                collision_angles.append(angle)

                vel -= 2 * np.dot(vel, normal) * normal
                pos = pin + normal * total_radius
                collisions.append(pos.copy())
                break

        traj.append(pos.copy())
        if pos[1] < y_bottom - 1.0:
            break

    return np.array(traj), np.array(collisions), np.array(collision_angles)

# ====== SIMULAÇÃO ======
all_free_paths = []
relaxation_times = []
all_angles = []
total_collisions = []

for _ in range(num_launches):
    for _ in range(100):
        x0 = x0_base + np.random.uniform(-x0_variation, x0_variation)
        result = simulate_collisions(x0)
        if result is not None:
            traj, collisions, angles = result
            if len(collisions) >= 2:
                for j in range(1, len(collisions)):
                    delta = collisions[j] - collisions[j - 1]
                    frac_delta = inv_cell @ delta
                    frac_delta = (frac_delta + np.array([Nx / 2, Ny / 2])) % [Nx, Ny] - np.array([Nx / 2, Ny / 2])
                    corrected_delta = cell.T @ frac_delta
                    d = np.linalg.norm(corrected_delta)
                    if d < 2 * (0.6 + ball_radius) * 5:
                        all_free_paths.append(d)
                        relaxation_times.append(d / thermal_velocity)
            all_angles.extend(angles)
            total_collisions.extend(collisions)
            break

# ====== DESCRITOR ======
descritor = {}

if all_free_paths:
    paths = np.array(all_free_paths)
    times = np.array(relaxation_times)
    mean_path = paths.mean()
    mean_time = times.mean()
    D = mean_path**2 / (4 * mean_time)

    descritor["mean_free_path"] = mean_path
    descritor["mean_relax_time"] = mean_time
    descritor["diffusivity"] = D
    descritor["std_path"] = paths.std()
    descritor["median_path"] = np.median(paths)
    descritor["skew_path"] = skew(paths)
    descritor["kurtosis_path"] = kurtosis(paths)

    # ====== AJUSTE GAUSSIANA BIMODAL ======
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(paths.reshape(-1, 1))
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_).flatten()
    weights = gmm.weights_

    descritor["gauss1_mean"] = means[0]
    descritor["gauss2_mean"] = means[1]
    descritor["gauss1_weight"] = weights[0]
    descritor["gauss2_weight"] = weights[1]

    # ====== GERA FIGURA DO AJUSTE FIT.PNG ======
    x = np.linspace(paths.min(), paths.max(), 1000)
    pdf1 = weights[0] * np.exp(-0.5 * ((x - means[0]) / stds[0])**2) / (stds[0] * np.sqrt(2 * np.pi))
    pdf2 = weights[1] * np.exp(-0.5 * ((x - means[1]) / stds[1])**2) / (stds[1] * np.sqrt(2 * np.pi))
    pdf_total = pdf1 + pdf2

    plt.figure(figsize=(10, 6), dpi=300)
    plt.hist(paths, bins=30, density=True, alpha=0.5, label="Histogram")
    plt.plot(x, pdf1, label="Gaussian 1", linestyle='--')
    plt.plot(x, pdf2, label="Gaussian 2", linestyle='--')
    plt.plot(x, pdf_total, label="GMM Total", linewidth=2)
    plt.xlabel("Free path (Å)", fontsize=14)
    plt.ylabel("Probability density", fontsize=14)
    plt.title("Gaussian Mixture Fit to Free Paths", fontsize=16)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(output_fit_figure)
    plt.close()

    # ====== FFT DOS ÂNGULOS ======
    hist_angles, _ = np.histogram(all_angles, bins=18, range=(-180, 180), density=True)
    p_angles = hist_angles + 1e-12
    entropy = -np.sum(p_angles * np.log(p_angles))
    max_entropy = np.log(len(hist_angles))
    normalized_entropy = entropy / max_entropy

    descritor["angular_entropy"] = entropy
    descritor["angular_entropy_norm"] = normalized_entropy

    fft_angles = np.fft.fft(hist_angles)
    fft_magnitudes = np.abs(fft_angles)

    for n in range(1, 10):
        descritor[f"{n}_fold_intensity"] = fft_magnitudes[n] / fft_magnitudes[0]

# ====== SALVA CSV ======
df_descritor = pd.DataFrame([descritor])
df_descritor.to_csv(output_descritor_csv, index=False)
print(f"Descritor salvo em: {output_descritor_csv}")
print(f"Figura salva em: {output_fit_figure}")
