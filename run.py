from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
import os

# ====== GENERAL CONFIGURATION ======
cif_file = "1.cif"
ball_diameter = 0.3
ball_radius = ball_diameter / 2
time_step = 0.1
num_steps = 1000
temperature = 300  # Kelvin
ball_mass = 1.0  # eV·ps²/Å²
x0_variation = 6.0
num_launches = 600
disorder = 0.0
replication = (2, 2, 1)
num_runs = 1

# ====== CONSTANTS ======
kB = 8.617333262145e-5  # eV/K
thermal_velocity = np.sqrt(2 * kB * temperature / ball_mass)
gravity_acceleration = -0.0

# ====== ATOMIC RADII RELATIVE TO CARBON ======
atomic_radii = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    # ... add more elements as needed ...
}
carbon_radius = atomic_radii['C']

print(f"Average thermal velocity: {thermal_velocity:.4f} Å/ps")
output = open("output.txt", "w")

for config in range(1, num_runs + 1):
    atoms = read(cif_file)
    cell = atoms.get_cell()[:2, :2]
    inv_cell = np.linalg.inv(cell.T)
    atoms_replicated = atoms.repeat(replication)
    pin_positions = atoms_replicated.get_positions()[:, :2]
    symbols = atoms_replicated.get_chemical_symbols()

    pin_radii = [atomic_radii.get(s, carbon_radius) / carbon_radius * 0.6 for s in symbols]
    pin_radii = np.array(pin_radii)

    if disorder > 0.0:
        num_pins = len(pin_positions)
        num_remove = int(num_pins * disorder)
        remove_indices = np.random.choice(num_pins, size=num_remove, replace=False)
        pin_positions = np.delete(pin_positions, remove_indices, axis=0)
        pin_radii = np.delete(pin_radii, remove_indices, axis=0)

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

        ay = gravity_acceleration
        traj = [pos.copy()]
        collisions = []
        collision_angles = []

        for _ in range(num_steps):
            vel[1] += ay * time_step
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

    all_free_paths = []
    relaxation_times = []
    all_angles = []
    total_collisions = []

    for _ in range(num_launches):
        attempts = 0
        while True:
            x0 = x0_base + np.random.uniform(-x0_variation, x0_variation)
            result = simulate_collisions(x0)

            attempts += 1
            if attempts > 100:
                break

            if result is not None:
                traj, collisions, angles = result

                if len(collisions) == 0:
                    continue

                if len(collisions) >= 2:
                    for j in range(1, len(collisions)):
                        delta = collisions[j] - collisions[j - 1]
                        frac_delta = inv_cell @ delta
                        frac_delta = (frac_delta + np.array([Nx / 2, Ny / 2])) % [Nx, Ny] - np.array([Nx / 2, Ny / 2])
                        corrected_delta = cell.T @ frac_delta
                        d = np.linalg.norm(corrected_delta)

                        if d < 2 * (0.6 + ball_radius) * 5:
                            all_free_paths.append(d)
                            relax_time = d / thermal_velocity
                            relaxation_times.append(relax_time)

                all_angles.extend(angles)
                total_collisions.extend(collisions)
                break

    if all_free_paths:
        mean_path = np.mean(all_free_paths)
        mean_time = np.mean(relaxation_times)
        D = mean_path**2 / (4 * mean_time)

        print(f"Configuration {config}: λ = {mean_path:.3f} Å | τ = {mean_time:.3f} ps | D = {D:.5f} Å²/ps")
        output.write(f"Configuration {config}: λ = {mean_path:.3f} Å | τ = {mean_time:.3f} ps | D = {D:.5f} Å²/ps\n")

        plt.figure(figsize=(10, 10), dpi=300)
        plt.hist(all_free_paths, bins=15, color='skyblue', edgecolor='black')
        plt.xlabel("Distance between collisions (Å)", fontsize=18)
        plt.ylabel("Frequency", fontsize=18)
        plt.title(f"Free Path", fontsize=18)
        plt.xticks(fontsize=18)
        plt.yticks(np.arange(0, plt.ylim()[1]+1, 1000),fontsize=18)
        plt.grid(True)        
        plt.tight_layout()
        plt.savefig(f"histogram_paths_config_{config}.png")
        plt.close()

        plt.figure(figsize=(10, 6), dpi=300)
        plt.hist(relaxation_times, bins=15, color='salmon', edgecolor='black')
        plt.xlabel("Time between collisions (ps)", fontsize=14)
        plt.ylabel("Frequency", fontsize=14)
        plt.title(f"Histogram - Relaxation Time - Config {config}", fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"histogram_times_config_{config}.png")
        plt.close()

        plt.figure(figsize=(10, 6), dpi=300)
        plt.hist(all_angles, bins=18, range=(-180, 180), color='lightgreen', edgecolor='black')
        plt.xlabel("Collision angle (degrees)", fontsize=14)
        plt.ylabel("Frequency", fontsize=14)
        plt.title(f"Histogram - Collision Angle", fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"histogram_angles_config_{config}.png")
        plt.close()

        hist_angles, _ = np.histogram(all_angles, bins=18, range=(-180, 180), density=True)
        p_angles = hist_angles + 1e-12
        entropy = -np.sum(p_angles * np.log(p_angles))
        max_entropy = np.log(len(hist_angles))
        normalized_entropy = entropy / max_entropy

        fft_angles = np.fft.fft(hist_angles)
        fft_magnitudes = np.abs(fft_angles)
        dominant_order = np.argmax(fft_magnitudes[1:9]) + 1
        order_intensity = fft_magnitudes[dominant_order] / fft_magnitudes[0]

        output.write(f"    Angular entropy: {entropy:.4f} (normalized: {normalized_entropy:.4f})\n")
        output.write(f"    Dominant rotational order: {dominant_order}-fold (intensity: {order_intensity:.4f})\n")
        print(f"    Angular entropy: {entropy:.4f} (normalized: {normalized_entropy:.4f})")
        print(f"    Dominant rotational order: {dominant_order}-fold (intensity: {order_intensity:.4f})")

        if normalized_entropy > 0.8:
            structure = "amorphous (high angular disorder)"
        elif order_intensity > 0.1:
            structure = f"ordered with {dominant_order}-fold symmetry"
        elif order_intensity > 0.05:
            structure = f"quasicrystalline with weak {dominant_order}-fold symmetry"
        else:
            structure = "no significant angular order"

        output.write(f"    Structure: {structure}\n")
        print(f"    Structure: {structure}")

        plt.figure(figsize=(10, 10), dpi=300)
        modes = np.arange(len(fft_magnitudes))
        plt.stem(modes[1:10], fft_magnitudes[1:10], basefmt=" ", linefmt='C0-', markerfmt='C0o')
        plt.xlabel("Rotational mode (n-fold)", fontsize=18)
        plt.ylabel("Intensity (|FFT|)", fontsize=18)
        plt.title(f"Angular Order Spectrum", fontsize=18)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"fourier_spectrum_angular_config_{config}.png")
        plt.close()

        if total_collisions:
            collisions_xy = np.array(total_collisions)
            plt.figure(figsize=(8, 8), dpi=300)
            plt.scatter(collisions_xy[:, 0], collisions_xy[:, 1], s=20, color='blue', alpha=0.6, label='Collisions')
            plt.scatter(pin_positions[:, 0], pin_positions[:, 1], s=25, color='red', alpha=0.4, label='Pins')
            plt.xlabel("x (Å)", fontsize=18)
            plt.ylabel("y (Å)", fontsize=18)
            plt.title(f"Graphene", fontsize=16)
            plt.xticks(fontsize=18)
            plt.yticks(fontsize=18)
            plt.legend(fontsize=18)
            plt.axis("equal")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(f"collisions_2D_config_{config}.png")
            plt.close()

    else:
        print(f"Configuration {config}: No free path recorded.")
        output.write(f"Configuration {config}: No free path recorded.\n")

output.close()
print("Simulations completed. Results saved in 'output.txt' and PNG images.")
