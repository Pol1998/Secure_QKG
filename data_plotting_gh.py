import numpy as np
import matplotlib.pyplot as plt
import os

# === Set the directory where the data files are located ===
# Change this path to the folder where your data files are stored
data_dir = "/path/to/your/data"  # <-- Replace this with your actual path

# === Plot 1: L vs N ===

# Define the dataset labels and corresponding filenames
datasets = {
    "Noiseless": "data_noiseless.txt",
    "Low Noise": "data_low_noise.txt",
    "Estimated noise": "data_estimated.txt",
    "High Noise": "data_high_noise.txt",
    "Noisy": "data_noisy.txt"
}

# Colors for each curve
colors = {
    "Noiseless": "green",
    "Low Noise": "darkorange",
    "Estimated noise": "blue",
    "High Noise": "purple",
    "Noisy": "red"
}

# Initialize plot
fig, ax = plt.subplots(figsize=(9, 6))

# Plot each dataset
for label, filename in datasets.items():
    filepath = os.path.join(data_dir, filename)
    N, L = np.loadtxt(filepath, unpack=True, skiprows=1)
    ax.scatter(N, L, label=label, color=colors[label],
               edgecolor='black', linewidth=0.2, s=10)

# Axis settings
ax.set_xlabel("N")
ax.set_ylabel("L")
ax.set_xlim(0.5, 100.5)
ax.set_ylim(-1, 101)
ax.set_xticks([1] + list(range(20, 101, 20)))
ax.grid(True)

# Make axis lines thinner
for spine in ax.spines.values():
    spine.set_linewidth(0.5)

ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(data_dir, "plot_L_vs_N.pdf"), format='pdf')
plt.show()

# === Plot 2: Cost (χ) and QBER vs N ===

# Load cost data
N_cost, cost_low = np.loadtxt(os.path.join(data_dir, "data_cost_low_noise.txt"), unpack=True, skiprows=1)
_, cost_estimated = np.loadtxt(os.path.join(data_dir, "data_cost_estimated.txt"), unpack=True, skiprows=1)
_, cost_high = np.loadtxt(os.path.join(data_dir, "data_cost_high_noise.txt"), unpack=True, skiprows=1)

# Load QBER data
N_qber, qber_low = np.loadtxt(os.path.join(data_dir, "data_qber_low_noise.txt"), unpack=True, skiprows=1)
_, qber_estimated = np.loadtxt(os.path.join(data_dir, "data_qber_estimated.txt"), unpack=True, skiprows=1)
_, qber_high = np.loadtxt(os.path.join(data_dir, "data_qber_high_noise.txt"), unpack=True, skiprows=1)

# Setup dual-axis plot
fig, ax1 = plt.subplots(figsize=(9, 6))

# Plot χ (cost) on left y-axis
line1, = ax1.plot(N_cost, cost_low, label=r'$\chi$ Low noise', color='green', linestyle='-')
line2, = ax1.plot(N_cost, cost_estimated, label=r'$\chi$ Estimated noise', color='blue', linestyle='-')
line3, = ax1.plot(N_cost, cost_high, label=r'$\chi$ High noise', color='red', linestyle='-')

ax1.set_xlabel("N")
ax1.set_ylabel(r"$\chi$", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.set_xlim(1, 100)
ax1.set_yticks(np.arange(-0.5, 0.01, 0.1))
ax1.set_ylim(-0.52, 0.02)

# Thinner axis lines
for spine in ax1.spines.values():
    spine.set_linewidth(0.5)

# Plot QBER on right y-axis
ax2 = ax1.twinx()
line4, = ax2.plot(N_qber, qber_low, label='QBER Low noise', color='green', linestyle='--')
line5, = ax2.plot(N_qber, qber_estimated, label='QBER Estimated noise', color='blue', linestyle='--')
line6, = ax2.plot(N_qber, qber_high, label='QBER High noise', color='red', linestyle='--')

ax2.set_ylabel("QBER", color='black')
ax2.tick_params(axis='y', labelcolor='black')
ax2.set_ylim(-0.01, 0.5)
ax2.set_yticks(np.linspace(0, 0.5, 6))
ax1.spines['bottom'].set_zorder(10)

# Combine legends from both axes
lines = [line1, line2, line3, line4, line5, line6]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.7), ncol=2, frameon=True)

# Final layout
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.grid(True)
plt.savefig(os.path.join(data_dir, "plot_cost_qber_vs_N.pdf"), format='pdf')
plt.show()
