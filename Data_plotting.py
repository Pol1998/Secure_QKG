##################
#Packages required
##################
import matplotlib.pyplot as plt
import matplotlib
import os

#######
#Plot 1
#######
#Note: Load saved data

# Labels and colors
labels = ["noiseless", "low-noise", "large-noise"]
colors_left = ["navy", "forestgreen", "darkred"]    # for L, S
colors_right = ["orange", "purple", "brown"]        # for Φ

# Plot setup
fig, ax_left = plt.subplots(figsize=(10, 6))
ax_right = ax_left.twinx()

# Plot each dataset
for i, filename in enumerate(files):
    data = load_data(filename)
    N = list(range(1, len(data) + 1))

    a_vals = [d[0] for d in data]
    b_vals = [d[1] for d in data]
    c_vals = [d[2] for d in data]

    ax_left.scatter(N, a_vals, color=colors_left[i], marker='o', s=20, label=rf"{labels[i]} $L$")
    ax_left.scatter(N, b_vals, color=colors_left[i], marker='x', s=20, label=rf"{labels[i]} $S$")
    ax_right.plot(N, c_vals, color=colors_right[i], linestyle='-', label=rf"{labels[i]} $\chi$")

# Axis labels
ax_left.set_xlabel(r"$N$")
ax_left.set_ylabel(r"$L$, $S$")
ax_right.set_ylabel(r"$\chi$")

# Show x-ticks only at 1, 25, 50, 75, 100 (if available)
max_N = len(data)
tick_positions = [i for i in [1, 25, 50, 75, 100] if i <= max_N]
ax_left.set_xticks(tick_positions)

# Combine legends from both axes and display it
lines_left, labels_left = ax_left.get_legend_handles_labels()
lines_right, labels_right = ax_right.get_legend_handles_labels()
combined_lines = lines_left + lines_right
combined_labels = labels_left + labels_right
ax_left.legend(combined_lines, combined_labels, loc='center', bbox_to_anchor=(0.25, 0.60))

#Note: Save or display plot

#######
#Plot 2
#######
# Simulate data loading (use your own file)
def load_data(filepath):
    data = []
    with open(filepath, "r") as f:
        for line in f:
            a, b = line.strip().split()
            data.append((int(a), float(b)))
    return data

# File path
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
file_path = os.path.join(desktop, "data_anomalous.txt")

# Load data
data = load_data(file_path)
a_vals_raw = [x[0] for x in data]
a_max = max(a_vals_raw)
a_vals = [a / a_max for a in a_vals_raw]  # Normalize η to [0, 1]
b_vals = [x[1] for x in data]             # χ values

# Create plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(a_vals, b_vals, color='blue', linewidth=2)
ax.set_xlabel(r"$\eta$")
ax.set_ylabel(r"$\chi$")
ax.grid(True)

# Force y-axis to include -0.5
ymin, ymax = ax.get_ylim()
if ymin > -0.5:
    ymin = -0.5
ax.set_ylim(ymin, ymax)

#Note: Save or display plot