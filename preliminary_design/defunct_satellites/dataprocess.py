import pandas as pd
import matplotlib.pyplot as plt
import numpy as np75
from pathlib import Path

# ============================================================
# 1. Load file
# ============================================================

file_path = "reaver_defunct_candidates_gcat.csv"
df = pd.read_csv(file_path)

# Ask user for concentration percentage
percentage = float(input("Enter percentage of most concentrated spacecraft to keep, e.g. 100, 75, 50: "))

if percentage <= 0 or percentage > 100:
    raise ValueError("Percentage must be between 0 and 100.")

# Output folder
output_dir = Path(f"orbital_distribution_plots_top_{int(percentage)}percent")
output_dir.mkdir(exist_ok=True)

# ============================================================
# 2. Orbital columns
# ============================================================

orbital_cols = [
    "SEMIMAJOR_AXIS",
    "ECCENTRICITY",
    "INCLINATION",
    "RA_OF_ASC_NODE",
    "ARG_OF_PERICENTER",
    "MEAN_ANOMALY",
    "APOAPSIS",
    "PERIAPSIS"
]

# Convert columns to numeric
for col in orbital_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# Remove rows with missing orbital data
df = df.dropna(subset=orbital_cols).copy()

print(f"\nOriginal number of valid objects: {len(df)}")

# ============================================================
# 3. Calculate local density
# ============================================================
# We use normalised orbital elements and calculate how close
# each spacecraft is to its nearest neighbours.
# Smaller average distance = more concentrated region.

density_features = [
    "SEMIMAJOR_AXIS",
    "ECCENTRICITY",
    "INCLINATION",
    "RA_OF_ASC_NODE",
    "ARG_OF_PERICENTER",
    "MEAN_ANOMALY"
]

X = df[density_features].copy()

# Handle angular variables properly using sine/cosine
# This avoids problems where 359 deg and 1 deg are actually close.
angle_cols = [
    "RA_OF_ASC_NODE",
    "ARG_OF_PERICENTER",
    "MEAN_ANOMALY"
]

X_expanded = pd.DataFrame(index=df.index)

for col in density_features:
    if col in angle_cols:
        angle_rad = np.deg2rad(X[col])
        X_expanded[col + "_sin"] = np.sin(angle_rad)
        X_expanded[col + "_cos"] = np.cos(angle_rad)
    else:
        X_expanded[col] = X[col]

# Normalise each feature
X_norm = (X_expanded - X_expanded.mean()) / X_expanded.std()
X_norm = X_norm.fillna(0)

# Convert to numpy array
X_array = X_norm.to_numpy()

# Compute distance matrix
# For this dataset size, this is simple and works fine.
diff = X_array[:, None, :] - X_array[None, :, :]
dist_matrix = np.sqrt(np.sum(diff**2, axis=2))

# Ignore distance to itself
np.fill_diagonal(dist_matrix, np.inf)

# Number of neighbours used for density estimate
k = max(5, int(0.05 * len(df)))  # 5% of dataset, minimum 5 neighbours

# Average distance to k nearest neighbours
nearest_distances = np.sort(dist_matrix, axis=1)[:, :k]
density_score = nearest_distances.mean(axis=1)

df["DENSITY_SCORE"] = density_score

# Lower density score = more concentrated
df = df.sort_values("DENSITY_SCORE", ascending=True)

# Keep selected percentage
n_keep = int(np.ceil(len(df) * percentage / 100))
df_selected = df.head(n_keep).copy()

print(f"Selected number of objects: {len(df_selected)}")
print(f"Keeping densest {percentage:.1f}% of the distribution")

# Save selected data
output_csv = f"reaver_top_{int(percentage)}percent_concentrated_objects.csv"
df_selected.to_csv(output_csv, index=False)

print(f"Selected objects saved as: {output_csv}")

# ============================================================
# 4. Basic statistics
# ============================================================

print("\nBasic orbital statistics for selected objects:")
print(df_selected[orbital_cols].describe())

# ============================================================
# 5. Plotting functions
# ============================================================

def plot_histogram(data, column, xlabel, bins=30):
    plt.figure(figsize=(8, 5))
    plt.hist(data[column], bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel("Number of objects")
    plt.title(f"Distribution of {xlabel} - Top {percentage:.0f}% concentration")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f"{column}_histogram.png", dpi=300)
    plt.show()


def plot_2d_histogram(data, x_col, y_col, xlabel, ylabel, bins=30):
    plt.figure(figsize=(8, 6))
    plt.hist2d(data[x_col], data[y_col], bins=bins)
    plt.colorbar(label="Number of objects")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"Concentration: {xlabel} vs {ylabel} - Top {percentage:.0f}%")
    plt.tight_layout()
    plt.savefig(output_dir / f"{x_col}_vs_{y_col}_2d.png", dpi=300)
    plt.show()


def plot_scatter(data, x_col, y_col, xlabel, ylabel):
    plt.figure(figsize=(8, 6))
    plt.scatter(data[x_col], data[y_col], s=15)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{xlabel} vs {ylabel} - Top {percentage:.0f}% concentration")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f"{x_col}_vs_{y_col}_scatter.png", dpi=300)
    plt.show()

# ============================================================
# 6. Histograms of Keplerian elements
# ============================================================

plot_histogram(df_selected, "SEMIMAJOR_AXIS", "Semi-major axis [km]")
plot_histogram(df_selected, "ECCENTRICITY", "Eccentricity [-]")
plot_histogram(df_selected, "INCLINATION", "Inclination [deg]")
plot_histogram(df_selected, "RA_OF_ASC_NODE", "RAAN [deg]")
plot_histogram(df_selected, "ARG_OF_PERICENTER", "Argument of pericentre [deg]")
plot_histogram(df_selected, "MEAN_ANOMALY", "Mean anomaly [deg]")
plot_histogram(df_selected, "APOAPSIS", "Apogee altitude [km]")
plot_histogram(df_selected, "PERIAPSIS", "Perigee altitude [km]")

# ============================================================
# 7. 2D concentration plots
# ============================================================

plot_2d_histogram(
    df_selected,
    "RA_OF_ASC_NODE",
    "INCLINATION",
    "RAAN [deg]",
    "Inclination [deg]"
)

plot_2d_histogram(
    df_selected,
    "SEMIMAJOR_AXIS",
    "INCLINATION",
    "Semi-major axis [km]",
    "Inclination [deg]"
)

plot_2d_histogram(
    df_selected,
    "RA_OF_ASC_NODE",
    "SEMIMAJOR_AXIS",
    "RAAN [deg]",
    "Semi-major axis [km]"
)

plot_2d_histogram(
    df_selected,
    "APOAPSIS",
    "PERIAPSIS",
    "Apogee altitude [km]",
    "Perigee altitude [km]"
)

# Scatter plots are useful to see individual spacecraft
plot_scatter(
    df_selected,
    "RA_OF_ASC_NODE",
    "INCLINATION",
    "RAAN [deg]",
    "Inclination [deg]"
)

plot_scatter(
    df_selected,
    "SEMIMAJOR_AXIS",
    "INCLINATION",
    "Semi-major axis [km]",
    "Inclination [deg]"
)

# ============================================================
# 8. Print strongest bins
# ============================================================

def print_top_bins(data, column, bins=20):
    values = data[column].dropna()
    counts, edges = np.histogram(values, bins=bins)

    top_indices = np.argsort(counts)[::-1][:5]

    print(f"\nTop concentration bins for {column}:")
    for i in top_indices:
        print(
            f"{edges[i]:.3f} to {edges[i+1]:.3f}: "
            f"{counts[i]} objects"
        )

for col in orbital_cols:
    print_top_bins(df_selected, col, bins=20)

print("\nPlots saved in:", output_dir)
print("Finished.")