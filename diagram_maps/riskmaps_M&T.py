import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from collections import defaultdict
import numpy as np
import math
import os

# ============================================================
# Risk data: ID, Pre-L, Pre-C, Post-L, Post-C
# ============================================================

risks = [
    ("TR-MT-01", 2, 4, 1, 3),
    ("TR-MT-02", 4, 2, 2, 2),
    ("TR-MT-03", 3, 3, 2, 2),
    ("TR-MT-04", 4, 3, 2, 3),
    ("TR-MT-5", 3, 4, 2, 3),
    ("TR-MT-6", 3, 3, 2, 3),
    ("TR-MT-07", 3, 5, 2, 5),
    ("TR-MT-08", 4, 3, 3, 3),
    ("TR-MT-09", 5, 4, 2, 4),
    ("TR-MT-10", 3, 5, 2, 5),

    ("TR-M-01", 3, 5, 1, 5),
    ("TR-M-02", 3, 4, 2, 4),
    ("TR-M-03", 4, 5, 2, 5),
    ("TR-M-04", 4, 5, 1, 5),
    ("TR-M-05", 4, 4, 2, 4),
    ("TR-M-06", 4, 4, 1, 4),

    ("TR-T-01", 3, 4, 2, 3),
    ("TR-T-02", 3, 4, 2, 3),
    ("TR-T-03", 3, 4, 2, 4),
    ("TR-T-04", 3, 4, 2, 4),
    ("TR-T-05", 3, 3, 2, 3),
    ("TR-T-06", 4, 4, 1, 4),
]

# ============================================================
# Gradual vivid risk colour map
# ============================================================

risk_cmap = LinearSegmentedColormap.from_list(
    "risk_cmap",
    [
        "#00C853",  # vivid green
        "#FFF176",  # bright yellow
        "#FF9800",  # orange
        "#D50000",  # red
    ]
)

norm = Normalize(vmin=1, vmax=25)

def shortened_id(risk_id):
    return (
        risk_id
        .replace("TR-MT-", "MT")
        .replace("TR-M-", "M")
        .replace("TR-T-", "T")
    )

def plot_risk_map(ax, points, subtitle):
    # Continuous smooth background
    x = np.linspace(0.5, 5.5, 400)
    y = np.linspace(0.5, 5.5, 400)
    X, Y = np.meshgrid(x, y)
    Z = X * Y

    ax.imshow(
        Z,
        origin="lower",
        extent=[0.5, 5.5, 0.5, 5.5],
        cmap=risk_cmap,
        norm=norm,
        interpolation="bilinear",
        alpha=0.9,
        aspect="equal",
    )

    # Grid lines
    for v in np.arange(0.5, 6.0, 1.0):
        ax.axvline(v, color="black", linewidth=0.8)
        ax.axhline(v, color="black", linewidth=0.8)

    # Score numbers in each cell
    for likelihood in range(1, 6):
        for consequence in range(1, 6):
            score = likelihood * consequence
            ax.text(
                consequence,
                likelihood,
                str(score),
                ha="center",
                va="center",
                fontsize=8,
                alpha=0.35,
                fontweight="bold",
                color="black",
            )

    # Group risks in same cell
    grouped = defaultdict(list)
    for risk_id, likelihood, consequence in points:
        grouped[(consequence, likelihood)].append(shortened_id(risk_id))

    # Add risk labels
    for (consequence, likelihood), labels in grouped.items():
        n = len(labels)

        if n == 1:
            offsets = [(0, 0)]
        else:
            radius = 0.30
            offsets = [
                (
                    radius * math.cos(2 * math.pi * i / n),
                    radius * math.sin(2 * math.pi * i / n),
                )
                for i in range(n)
            ]

        for label, (dx, dy) in zip(labels, offsets):
            ax.text(
                consequence + dx,
                likelihood + dy,
                label,
                ha="center",
                va="center",
                fontsize=7,
                fontweight="bold",
                bbox=dict(
                    boxstyle="round,pad=0.16",
                    facecolor="white",
                    edgecolor="black",
                    linewidth=0.45,
                    alpha=0.95,
                ),
            )

    ax.set_xlim(0.5, 5.5)
    ax.set_ylim(0.5, 5.5)
    ax.set_xticks(range(1, 6))
    ax.set_yticks(range(1, 6))
    ax.set_xlabel("Consequence")
    ax.set_ylabel("Likelihood")
    ax.set_title(subtitle, fontsize=11, fontweight="bold")
    ax.set_aspect("equal")

# ============================================================
# Output figures
# ============================================================

output_dir = "risk_maps_output"
os.makedirs(output_dir, exist_ok=True)

pre_points = [(risk_id, pre_l, pre_c) for risk_id, pre_l, pre_c, post_l, post_c in risks]
post_points = [(risk_id, post_l, post_c) for risk_id, pre_l, pre_c, post_l, post_c in risks]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.6), constrained_layout=True)

plot_risk_map(axes[0], pre_points, "Pre-mitigation")
plot_risk_map(axes[1], post_points, "Post-mitigation")

# Colour bar
sm = plt.cm.ScalarMappable(cmap=risk_cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, shrink=0.82, pad=0.03)
cbar.set_label("Risk score", rotation=270, labelpad=15)

combined_path = os.path.join(output_dir, "risk_maps_midterm.png")
plt.savefig(combined_path, dpi=300, bbox_inches="tight")
plt.show()

print(f"Saved: {combined_path}")