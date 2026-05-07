import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import defaultdict
from adjustText import adjust_text

# Risk data: (ID, Pre-L, Pre-C, Post-L, Post-C)
risks = [
    ("TR-01", 4, 5, 2, 5),
    ("TR-02", 3, 5, 2, 4),
    ("TR-03", 3, 5, 2, 5),
    ("TR-04", 3, 4, 2, 4),
    ("TR-05", 3, 4, 2, 4),
    ("TR-06", 2, 4, 1, 4),
    ("TR-07", 2, 5, 1, 5),
    ("TR-08", 2, 5, 1, 5),
    ("TR-09", 2, 4, 1, 4),
    ("TR-10", 3, 3, 2, 3),
    ("TR-11", 3, 4, 2, 4),
    ("TR-12", 3, 4, 2, 4),
    ("TR-13", 4, 4, 2, 4),
    ("TR-14", 2, 5, 1, 5),
    ("TR-15", 3, 4, 2, 4),
    ("TR-16", 3, 4, 2, 4),
    ("TR-17", 2, 4, 1, 4),
    ("TR-18", 3, 4, 2, 4),
    ("TR-19", 3, 4, 2, 4),
    ("TR-20", 3, 3, 2, 3),
    ("TR-21", 2, 4, 1, 4),
    ("TR-22", 2, 3, 1, 3),
    ("TR-23", 2, 4, 1, 4),
    ("TR-24", 2, 4, 1, 4),
    ("TR-25", 1, 5, 1, 5),
    ("TR-26", 3, 4, 2, 3),
    ("TR-27", 2, 4, 1, 4),
    ("TR-28", 2, 4, 1, 4),
    ("TR-29", 2, 5, 2, 5),
    ("TR-30", 2, 4, 1, 4),
]


def group_by_coord(risk_list, l_idx, c_idx):
    """Return dict {(L, C): [IDs]} sorted by ID."""
    clusters = defaultdict(list)
    for r in risk_list:
        key = (r[l_idx], r[c_idx])
        clusters[key].append(r[0])
    for key in clusters:
        clusters[key].sort()
    return clusters


def format_label(id_list, ids_per_row=3):
    rows = [
        " ".join(id_list[i : i + ids_per_row])
        for i in range(0, len(id_list), ids_per_row)
    ]
    return "\n".join(rows)


def make_risk_background(ax):
    res = 400
    x = np.linspace(0, 5, res)
    y = np.linspace(0, 5, res)
    X, Y = np.meshgrid(x, y)
    Z = X * Y

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "risk",
        [
            (0.00, "#56b400"),
            (0.30, "#c8d400"),
            (0.50, "#f5c400"),
            (0.70, "#f07800"),
            (1.00, "#c80000"),
        ],
    )
    ax.imshow(
        Z,
        extent=[0, 5, 0, 5],
        origin="lower",
        cmap=cmap,
        vmin=0,
        vmax=25,
        aspect="auto",
        alpha=0.88,
        zorder=0,
    )


def draw_grid(ax):
    for v in range(1, 6):
        ax.axhline(v, color="white", linewidth=0.7, zorder=1)
        ax.axvline(v, color="white", linewidth=0.7, zorder=1)


def plot_riskmap(ax, clusters, title):
    make_risk_background(ax)
    draw_grid(ax)

    lx = [k[0] for k in clusters]
    cy = [k[1] for k in clusters]

    ax.scatter(
        lx, cy,
        s=55,
        color="#1a6fa8",
        edgecolors="white",
        linewidths=0.8,
        zorder=5,
    )

    texts = []
    for (lv, cv), id_list in clusters.items():
        ids_per_row = 2 if len(id_list) > 6 else 3
        label = format_label(id_list, ids_per_row=ids_per_row)
        t = ax.text(
            lv, cv,
            label,
            fontsize=6.0,
            fontweight="bold",
            color="#111111",
            ha="center",
            va="bottom",
            zorder=7,
            linespacing=1.25,
            bbox=dict(
                boxstyle="round,pad=0.18",
                fc="white",
                ec="none",
                alpha=0.72,
            ),
        )
        texts.append(t)

    adjust_text(
        texts,
        x=np.array(lx, dtype=float),
        y=np.array(cy, dtype=float),
        ax=ax,
        expand_points=(2.0, 2.0),
        expand_text=(1.5, 1.5),
        force_points=(1.2, 1.2),
        force_text=(0.8, 0.8),
        avoid_points=True,
        avoid_self=True,
        only_move={"points": "xy", "text": "xy"},
        arrowprops=None,
        lim=1500,
    )

    ax.set_xlim(0, 5)
    ax.set_ylim(0, 5.8)
    ax.set_xticks(range(6))
    ax.set_yticks(range(6))
    ax.tick_params(labelsize=10)
    ax.set_xlabel("Likelihood [L]", fontsize=12)
    ax.set_ylabel("Consequence [C]", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)


pre_clusters  = group_by_coord(risks, l_idx=1, c_idx=2)
post_clusters = group_by_coord(risks, l_idx=3, c_idx=4)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 9))
fig.patch.set_facecolor("white")

plot_riskmap(ax1, pre_clusters,  "Risk Map Before Mitigation")
plot_riskmap(ax2, post_clusters, "Risk Map After Mitigation")

plt.tight_layout(pad=3.0)
plt.savefig("riskmap_output.png", dpi=180, bbox_inches="tight",
            facecolor="white")
plt.show()
print("Saved: riskmap_output.png")
