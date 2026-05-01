import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch

# Create figure
fig, ax = plt.subplots(figsize=(13, 8))

# Axis limits
ax.set_xlim(0, 5)
ax.set_ylim(0, 5)

# Background quadrants
ax.add_patch(Rectangle((0, 3), 2.5, 2, facecolor="#d8efd3", edgecolor="#4caf50", linewidth=1.5))
ax.add_patch(Rectangle((2.5, 3), 2.5, 2, facecolor="#fff0c2", edgecolor="#ffa500", linewidth=1.5))
ax.add_patch(Rectangle((0, 0), 2.5, 3, facecolor="#f8caca", edgecolor="#ef5350", linewidth=1.5))
ax.add_patch(Rectangle((2.5, 0), 2.5, 3, facecolor="#cfe0ff", edgecolor="#3f7ee8", linewidth=1.5))

# Function to add stakeholder boxes
def stakeholder_box(x, y, title, subtitle, color):
    box = FancyBboxPatch(
        (x, y), 1.45, 0.55,
        boxstyle="round,pad=0.06,rounding_size=0.06",
        facecolor="white",
        edgecolor=color,
        linewidth=1.1
    )
    ax.add_patch(box)

    ax.text(x + 0.18, y + 0.28, "●", fontsize=28, color=color,
            ha="center", va="center")
    ax.text(x + 0.40, y + 0.36, title, fontsize=11, fontweight="bold",
            ha="left", va="center")
    ax.text(x + 0.40, y + 0.17, subtitle, fontsize=10,
            ha="left", va="center")

# Quadrant headings
ax.text(1.25, 4.75, "KEEP SATISFIED", ha="center", va="center",
        fontsize=16, fontweight="bold")
ax.text(1.25, 4.52, "High Interest – Low Influence", ha="center",
        va="center", fontsize=12, style="italic")

ax.text(3.75, 4.75, "MANAGE CLOSELY", ha="center", va="center",
        fontsize=16, fontweight="bold")
ax.text(3.75, 4.52, "High Interest – High Influence", ha="center",
        va="center", fontsize=12, style="italic")

ax.text(1.25, 2.75, "MONITOR", ha="center", va="center",
        fontsize=16, fontweight="bold")
ax.text(1.25, 2.52, "Low Interest – Low Influence", ha="center",
        va="center", fontsize=12, style="italic")

ax.text(3.75, 2.75, "KEEP INFORMED", ha="center", va="center",
        fontsize=16, fontweight="bold")
ax.text(3.75, 2.52, "Low Interest – High Influence", ha="center",
        va="center", fontsize=12, style="italic")

# Stakeholder boxes
stakeholder_box(0.60, 4.05, "Regulatory Bodies",
                "ITU, National Space\nAuthorities", "#4caf50")

stakeholder_box(0.60, 3.25, "National Government\n& Policy Makers",
                "Non-primary space agencies", "#4caf50")

stakeholder_box(3.10, 4.05, "Space Agencies",
                "ESA, NASA, JAXA", "#ffb300")

stakeholder_box(3.10, 3.25, "Satellite Operators",
                "Commercial & government\nsatellite owners", "#ff9800")

stakeholder_box(0.60, 0.95, "General Public",
                "Society at large,\ntaxpayers", "#e53935")

stakeholder_box(3.10, 1.65, "Scientific Community",
                "Universities, research\ninstitutes", "#3f7ee8")

stakeholder_box(3.10, 0.85, "Insurance Companies",
                "Space insurance\nproviders", "#3f7ee8")

# Axis labels
ax.set_xlabel("Influence on Project Outcome", fontsize=16, fontweight="bold", labelpad=18)
ax.set_ylabel("Interest in Project Outcome", fontsize=16, fontweight="bold", labelpad=18)

# Ticks
ax.set_xticks([1, 2, 3, 4, 5])
ax.set_yticks([1, 2, 3, 4, 5])
ax.tick_params(axis="both", labelsize=13, length=0)

# Add Low/High labels
ax.text(0.02, -0.32, "Low\n(1)", ha="left", va="top", fontsize=12)
ax.text(4.98, -0.32, "High\n(5)", ha="right", va="top", fontsize=12)
ax.text(-0.30, 0.02, "Low\n(1)", ha="left", va="bottom", fontsize=12)
ax.text(-0.30, 4.98, "High\n(5)", ha="left", va="top", fontsize=12)

# Clean frame
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)

# Keep square proportions
ax.set_aspect("auto")

# Save for Overleaf
plt.savefig("stakeholder_interest_influence_map.pdf", bbox_inches="tight")
plt.savefig("stakeholder_interest_influence_map.png", dpi=300, bbox_inches="tight")

plt.show()