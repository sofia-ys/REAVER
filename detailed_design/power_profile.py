import os
import pandas as pd
import matplotlib.pyplot as plt

os.makedirs("files", exist_ok=True)

modes = [
    "Safe Mode",
    "Commissioning",
    "LR Rendezvous",
    "Proximity Ops",
    "Capture",
    "Standby",
]

# -----------------------------------------------------------------------------
# 1. MOTHERSHIP POWER PROFILE
# -----------------------------------------------------------------------------
mothership_data = {
    "TCS":        [0.00,   0.00,   0.00,   110.00, 0.00,   0.00],
    "Payload":    [0.00,   0.00,   0.00,   381.75, 381.75, 3.36],
    "C&DH":       [132.00, 132.00, 132.00, 132.00, 132.00, 132.00],
    "Propulsion": [0.00,   0.00,   212.52, 480.00, 480.00, 0.00],
    "AOCS":       [49.98,  49.98,  49.98,  49.98,  49.98,  49.98],
    "TT&C":       [8.40,   8.40,   8.40,   50.40,  50.40,  8.40],
}

df_ms = pd.DataFrame(mothership_data, index=modes)
df_ms = df_ms.loc[:, (df_ms != 0).any(axis=0)]
ms_cbe = df_ms.sum(axis=1)
df_ms["20% Contingency"] = ms_cbe * 0.20

fig1, ax1 = plt.subplots(figsize=(8, 6))
df_ms.plot(kind="bar", stacked=True, ax=ax1, colormap="Set2", edgecolor="black", linewidth=0.5)

for patch in ax1.containers[-1].patches:
    patch.set_hatch("///")
    patch.set_edgecolor("crimson")
    patch.set_facecolor("#ffebee")

ax1.set_ylabel("Power (W)", fontsize=12)
ax1.set_xlabel("Operational Mode", fontsize=12)
ax1.grid(axis="y", linestyle="--", alpha=0.5, zorder=0)
ax1.set_axisbelow(True)
handles1, labels1 = ax1.get_legend_handles_labels()
ax1.legend(handles1, labels1, loc="upper left", title="Subsystems")
plt.xticks(rotation=15, ha="right")
plt.tight_layout()
plt.savefig("files/mothership_power_profile.png", dpi=300)
plt.close(fig1)

# -----------------------------------------------------------------------------
# 2. TUG POWER PROFILE
# -----------------------------------------------------------------------------
tug_data = {
    "TCS":        [0.00,  0.00,  22.00,   0.00,  0.00,  0.00],
    "Payload":    [0.00,  0.00,  0.00,    0.00,  62.10, 0.00],
    "C&DH":       [31.90, 31.90, 31.90,   31.90, 31.90, 31.90],
    "Propulsion": [0.00,  0.00,  1172.01, 16.38, 16.38, 0.00],
    "AOCS":       [54.71, 54.71, 54.71,   54.71, 54.71, 54.71],
    "TT&C":       [4.20,  4.20,  4.20,    8.40,  8.40,  4.20],
}

df_tug = pd.DataFrame(tug_data, index=modes)
df_tug = df_tug.loc[:, (df_tug != 0).any(axis=0)]
tug_cbe = df_tug.sum(axis=1)
df_tug["20% Contingency"] = tug_cbe * 0.20

fig2, ax2 = plt.subplots(figsize=(8, 6))
df_tug.plot(kind="bar", stacked=True, ax=ax2, colormap="Set2", edgecolor="black", linewidth=0.5)

for patch in ax2.containers[-1].patches:
    patch.set_hatch("///")
    patch.set_edgecolor("crimson")
    patch.set_facecolor("#ffebee")

ax2.set_ylabel("Power (W)", fontsize=12)
ax2.set_xlabel("Operational Mode", fontsize=12)
ax2.grid(axis="y", linestyle="--", alpha=0.5, zorder=0)
ax2.set_axisbelow(True)
handles2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(handles2, labels2, loc="upper left", title="Subsystems")
plt.xticks(rotation=15, ha="right")
plt.tight_layout()
plt.savefig("files/tug_power_profile.png", dpi=300)
plt.close(fig2)

print("Power profile plots saved to files/")