"""
REAVER RPO Simulation — Input Parameters
========================================
Input parameters for the close-range RPO ΔV simulation (Parts 1 & 2 of the
Rev 2 build checklist).  These are *inputs* to the simulation.

The simulation OUTPUT — the five ΔV constants
    DV_RPO_DEBRIS, DV_RPO_DETUMBLE, DV_RPO_RH_MEET, DV_RPO_RH_DOCK, DV_RPO_RH
is printed by rpo_run.py for pasting into config.py (Part 3).

Anything tagged  ``# TBD MIRON``  or  ``# PLACEHOLDER``  is a flagged
assumption to be refined when the design team provides the real value.

References:
  Conings & Mooij, "Integrated GNC System Design for Active Debris Removal",
  AIAA-2025-0085.  Table 1 (guidance), Table 2 (control), Table 4 (VBN).
"""
import numpy as np
from config import MU, RH_SMA, MS_DRY, TUG_DRY, MS_VEX, DAY

# ── Target orbit mean motion (drives all CW dynamics) ─────────────────────────
# At GEO this is ~15x smaller than the LEO value (0.063 deg/s) used in the paper.
N_MEAN = float(np.sqrt(MU / RH_SMA**3))          # rad/s  ≈ 7.1e-5

# ── Navigation sensor handover (Rev 2 RPO Logistics) ──────────────────────────
R_CUTOFF = 150.0        # m   RF/UKF -> LiDAR/VBN handover (100-200 m band).
                        #     This is the START boundary of the RPO simulation.

# ── Close-range geometry (Conings Table 1, mapped to REAVER) ─────────────────
H_BAR_STANDOFF   = 20.0     # m   H-bar hold-point standoff
R_KOS1           = 5.0      # m   keep-out sphere 1 = interaction distance (arm reach)  # TBD MIRON
R_KOS2           = 1.35     # m   keep-out sphere 2 (final approach)
APPROACH_CONE_DEG = 30.0    # deg approach-cone half-angle (Conings Eq. 17)
MEV_LOCK_DIST    = 1.35     # m   MEV engagement distance from LAE nozzle           # TBD MIRON

# ── Translational limits (Conings Table 1) ────────────────────────────────────
V_MAX  = 1.0        # m/s   max approach velocity (low-impact docking requirement)
T_MAX  = 220.0      # N     max commanded thrust (MS thruster level)              # PLACEHOLDER

# ── Attitude-control actuators (Conings Table 2) ──────────────────────────────
RCS_THR = 55.0      # N     RCS thruster nominal force
RCS_ARM = 1.0       # m     effective moment arm (≈ MS half-width)                # TBD MIRON
DW_SYNC_DPS = 0.6   # deg/s P->INDI switch threshold (Conings Eq. 19)

# ── GEO perturbations (replace LEO drag/magnetic from the paper) ──────────────
P_SOL    = 4.56e-6  # N/m^2 solar-radiation pressure at 1 AU
C_R      = 1.5      #       reflectivity coefficient (MLI-covered bus)
AOM_MS   = 0.015    # m^2/kg area-to-mass ratio of the mothership                 # PLACEHOLDER
A_TRIAX  = 3.0e-8   # m/s^2 representative Earth-triaxiality (J22) tangential accel

# ── Tumble state (REQ-MIS-06) ─────────────────────────────────────────────────
OMEGA_REQ_DPS    = 6.0                 # deg/s  1 rpm requirement ceiling
TUMBLE_SWEEP_DPS = [0.0, 3.0, 6.0]     # deg/s  parametric sweep (Conings Table 5)

# ── Masses ────────────────────────────────────────────────────────────────────
M_DEBRIS_WC   = 2700.0                       # kg  EUTE 12 WEST A (worst case)
M_TUG         = float(TUG_DRY)               # kg  tug dry mass
M_TUGDEB_WC   = M_TUG + M_DEBRIS_WC          # kg  tug+debris combined (Phase A target)
M_COMBINED_RH = float(MS_DRY) + M_TUGDEB_WC  # kg  MS+tug+debris (Phase B chaser)

# ── Body geometry for inertia tensors (box models) ────────────────────────────
MS_DIMS            = (2.0, 2.0, 2.5)    # m   MS bounding box (x,y,z)              # PLACEHOLDER
TUG_DIMS           = (1.0, 1.0, 1.0)    # m   tug bounding box                    # PLACEHOLDER
DEBRIS_BUS_DIMS    = (2.5, 2.5, 3.0)    # m   debris central bus                  # PLACEHOLDER
DEBRIS_PANEL_SPAN  = 20.0               # m   solar-panel tip-to-tip span         # PLACEHOLDER
DEBRIS_PANEL_FRAC  = 0.15               #     fraction of debris mass in panels
ARM_LEN_EXTENDED   = 5.0                # m   robotic-arm length when extended     # TBD MIRON

# ── Dwell / timing (Rev 2: feeds the T_OPS update) ───────────────────────────
T_DWELL_ANALYSIS = 1.0 * 3600.0    # s   debris-analysis dwell (sensor scan)      # PLACEHOLDER
T_HOLD_RH        = 0.5 * 3600.0    # s   MS hold at RH between tug arrivals        # PLACEHOLDER

# ── Risk margins ──────────────────────────────────────────────────────────────
ABORT_MARGIN = 0.50    #  +50% abort/retry margin on debris RPO (risk TR-MT-08)

__all__ = [
    'N_MEAN', 'R_CUTOFF', 'H_BAR_STANDOFF', 'R_KOS1', 'R_KOS2',
    'APPROACH_CONE_DEG', 'MEV_LOCK_DIST', 'V_MAX', 'T_MAX',
    'RCS_THR', 'RCS_ARM', 'DW_SYNC_DPS',
    'P_SOL', 'C_R', 'AOM_MS', 'A_TRIAX',
    'OMEGA_REQ_DPS', 'TUMBLE_SWEEP_DPS',
    'M_DEBRIS_WC', 'M_TUG', 'M_TUGDEB_WC', 'M_COMBINED_RH',
    'MS_DIMS', 'TUG_DIMS', 'DEBRIS_BUS_DIMS', 'DEBRIS_PANEL_SPAN',
    'DEBRIS_PANEL_FRAC', 'ARM_LEN_EXTENDED',
    'T_DWELL_ANALYSIS', 'T_HOLD_RH', 'ABORT_MARGIN',
]
