"""
REAVER RPO Simulation — Input Parameters
========================================
Input parameters for the close-range RPO ΔV simulation (Parts 1 & 2 of the
Rev 2 build checklist).  These are *inputs* to the simulation.

The simulation OUTPUT — the five ΔV constants
    DV_RPO_DEBRIS, DV_RPO_DETUMBLE, DV_RPO_TUG_MEET, DV_RPO_RH_DOCK, DV_RPO_RH
is printed by rpo_run.py for pasting into config.py (Part 3).

Anything tagged  ``# TBD MIRON``  or  ``# PLACEHOLDER``  is a flagged
assumption to be refined when the design team provides the real value.

References:
  Conings & Mooij, "Integrated GNC System Design for Active Debris Removal",
  AIAA-2025-0085.  Table 1 (guidance), Table 2 (control), Table 4 (VBN).
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
import numpy as np
from config import MU, RH_SMA, MS_DRY, TUG_DRY, MS_VEX, DAY

# ── Target orbit mean motion (drives all CW dynamics) ─────────────────────────
# At GEO this is ~15x smaller than the LEO value (0.063 deg/s) used in the paper.
N_MEAN = float(np.sqrt(MU / RH_SMA**3))          # rad/s  ≈ 7.1e-5

# ── Navigation sensor handover (Rev 2 RPO Logistics) ──────────────────────────
R_CUTOFF = 150.0        # m   RF/UKF -> LiDAR/VBN handover (100-200 m band).
                        #     This is the START boundary of the RPO simulation.

# ── Close-range geometry (Conings Table 1, mapped to REAVER) ─────────────────
# Phase IDs from the RPO Logistics sequence diagram:
#   Approach + Interaction = MD-MX-PRX, Capture & control = MD-TUG-CAP,
#   Debris retrieval / retreat = MD-MS-LRR.
# Hold + inspection are done on V-bar (the only naturally fixed CW hold point:
# a body at rest at an along-track offset stays put, whereas an H-bar offset
# oscillates through the orbital plane).  The capture corridor is along the
# IDENTIFIED tumble axis (modelled as +H-bar); the MS only transitions there
# after the V-bar inspection has measured the spin rate + tumble axis.
INSPECT_STANDOFF = 20.0     # m   V-bar inspection hold standoff (KOS1-clear, fixed hold)
H_BAR_STANDOFF   = 20.0     # m   tumble-axis (H-bar) capture standoff
ARM_LEN_EXTENDED = 4.627      # m   robotic-arm length when extended

# ── KOS1: operational definition (NOT a mission constant) ────────────────────
# KOS1 is set PER TARGET at the cut-off distance (150-200 m) based on the
# actual debris solar-panel span uplinked by GO to the MS at the RF->LiDAR
# navigation handover.  GO has characterised the target before RPO begins and
# sends the confirmed panel span so the MS can configure its hold distance
# correctly for that specific debris.
#
# KOS1 = panel_span / 2  (half the tip-to-tip solar array wingspan)
# This is the bounding-sphere radius of the tumbling debris and ensures the
# MS stays outside the panel sweep zone before spin synchronisation.
#
# Known values (from public satellite databases):
#   COMSATBW-1     (BSS-702)           : span = 17.2 m  -> KOS1 = 17.2/2 * 1.20 = 10.3 m
#   EUTE 12 WEST A (Astra 1B, HS-601)  : span = 26.2 m  -> KOS1 = TBD (not in sequence)
#   All others                          : span = TBD (GO uplink at cut-off)
#
# Simulation reference target: COMSATBW-1 (confirmed span, part of optimal sequence).
# KOS1 = (panel_span / 2) * KOS1_MARGIN — margin applied in rpo_debris_sim.py.
DEBRIS_PANEL_SPAN_WC = 17.2   # m   COMSATBW-1 (BSS-702) confirmed span
KOS1_MARGIN          = 1.20   #     20% safety margin on panel half-span

R_SPIN_STANDOFF  = 2.0      # m   spin-up / capture standoff (Mechanics of Space Debris
                            #     Removal review, pg.26: spin-up & align at 2 m)
APPROACH_CONE_DEG = 5.0     # deg approach-cone half-angle — terminal docking corridor
                            #     abort if lateral offset > KOS2*tan(5°) = 0.49 m at KOS2
                            #     (MEV-1/Intelsat-901 analog, LAE nozzle ~0.35 m diam)
R_RETREAT        = 20.0     # m   safe standoff the MS retreats to after tug release
                            #     (sequence-diagram step 4c, before long-range departure)

# ── Translational limits (Conings Table 1) ────────────────────────────────────
V_MAX  = 0.05        # m/s   max approach velocity (low-impact docking requirement)

# ── Attitude-control actuators ────────────────────────────────────────────────
# MS has an octagonal cross-section with 24 RCS thrusters total:
#   12 on the top  edges of the octagon (4 active edges, 3 thrusters per edge)
#   12 on the bottom edges              (4 active edges, 3 thrusters per edge)
# For detumble the torque couple uses ONE opposed edge pair (conservative):
#   tau = 2 * N_RCS_PER_EDGE * F_RCS * RCS_ARM
#       = 2 * 3 * 4.0 * 1.05 = 25.2 N.m
# Propellant formula m_p = H / (RCS_ARM * vex) is independent of thruster count.
F_RCS          = 4.0    # N     single RCS thruster force (confirmed)
N_RCS_TOTAL    = 24     #       total RCS thrusters on MS (confirmed)
N_RCS_PER_EDGE = 3      #       thrusters per active edge (confirmed)
# RCS_THR = force contribution from one edge of the attitude couple (used in torque formula)
RCS_THR        = float(N_RCS_PER_EDGE * F_RCS)   # = 12.0 N
RCS_ARM        = 1.05   # m     minimum effective moment arm (confirmed)
DW_SYNC_DPS    = 0.6    # deg/s P->INDI switch threshold (Conings Eq. 19)

# ── GEO perturbations (replace LEO drag/magnetic from the paper) ──────────────
P_SOL    = 4.56e-6  # N/m^2 solar-radiation pressure at 1 AU
C_R      = 1.5      #       reflectivity coefficient (MLI-covered bus)
AOM_MS   = 0.020    # m^2/kg area-to-mass ratio of the mothership                 # PLACEHOLDER
A_TRIAX  = 3.0e-8   # m/s^2 representative Earth-triaxiality (J22) tangential accel

# ── Tumble state (REQ-MIS-06) ─────────────────────────────────────────────────
OMEGA_REQ_DPS    = 6.0                 # deg/s  1 rpm requirement ceiling
TUMBLE_SWEEP_DPS = [0.0, 3.0, 6.0]     # deg/s  parametric sweep (Conings Table 5)

# ── Masses ────────────────────────────────────────────────────────────────────
M_DEBRIS_WC   = 2440.0                       # kg  COMSATBW-1 (reference target, confirmed data)
M_TUG         = float(TUG_DRY)               # kg  tug dry mass
M_TUGDEB_WC   = M_TUG + M_DEBRIS_WC          # kg  tug+debris combined (Phase A target)
M_COMBINED_RH = float(MS_DRY) + M_TUGDEB_WC  # kg  MS+tug+debris (Phase B chaser)

# ── Body geometry for inertia tensors (box models) ────────────────────────────
MS_DIMS            = (3.0, 3.0, 3.5)    # m   MS bounding box (x,y,z)  (confirmed)
TUG_DIMS           = (0.75, 0.75, 1.2)   # m   tug bounding box (x,y,z) (confirmed)
DEBRIS_BUS_DIMS    = (2.5, 2.5, 3.0)    # m   debris central bus                  # PLACEHOLDER
DEBRIS_PANEL_FRAC  = 0.15               #     fraction of debris mass in panels

# KOS2: inner sphere — MS CoM distance from debris CoM when tug MEV locks LAE.
# = arm length (TBD MIRON) + tug half-length.  Fixed geometry, not target-specific.
# KOS2 < KOS1 by construction (arm < panel half-span for any reasonable arm length).
R_KOS2        = ARM_LEN_EXTENDED + TUG_DIMS[2] / 2.0  # m (derived)  # TBD MIRON
MEV_LOCK_DIST = R_KOS2                                  # m same position by definition

# KOS2 for the TUG MEET (Phase A): the MS reaches its BARE arm out to grab the
# cooperative tug (no tug stowed on the MS arm yet), so MS-CoM to tug-CoM is the
# arm length ONLY — the tug half-length is NOT added here (contrast R_KOS2 above,
# where the MS carries a tug whose MEV reaches the debris).
R_KOS2_MEET   = ARM_LEN_EXTENDED                       # m (arm only)  # TBD MIRON

# ── RH docking interface (Phase B: DOCK) ─────────────────────────────────────
# The MS docks to the RH with a dedicated docking port mounted on TOP of the
# bus — NOT the robotic arm.  Throughout the RH approach the arm stays EXTENDED,
# holding the captured tug+debris out to the side, so the chaser is a large
# off-balance body docking nose-first via its top port.  Two keep-out spheres,
# centred on the RH docking port:
DOCK_PORT_HEIGHT = 0.40   # m   MS top-mounted docking-port height (confirmed)

# KOS2_DOCK (inner / hard-dock standoff): MS CoM -> RH port at port-to-port
# contact = half the MS height along the dock axis + the port height.
R_KOS2_DOCK = MS_DIMS[2] / 2.0 + DOCK_PORT_HEIGHT      # = 1.75 + 0.40 = 2.15 m

# KOS1_DOCK (outer / alignment + go-no-go hold): the MS must hold far enough out
# that, while it slews to line its top port up with the RH port, the EXTENDED
# arm + tug + debris it is carrying clears the RH structure.  Suggested sizing =
# the chaser's extended envelope (arm reach + tug + a debris-bus clearance) with
# the standard KOS margin.  Assumes the debris solar panels are oriented along
# the approach axis; if they sweep broadside this must grow to the panel span.
# RH's own structural bounding sphere is assumed to fit inside this envelope.
R_KOS1_DOCK = (ARM_LEN_EXTENDED + TUG_DIMS[2] + DEBRIS_BUS_DIMS[2]) * KOS1_MARGIN  # ~11.0 m  # TBD MIRON

# ── Dwell / timing (Rev 2: feeds the T_OPS update) ───────────────────────────
# T_DWELL_ANALYSIS: time at interaction distance for debris characterisation
# AND the full ground relay cycle:
#   - Onboard C&DH sensor processing  (~30 min)
#   - Telemetry downlink + ground review (~6-7 h)
#   - Go/no-go uplink back to MS        (~15 min)
# Light-time delay at GEO is negligible (~120 ms one-way); bottleneck is ground ops.
# 8 hours confirmed by mission ops input.
T_DWELL_ANALYSIS = 8.0 * 3600.0    # s   sensor scan + ground relay cycle (confirmed)
T_HOLD_RH        = 0.5 * 3600.0    # s   MS stationkeeping at RH between tug arrivals
                                    #     (actual wait time comes from optimizer tug
                                    #     spiral timing, not a fixed value)
# Arm extension is mechanically limited (deploy the boom + grapple-lock onto the
# tug), NOT momentum-limited — the RCS impulse is a couple of seconds but the
# physical deploy + capture-lock takes ~1.5 min, which sets the phase duration.
T_ARM_EXTEND     = 90.0            # s   arm deploy + lock-on-to-tug wall-clock time

# ── Risk margins ──────────────────────────────────────────────────────────────
ABORT_MARGIN = 0.50    #  +50% abort/retry margin on debris RPO (risk TR-MT-08)

__all__ = [
    'N_MEAN', 'R_CUTOFF', 'INSPECT_STANDOFF', 'H_BAR_STANDOFF',
    'ARM_LEN_EXTENDED', 'DEBRIS_PANEL_SPAN_WC', 'KOS1_MARGIN', 'R_KOS2', 'R_KOS2_MEET', 'MEV_LOCK_DIST',
    'DOCK_PORT_HEIGHT', 'R_KOS2_DOCK', 'R_KOS1_DOCK',
    'R_SPIN_STANDOFF', 'APPROACH_CONE_DEG', 'R_RETREAT',
    'V_MAX', 'F_RCS', 'N_RCS_TOTAL', 'N_RCS_PER_EDGE', 'RCS_THR', 'RCS_ARM', 'DW_SYNC_DPS',
    'P_SOL', 'C_R', 'AOM_MS', 'A_TRIAX',
    'OMEGA_REQ_DPS', 'TUMBLE_SWEEP_DPS',
    'M_DEBRIS_WC', 'M_TUG', 'M_TUGDEB_WC', 'M_COMBINED_RH',
    'MS_DIMS', 'TUG_DIMS', 'DEBRIS_BUS_DIMS', 'DEBRIS_PANEL_FRAC',
    'T_DWELL_ANALYSIS', 'T_HOLD_RH', 'T_ARM_EXTEND', 'ABORT_MARGIN',
]
