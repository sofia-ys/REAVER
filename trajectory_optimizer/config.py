"""
REAVER Mission Configuration
=============================
Single source of truth for all mission, spacecraft, and orbital parameters.
Edit here; reaver_optimizer.py, rh_raan_sweep.py, and nominal_mission.py
all import from this module via  `from config import *`.
"""
import numpy as np

# ── Physical constants ────────────────────────────────────────────────────────
MU  = 3.986004418e14   # m³/s²  Earth gravitational parameter
G0  = 9.80665          # m/s²   standard gravity
DAY = 86400.0          # s      seconds per day
D2R = np.pi / 180.0   # rad/°  degrees to radians

# ── Mothership ────────────────────────────────────────────────────────────────
MS_DRY        = 956      # kg   dry mass
MS_ISP        = 253.0    # s    specific impulse (monoprop chemical)
MS_VEX        = MS_ISP * G0   # m/s  effective exhaust velocity
MS_RCS_MARGIN = 0.10     # -    10 % of orbital prop reserved for RCS / proximity ops

# ── Tug ───────────────────────────────────────────────────────────────────────
TUG_DRY = 221.4   # kg   dry mass per tug
TUG_ISP = 1600.0       # s    specific impulse (electric propulsion)
TUG_VEX = TUG_ISP * G0
TUG_THR = 0.07        # N    thrust per tug

# ── Mission ───────────────────────────────────────────────────────────────────
N_PHASE_REV = 18       # revolutions on phasing orbit to close 90° phase gap
T_OPS       = 5.0      # days  proximity operations time per debris capture
MAX_DAYS    = 365.0    # days  mission completion constraint
SOFT_MASS   = 2000.0   # kg   debris mass flag threshold for reporting

# ── Recycling Hub ─────────────────────────────────────────────────────────────
RH_SMA  = 42878.0e3    # m    SMA  (super-synchronous graveyard orbit)
RH_INC  = 7.0          # deg  inclination
RH_RAAN = 64.0         # deg  RAAN  (optimal value from rh_raan_sweep.py)

# ── RPO ΔV constants (close-range proximity ops) ─────────────────────────────
# Produced by the RPO simulation (run:  python trajectory_optimizer/rpo_run.py).
# Kept split: each has an independent physical driver and is updated separately.
DV_RPO_DEBRIS   = 0.65   # m/s  MS alone, COMSATBW-1 @ 1 rpm, V-bar inspect + tumble-axis capture (+50% abort)
DV_RPO_DETUMBLE = 1.33   # m/s  combined body momentum dump (AOCS sizing, REQ-ACS-M4)
DV_RPO_TUG_MEET = 0.14   # m/s  MS alone, cooperative close approach to tug+debris (30 m -> KOS2 = arm)
DV_RPO_RH_DOCK  = 0.27   # m/s  MS+tug+debris dock to RH via top port (KOS spheres, CoM offset)
DV_RPO_RH       = DV_RPO_TUG_MEET + DV_RPO_RH_DOCK   # m/s  total per RH cycle (tug-meet + RH-dock)

# ── Debris catalogue ──────────────────────────────────────────────────────────
# Columns:
#   norad_id, name, mass_kg,
#   sma_km, inc_deg, raan_deg,
#   eccentricity, arg_of_perigee_deg, mean_anomaly_deg,
#   rev_at_epoch, period_min, apoapsis_km, periapsis_km,
#   epoch (ISO 8601 UTC), launch_date
_RAW = [(443, 'ECHOSTAR 1', 1902.87, 42544.807, 8.5674, 63.1596), (488, 'AMC-3 (GE-3)', 1585.73, 42165.071, 8.092, 64.3488), (505, 'NSS 806 (INTELSAT 806)', 2135.78, 42502.032, 7.0712, 69.9995), (513, 'INTELSAT 805', 1932.61, 42388.082, 7.8332, 66.0575), (579, 'AMC-7 (GE-7)', 1935.0, 42486.971, 6.8592, 70.5256), (590, 'AMC-8 (GE-8)', 2015.0, 42482.105, 6.0647, 72.8264), (628, 'EUTE 12 WEST A (AB 1)', 2700.0, 42728.199, 7.2278, 68.7938), (640, 'GALAXY 12', 1760.0, 42528.577, 6.9338, 69.5405), (663, 'AMC-10 (GE-10)', 2315.0, 42490.063, 6.6355, 69.276), (672, 'AMC-11 (GE-11)', 2315.0, 42499.585, 0.2335, 80.2545), (697, 'GALAXY 14', 2087.0, 42549.012, 4.9754, 75.3425), (699, 'GALAXY 15', 2033.0, 42542.743, 3.481, 77.4994), (731, 'AMC-18', 2081.0, 42565.455, 2.5582, 84.3917), (748, 'INTELSAT 11 (PAS 11)', 2450.0, 42510.914, 3.5105, 78.9272), (755, 'HORIZONS 2', 2304.0, 42164.58, 2.5479, 80.8867), (774, 'AMC-21', 2473.0, 42522.009, 2.8956, 80.3233), (786, 'NSS 9', 2290.0, 42164.923, 1.2073, 82.8539), (804, 'COMSATBW-1', 2440.0, 42164.245, 0.0696, 89.6878), (834, 'HYLAS 1', 2570.0, 42164.346, 5.5899, 73.2864), (884, 'METEOSAT 10 (MSG 3)', 2035.0, 42164.226, 4.6216, 61.1229), (967, 'METEOSAT 11 (MSG 4)', 2043.0, 42164.271, 3.1471, 71.2262)]

N_DEB = len(_RAW)

# ── Basic arrays (used by optimizer and transfer table) ──────────────────────
IDS   = np.array([r[0] for r in _RAW])
NAMES = [r[1] for r in _RAW]
MASS  = np.array([r[2] for r in _RAW])           # kg
SMA   = np.array([r[3] for r in _RAW]) * 1e3     # m
INC   = np.array([r[4] for r in _RAW])           # deg
RAAN  = np.array([r[5] for r in _RAW])           # deg

# ── Extended orbital elements (needed for visualisation) ─────────────────────
#ECC   = np.array([r[6]  for r in _RAW])          # eccentricity
#AOP   = np.array([r[7]  for r in _RAW])          # argument of perigee [deg]
#MA    = np.array([r[8]  for r in _RAW])          # mean anomaly at epoch [deg]
#PERIOD    = np.array([r[10] for r in _RAW])      # orbital period [min]
# APOAPSIS  = np.array([r[11] for r in _RAW])      # apoapsis altitude [km]
# PERIAPSIS = np.array([r[12] for r in _RAW])      # periapsis altitude [km]
# EPOCH       = [r[13] for r in _RAW]              # TLE epoch (ISO 8601 UTC)
# LAUNCH_DATE = [r[14] for r in _RAW]              # launch date strings

# ── Augment arrays with Recycling Hub as index N_DEB (= 16) ──────────────────
RH_IDX   = N_DEB
SMA_ALL  = np.append(SMA,  RH_SMA)
INC_ALL  = np.append(INC,  RH_INC)
RAAN_ALL = np.append(RAAN, RH_RAAN)
MASS_ALL = np.append(MASS, 0.0)
# ECC_ALL  = np.append(ECC,  0.0)    # RH assumed circular
# AOP_ALL  = np.append(AOP,  0.0)
# MA_ALL   = np.append(MA,   0.0)

# Explicit export list so that  `from config import *`  includes _RAW
__all__ = [
    'MU', 'G0', 'DAY', 'D2R',
    'MS_DRY', 'MS_ISP', 'MS_VEX', 'MS_RCS_MARGIN',
    'TUG_DRY', 'TUG_ISP', 'TUG_VEX', 'TUG_THR',
    'N_PHASE_REV', 'T_OPS', 'MAX_DAYS', 'SOFT_MASS',
    'RH_SMA', 'RH_INC', 'RH_RAAN',
    'DV_RPO_DEBRIS', 'DV_RPO_DETUMBLE', 'DV_RPO_TUG_MEET',
    'DV_RPO_RH_DOCK', 'DV_RPO_RH',
    '_RAW', 'N_DEB', 'IDS', 'NAMES', 'MASS', 'SMA', 'INC', 'RAAN',
    'RH_IDX', 'SMA_ALL', 'INC_ALL', 'RAAN_ALL', 'MASS_ALL',
    'MS_THR', 'MS_BURN_S', 'MS_WET_ESTIMATE',
]

# 'ECC', 'AOP', 'MA', 'PERIOD', 'APOAPSIS', 'PERIAPSIS', 'EPOCH', 'LAUNCH_DATE'
# 'ECC_ALL', 'AOP_ALL', 'MA_ALL',

# ── Mothership finite-burn model ──────────────────────────────────────────────
MS_THR   = 64.0        # N, thruster force
MS_BURN_S = 38 * 60    # s, max single firing duration
MS_WET_ESTIMATE = 4000.0  # kg, representative mid-mission mass for table pre-computation