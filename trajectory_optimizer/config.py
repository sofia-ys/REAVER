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
MS_DRY = 1523.406      # kg   dry mass
MS_ISP = 253.0         # s    specific impulse (monoprop chemical)
MS_VEX = MS_ISP * G0  # m/s  effective exhaust velocity

# ── Tug ───────────────────────────────────────────────────────────────────────
TUG_DRY = 284.36796        # kg   dry mass per tug
TUG_ISP = 1600.0       # s    specific impulse (electric propulsion)
TUG_VEX = TUG_ISP * G0
TUG_THR = 0.065        # N    thrust per tug

# ── Mission ───────────────────────────────────────────────────────────────────
N_PHASE_REV = 15       # revolutions on phasing orbit to close 90° phase gap
T_OPS       = 10.0     # days  proximity operations time per debris capture
MAX_DAYS    = 365.0    # days  mission completion constraint
SOFT_MASS   = 2000.0   # kg   debris mass flag threshold for reporting

# ── Recycling Hub ─────────────────────────────────────────────────────────────
RH_SMA  = 42878.0e3    # m    SMA  (super-synchronous graveyard orbit)
RH_INC  = 7.0          # deg  inclination
RH_RAAN = 35.0         # deg  RAAN  (optimal value from rh_raan_sweep.py)

# ── Debris catalogue ──────────────────────────────────────────────────────────
#   (norad_id, name, mass_kg, sma_km, inc_deg, raan_deg)
_RAW = [
    (489, 'NSS 5 (Intelsat 803)',     2060.46, 42494.423, 10.4813,  48.1496),
    (495, 'Sirius 2 (GE-1E)',        1762.14, 42410.084, 11.9727,  34.1101),
    (518, 'EUTE W2',                  1793.86, 42451.133, 11.6171,  38.4600),
    (530, 'Skynet 4E',                1490.00, 42516.912, 12.0550,   9.2220),
    (556, 'EUTE 16C (SESAT 1)',       2600.00, 42537.778, 10.8393,  46.2495),
    (593, 'Skynet 4F',                1489.00, 42482.488, 12.1030,  18.6123),
    (628, 'EUTE 12 West A (AB 1)',    2700.00, 42728.199,  7.2278,  68.7938),
    (653, 'EUTE 33A (EB 3)',          1552.00, 42558.620,  9.6938,  55.8840),
    (660, 'SESAT 2 (Express AM-22)',  2542.00, 42420.904,  8.2239,  64.3413),
    (705, 'Meteosat 9 (MSG 2)',       2054.00, 42163.863,  9.3136,  54.3880),
    (741, 'EUTE 8A (Sinosat 3)',      2320.00, 42716.964,  9.5162,  48.3713),
    (804, 'COMSATBW-1',              2440.00, 42164.245,  0.0696,  89.6878),
    (819, 'COMSATBW-2',              2440.00, 42164.548,  0.0344, 111.5455),
    (834, 'HYLAS 1',                 2570.00, 42164.346,  5.5899,  73.2864),
    (884, 'Meteosat 10 (MSG 3)',      2035.00, 42164.226,  4.6216,  61.1229),
    (967, 'Meteosat 11 (MSG 4)',      2043.00, 42164.271,  3.1471,  71.2262),
]

N_DEB    = len(_RAW)
IDS      = np.array([r[0] for r in _RAW])
NAMES    = [r[1] for r in _RAW]
MASS     = np.array([r[2] for r in _RAW])          # kg
SMA      = np.array([r[3] for r in _RAW]) * 1e3    # m
INC      = np.array([r[4] for r in _RAW])          # deg
RAAN     = np.array([r[5] for r in _RAW])          # deg

# Augment arrays with Recycling Hub as index N_DEB (= 16)
RH_IDX   = N_DEB
SMA_ALL  = np.append(SMA,  RH_SMA)
INC_ALL  = np.append(INC,  RH_INC)
RAAN_ALL = np.append(RAAN, RH_RAAN)
MASS_ALL = np.append(MASS, 0.0)

# Explicit export list so that  `from config import *`  includes _RAW
__all__ = [
    'MU', 'G0', 'DAY', 'D2R',
    'MS_DRY', 'MS_ISP', 'MS_VEX',
    'TUG_DRY', 'TUG_ISP', 'TUG_VEX', 'TUG_THR',
    'N_PHASE_REV', 'T_OPS', 'MAX_DAYS', 'SOFT_MASS',
    'RH_SMA', 'RH_INC', 'RH_RAAN',
    '_RAW', 'N_DEB', 'IDS', 'NAMES', 'MASS', 'SMA', 'INC', 'RAAN',
    'RH_IDX', 'SMA_ALL', 'INC_ALL', 'RAAN_ALL', 'MASS_ALL',
]