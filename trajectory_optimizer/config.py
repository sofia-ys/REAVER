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
TUG_DRY = 284.36796    # kg   dry mass per tug
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
# Columns:
#   norad_id, name, mass_kg,
#   sma_km, inc_deg, raan_deg,
#   eccentricity, arg_of_perigee_deg, mean_anomaly_deg,
#   rev_at_epoch, period_min, apoapsis_km, periapsis_km,
#   epoch (ISO 8601 UTC), launch_date
_RAW = [
    (489, 'NSS 5 (Intelsat 803)',    2060.46, 42494.423, 10.4813,  48.1496,
     0.00037807, 252.7691, 298.2382, 10473, 1452.973, 36132.354, 36100.222,
     '2026-05-26T00:14:15.561312', '1997-09-23'),
    (495, 'Sirius 2 (GE-1E)',        1762.14, 42410.143, 11.9727,  34.1101,
     0.00034646,  39.4548, 144.5444, 10390, 1448.653, 36046.701, 36017.314,
     '2026-05-26T03:20:08.176992', '1997-11-11'),
    (518, 'EUTE W2',                  1793.86, 42451.133, 11.6171,  38.4600,
     0.00022557, 215.5396, 321.0956, 10063, 1450.754, 36082.573, 36063.422,
     '2026-05-25T16:08:47.570784', '1998-10-05'),
    (530, 'Skynet 4E',                1490.00, 42516.912, 12.0550,   9.2220,
     0.00015206, 283.6254, 277.8551,  5093, 1454.127, 36145.242, 36132.312,
     '2026-05-26T01:15:38.637792', '1999-02-26'),
    (556, 'EUTE 16C (SESAT 1)',       2600.00, 42537.778, 10.8393,  46.2495,
     0.00032070, 148.6955,  26.4868,  9512, 1455.197, 36173.285, 36146.001,
     '2026-05-25T20:58:02.519328', '2000-04-17'),
    (593, 'Skynet 4F',                1489.00, 42482.488, 12.1030,  18.6123,
     0.00042097, 227.2613, 322.7663,  9239, 1452.361, 36122.237, 36086.470,
     '2026-05-26T07:47:52.035648', '2001-02-07'),
    (628, 'EUTE 12 West A (AB 1)',    2700.00, 42728.199,  7.2278,  68.7938,
     0.00074117,  13.7566,  94.6286,  8643, 1464.980, 36381.733, 36318.395,
     '2026-05-26T07:08:07.361952', '2002-08-28'),
    (653, 'EUTE 33A (EB 3)',          1552.00, 42558.620,  9.6938,  55.8840,
     0.00033649, 209.4150, 171.6507,  4484, 1456.267, 36194.806, 36166.165,
     '2026-05-25T20:19:13.590912', '2003-09-27'),
    (660, 'SESAT 2 (Express AM-22)',  2542.00, 42420.904,  8.2239,  64.3413,
     0.00092447, 134.6656,  56.8443,  8181, 1449.204, 36081.986, 36003.552,
     '2026-05-26T01:57:38.623392', '2003-12-28'),
    (705, 'Meteosat 9 (MSG 2)',       2054.00, 42163.863,  9.3136,  54.3880,
     0.00015392, 326.1315, 216.2429,   677, 1436.053, 35792.218, 35779.238,
     '2026-05-25T20:30:38.697984', '2005-12-21'),
    (741, 'EUTE 8A (Sinosat 3)',      2320.00, 42716.964,  9.5162,  48.3713,
     0.00144023, 102.0507,  79.4411,  4121, 1464.402, 36400.352, 36277.307,
     '2026-05-25T17:48:38.858688', '2007-05-31'),
    (804, 'COMSATBW-1',              2440.00, 42164.245,  0.0696,  89.6878,
     0.00033791, 340.9467, 321.1036,  6108, 1436.072, 35800.358, 35771.862,
     '2026-05-26T05:39:38.566656', '2009-09-30'),
    (819, 'COMSATBW-2',              2440.00, 42164.548,  0.0344, 111.5455,
     0.00029020, 325.9134, 261.5690,  5873, 1436.087, 35798.649, 35774.177,
     '2026-05-26T05:28:18.863040', '2010-05-21'),
    (834, 'HYLAS 1',                  2570.00, 42164.346,  5.5899,  73.2864,
     0.00013473, 314.0777, 255.5018,   666, 1436.077, 35791.891, 35780.530,
     '2026-05-26T01:41:38.730624', '2010-11-26'),
    (884, 'Meteosat 10 (MSG 3)',      2035.00, 42164.226,  4.6216,  61.1229,
     0.00019402, 327.6195, 244.7973,  5063, 1436.071, 35794.271, 35777.910,
     '2026-05-26T01:59:15.112320', '2012-07-05'),
    (967, 'Meteosat 11 (MSG 4)',      2043.00, 42164.271,  3.1471,  71.2262,
     0.00016638,   0.2149, 263.4215,   699, 1436.073, 35793.151, 35779.120,
     '2026-05-26T05:27:11.471040', '2015-07-15'),
]

N_DEB = len(_RAW)

# ── Basic arrays (used by optimizer and transfer table) ──────────────────────
IDS   = np.array([r[0] for r in _RAW])
NAMES = [r[1] for r in _RAW]
MASS  = np.array([r[2] for r in _RAW])           # kg
SMA   = np.array([r[3] for r in _RAW]) * 1e3     # m
INC   = np.array([r[4] for r in _RAW])           # deg
RAAN  = np.array([r[5] for r in _RAW])           # deg

# ── Extended orbital elements (needed for visualisation) ─────────────────────
ECC   = np.array([r[6]  for r in _RAW])          # eccentricity
AOP   = np.array([r[7]  for r in _RAW])          # argument of perigee [deg]
MA    = np.array([r[8]  for r in _RAW])          # mean anomaly at epoch [deg]
PERIOD    = np.array([r[10] for r in _RAW])      # orbital period [min]
APOAPSIS  = np.array([r[11] for r in _RAW])      # apoapsis altitude [km]
PERIAPSIS = np.array([r[12] for r in _RAW])      # periapsis altitude [km]
EPOCH       = [r[13] for r in _RAW]              # TLE epoch (ISO 8601 UTC)
LAUNCH_DATE = [r[14] for r in _RAW]              # launch date strings

# ── Augment arrays with Recycling Hub as index N_DEB (= 16) ──────────────────
RH_IDX   = N_DEB
SMA_ALL  = np.append(SMA,  RH_SMA)
INC_ALL  = np.append(INC,  RH_INC)
RAAN_ALL = np.append(RAAN, RH_RAAN)
MASS_ALL = np.append(MASS, 0.0)
ECC_ALL  = np.append(ECC,  0.0)    # RH assumed circular
AOP_ALL  = np.append(AOP,  0.0)
MA_ALL   = np.append(MA,   0.0)

# Explicit export list so that  `from config import *`  includes _RAW
__all__ = [
    'MU', 'G0', 'DAY', 'D2R',
    'MS_DRY', 'MS_ISP', 'MS_VEX',
    'TUG_DRY', 'TUG_ISP', 'TUG_VEX', 'TUG_THR',
    'N_PHASE_REV', 'T_OPS', 'MAX_DAYS', 'SOFT_MASS',
    'RH_SMA', 'RH_INC', 'RH_RAAN',
    '_RAW', 'N_DEB', 'IDS', 'NAMES', 'MASS', 'SMA', 'INC', 'RAAN',
    'ECC', 'AOP', 'MA', 'PERIOD', 'APOAPSIS', 'PERIAPSIS', 'EPOCH', 'LAUNCH_DATE',
    'RH_IDX', 'SMA_ALL', 'INC_ALL', 'RAAN_ALL', 'MASS_ALL',
    'ECC_ALL', 'AOP_ALL', 'MA_ALL',
]
