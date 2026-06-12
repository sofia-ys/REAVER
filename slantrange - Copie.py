"""
REAVER Mission — Worst-Case Slant Range Calculator
===================================================
Computes slant range and elevation angle for each spacecraft × Estrack
ground station combination, propagating the inclined orbit over one full
figure-8 cycle to find the true worst case (longest range, lowest elevation).

Spacecraft catalogue includes:
  - All 16 debris targets (from REAVER trajectory optimiser notes, 2026-05-26)
  - Recycling Hub orbit (36,500 km altitude, 7° inc, RAAN = 35°, circular)
    Source: REAVER trajectory optimiser astrodynamics notes 03/06/2026

Estrack stations: New Norcia (NNO), Cebreros (CEB), Kourou (KRU).

Method
------
For an inclined near-GEO orbit the sub-satellite point traces a figure-8
(analemma). We propagate the argument of latitude u = ω + M over one
sidereal day (the figure-8 period) in fine time steps, compute the
instantaneous sub-satellite lat/lon, then compute the 3-D slant range and
elevation angle to each ground station using the dot-product formula.

No external libraries beyond numpy are required.
"""

import numpy as np

# ── Constants ─────────────────────────────────────────────────────────────────
R_E   = 6371.0          # Earth mean radius, km
DEG   = np.pi / 180.0   # deg → rad
MIN_EL = 5.0            # Estrack minimum operational elevation, deg

# ── Estrack ground stations (lat °N, lon °E) ──────────────────────────────────
STATIONS = {
    "New Norcia (NNO)": {"lat": -31.048, "lon": 116.191},
    "Cebreros  (CEB)":  {"lat":  40.452, "lon":  -4.367},
    "Kourou    (KRU)":  {"lat":   5.252, "lon": -52.805},
}

# ── REAVER spacecraft catalogue ───────────────────────────────────────────────
# Fields: (IDX, name, mass_kg, sma_km, inc_deg, raan_deg, ecc, omega_deg, M_deg)
# Source: REAVER trajectory optimiser notes, epoch 2026-05-26
#
# Recycling Hub orbital elements:
#   altitude = 36,500 km  →  sma = 6371 + 36500 = 42871 km
#   inc = 7°, RAAN = 35° (optimal per astrodynamics notes 03/06/2026)
#   circular orbit → ecc = 0, omega = 0, M sweeps full 360°
SPACECRAFT = [(443, 'ECHOSTAR 1', 1902.87, 42544.807, 8.5674, 63.1596, 0.00027973, 335.7476, 221.2268), (488, 'AMC-3 (GE-3)', 1585.73, 42165.071, 8.092, 64.3488, 0.00030808, 1.0184, 218.1718), (505, 'NSS 806 (INTELSAT 806)', 2135.78, 42502.032, 7.0712, 69.9995, 0.00036225, 87.6677, 23.711), (513, 'INTELSAT 805', 1932.61, 42388.082, 7.8332, 66.0575, 0.00100564, 167.5631, 53.8461), (579, 'AMC-7 (GE-7)', 1935.0, 42486.971, 6.8592, 70.5256, 0.00041903, 95.2883, 74.9744), (590, 'AMC-8 (GE-8)', 2015.0, 42482.105, 6.0647, 72.8264, 0.00053041, 261.2009, 291.2849), (628, 'EUTE 12 WEST A (AB 1)', 2700.0, 42728.199, 7.2278, 68.7938, 0.00074117, 13.7566, 94.6286), (640, 'GALAXY 12', 1760.0, 42528.577, 6.9338, 69.5405, 0.00047444, 303.0147, 246.7094), (663, 'AMC-10 (GE-10)', 2315.0, 42490.063, 6.6355, 69.276, 0.00018023, 45.0329, 155.2031), (672, 'AMC-11 (GE-11)', 2315.0, 42499.585, 0.2335, 80.2545, 0.00017361, 35.6248, 245.859), (697, 'GALAXY 14', 2087.0, 42549.012, 4.9754, 75.3425, 0.00046711, 337.5921, 265.2236), (699, 'GALAXY 15', 2033.0, 42542.743, 3.481, 77.4994, 0.00139388, 139.9163, 1.189), (731, 'AMC-18', 2081.0, 42565.455, 2.5582, 84.3917, 0.00054553, 121.0765, 4.3151), (748, 'INTELSAT 11 (PAS 11)', 2450.0, 42510.914, 3.5105, 78.9272, 8.535e-05, 77.9917, 39.1927), (755, 'HORIZONS 2', 2304.0, 42164.58, 2.5479, 80.8867, 0.00039736, 344.9402, 205.807), (774, 'AMC-21', 2473.0, 42522.009, 2.8956, 80.3233, 0.00126492, 356.7648, 130.5413), (786, 'NSS 9', 2290.0, 42164.923, 1.2073, 82.8539, 0.00022645, 330.7351, 122.6128), (804, 'COMSATBW-1', 2440.0, 42164.245, 0.0696, 89.6878, 0.00033791, 340.9467, 321.1036), (834, 'HYLAS 1', 2570.0, 42164.346, 5.5899, 73.2864, 0.00013473, 314.0777, 255.5018), (884, 'METEOSAT 10 (MSG 3)', 2035.0, 42164.226, 4.6216, 61.1229, 0.00019402, 327.6195, 244.7973), (967, 'METEOSAT 11 (MSG 4)', 2043.0, 42164.271, 3.1471, 71.2262, 0.00016638, 0.2149, 263.4215),
    # ── Recycling Hub ─────────────────────────────────────────────────────────
    ("RH", "Recycling Hub",              0.00, 42871.000,  7.0000,  35.0000, 0.00000000,   0.0000,   0.0000),
]

# ── Helper: satellite position vector (ECI-like, Earth-fixed for GEO) ─────────
def sat_position_eci(sma, inc, raan, ecc, omega, M, n_steps=720):
    """
    Propagate one full orbital revolution in n_steps steps.
    Returns array of unit position vectors (3, n_steps) in
    Earth-Centred Earth-Fixed frame (valid for slow-precessing GEO orbits
    over a single day — precession is negligible on this timescale).
    """
    inc_r  = inc  * DEG
    raan_r = raan * DEG
    om_r   = omega * DEG

    # Solve Kepler's equation E - e*sin(E) = M for each mean anomaly step
    M_arr = np.linspace(0, 2*np.pi, n_steps, endpoint=False) + M * DEG
    E = M_arr.copy()
    for _ in range(50):                          # Newton–Raphson
        E -= (E - ecc * np.sin(E) - M_arr) / (1 - ecc * np.cos(E))

    # True anomaly
    nu = 2 * np.arctan2(np.sqrt(1+ecc)*np.sin(E/2),
                        np.sqrt(1-ecc)*np.cos(E/2))

    # Radius
    r = sma * (1 - ecc * np.cos(E))             # km

    # Position in perifocal frame
    x_pf = r * np.cos(nu)
    y_pf = r * np.sin(nu)

    # Rotation matrices: perifocal → ECI (313 Euler: RAAN, i, omega)
    cos_O, sin_O = np.cos(raan_r), np.sin(raan_r)
    cos_i, sin_i = np.cos(inc_r),  np.sin(inc_r)
    cos_w, sin_w = np.cos(om_r),   np.sin(om_r)

    # Standard perifocal-to-ECI rotation (Curtis, Orbital Mechanics, eq 4.46)
    Xx =  cos_O*cos_w - sin_O*sin_w*cos_i
    Xy = -cos_O*sin_w - sin_O*cos_w*cos_i
    Yx =  sin_O*cos_w + cos_O*sin_w*cos_i
    Yy = -sin_O*sin_w + cos_O*cos_w*cos_i
    Zx =  sin_w*sin_i
    Zy =  cos_w*sin_i

    X = Xx*x_pf + Xy*y_pf
    Y = Yx*x_pf + Yy*y_pf
    Z = Zx*x_pf + Zy*y_pf

    return np.vstack([X, Y, Z]), r   # shape (3, n_steps), (n_steps,)

# ── Helper: ground station ECEF unit vector ────────────────────────────────────
def gs_vector(lat_deg, lon_deg):
    lat = lat_deg * DEG
    lon = lon_deg * DEG
    return R_E * np.array([np.cos(lat)*np.cos(lon),
                            np.cos(lat)*np.sin(lon),
                            np.sin(lat)])          # km

# ── Main computation ───────────────────────────────────────────────────────────
print(f"\n{'REAVER — Worst-Case Slant Range per Spacecraft × Estrack Station':^90}")
print(f"{'Minimum elevation threshold':^90}: {MIN_EL}°")
print("="*120)
print(f"{'Spacecraft':<28} {'Station':<22} {'Slant Range [km]':>18} {'Elev at worst [°]':>19} {'Above 5°?':>12}")
print("-"*120)

global_worst        = {"range": 0, "name": "", "station": "", "el": 0}
global_worst_debris = {"range": 0, "name": "", "station": "", "el": 0}
global_worst_rh     = {"range": 0, "name": "", "station": "", "el": 0}

for (idx, name, mass, sma, inc, raan, ecc, omega, M) in SPACECRAFT:

    # Print section header to separate debris from recycling hub
    if name == "Recycling Hub":
        print("-"*120)
        print(f"  {'── Recycling Hub ──'}")
        print("-"*120)

    pos_arr, r_arr = sat_position_eci(sma, inc, raan, ecc, omega, M)

    for st_name, st in STATIONS.items():
        g_vec = gs_vector(st["lat"], st["lon"])
        g_hat = g_vec / np.linalg.norm(g_vec)

        d_vec = pos_arr - g_vec[:, np.newaxis]
        d_mag = np.linalg.norm(d_vec, axis=0)

        sin_el = np.einsum('i,ij->j', g_hat, d_vec) / d_mag
        el_deg = np.degrees(np.arcsin(np.clip(sin_el, -1, 1)))

        mask = el_deg >= MIN_EL
        if mask.any():
            worst_idx = np.argmax(d_mag[mask])
            d_worst   = d_mag[mask][worst_idx]
            el_worst  = el_deg[mask][worst_idx]
            above     = "YES"
        else:
            d_worst  = d_mag[np.argmax(el_deg)]
            el_worst = el_deg.max()
            above    = "NO — never visible"

        print(f"{name:<28} {st_name:<22} {d_worst:>18.1f} {el_worst:>19.2f} {above:>12}")

        if above == "YES":
            if d_worst > global_worst["range"]:
                global_worst.update({"range": d_worst, "name": name,
                                     "station": st_name, "el": el_worst})
            if name == "Recycling Hub" and d_worst > global_worst_rh["range"]:
                global_worst_rh.update({"range": d_worst, "name": name,
                                        "station": st_name, "el": el_worst})
            if name != "Recycling Hub" and d_worst > global_worst_debris["range"]:
                global_worst_debris.update({"range": d_worst, "name": name,
                                            "station": st_name, "el": el_worst})

print("="*120)
print(f"\n▶  WORST-CASE DEBRIS (longest slant range above {MIN_EL}°):")
print(f"   Spacecraft : {global_worst_debris['name']}")
print(f"   Station    : {global_worst_debris['station']}")
print(f"   Range      : {global_worst_debris['range']:.1f} km")
print(f"   Elevation  : {global_worst_debris['el']:.2f}°")

print(f"\n▶  WORST-CASE RECYCLING HUB (longest slant range above {MIN_EL}°):")
print(f"   Spacecraft : {global_worst_rh['name']}")
print(f"   Station    : {global_worst_rh['station']}")
print(f"   Range      : {global_worst_rh['range']:.1f} km")
print(f"   Elevation  : {global_worst_rh['el']:.2f}°")

print(f"\n▶  OVERALL WORST-CASE across entire mission:")
print(f"   Spacecraft : {global_worst['name']}")
print(f"   Station    : {global_worst['station']}")
print(f"   Range      : {global_worst['range']:.1f} km")
print(f"   Elevation  : {global_worst['el']:.2f}°")
print(f"\n   → Use {global_worst['range']:.0f} km in your FSPL calculation.\n")