"""
REAVER Mothership – Propellant Slosh, Tug Release, CG Shift & Torque Analysis
=============================================================================

Purpose
-------
Reduced-order engineering analysis to estimate the worst-case lateral centre of
mass / centre of gravity (CG) shift of the REAVER mothership and the resulting
main-burn torque disturbance.

This script combines:
    1. Chemical propellant depletion during the mothership burn sequence.
    2. Static lateral CG shift from unequal chemical tank draining.
    3. Conservative propellant slosh CG envelope using spherical tank geometry.
    4. Dynamic first-mode slosh response using a spring-mass-damper model.
    5. Unequal wet tug masses caused by tailored xenon loading.
    6. Tug release and circular rail rebalancing.

Modelling philosophy
--------------------
This is NOT CFD and is NOT a final flight slosh model. It is a preliminary,
report-ready engineering model inspired by NASA SP-106 style reduced-order
sloshing analysis. It is intended to bound the CG envelope and convert it into
an induced torque:

    tau = |Delta r_CG x T|

The most important output is therefore the maximum CG shift and torque, not a
perfect prediction of the liquid free-surface shape.

Author: REAVER DSE Team
Date: 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

# =============================================================================
# 1. GLOBAL USER INPUTS
# =============================================================================

# ----------------------------- Spacecraft ------------------------------------
M_DRY_KG = 1120.0                         # mothership dry mass without prop/tugs
CG_DRY_M = np.array([0.0, 0.0, 1.5])       # dry CG in body frame [m]

# ----------------------------- Chemical propulsion ---------------------------
G0 = 9.80665                              # standard gravity [m/s^2]
ISP_CHEM_S = 253.0                        # HPGP / LMP-103S assumed Isp [s]
MAIN_THRUST_N = 4 * 22.0                  # 4 x 22 N main thrusters [N]
THRUST_DIR = np.array([0.0, 0.0, 1.0])     # nominal thrust direction along +z

# Chemical propellant physically loaded in the mothership tanks.
# Use 1780 kg if including main transfer propellant, RCS/proximity allocation,
# residual and unusable allowance in the physical tank fill state.
M_CHEM_INITIAL_KG = 1550.0

# ----------------------------- Chemical tanks --------------------------------
N_TANKS = 4
RHO_LMP103S_KG_M3 = 1240.0
V_TANK_M3 = 0.3679                        # each selected spherical PMD tank
Z_TANK_M = 2.0                            # tank centre z-position [editable]

# Four tank centres in a square of side length 1 m, centred on spacecraft axis.
TANK_POSITIONS_M = np.array([
    [ 0.5,  0.5, Z_TANK_M],
    [-0.5,  0.5, Z_TANK_M],
    [-0.5, -0.5, Z_TANK_M],
    [ 0.5, -0.5, Z_TANK_M],
], dtype=float)

# Tank-drain imbalance cases. These are operational imbalance envelopes only.
# No single-tank failure is included because the aim is a simple worst-case
# operational CG envelope, not a fault-tree analysis.
TANK_IMBALANCE_CASES = {
    "equal draining": 0.00,
    "2% pair imbalance": 0.02,
    "5% pair imbalance": 0.05,
    "10% pair imbalance": 0.10,
}
DEFAULT_TANK_IMBALANCE_FOR_COMBINED = 0.05

# ----------------------------- PMD/slosh assumptions -------------------------
# Slosh reduction factor for a PMD tank. 1.0 means bare/free-surface tank;
# lower values represent reduced mobile liquid and centroid motion due to PMD.
PMD_CASES = {
    "no PMD / free surface": 1.00,
    "moderate PMD": 0.50,
    "strong PMD": 0.25,
}
DEFAULT_PMD_CASE = "moderate PMD"

# Dynamic first-mode slosh parameters.
ZETA_SLOSH = 0.02                         # damping ratio
OMEGA_REF_RAD_S = 1.5                     # fundamental slosh frequency at phi=0.5
ALPHA_FUNDAMENTAL = 0.45                  # participating mobile mass fraction
PENDULUM_ARM_FRACTION = 0.55              # effective slosh displacement scale

MODE_PARTICIPATION = {
    "fundamental lateral": 1.00,
    "rotary / non-planar": 0.30,
    "2nd lateral mode": 0.20,
    "vertical / axial": 0.05,
}

# ----------------------------- Tug system ------------------------------------
TUG_DRY_MASS_KG = 191.0
R_TUG_RING_M = 1.0                        # tugs lie on circular rail radius [m]
Z_TUG_RING_M = 3.0                        # tug rail z-position [editable]

# Target xenon requirements from the tug propellant margin table [kg].
# These are treated as final loaded Xe budgets by default.
TARGET_SEQUENCE = [
    {"target": "COMSATBW-1",        "xe_req_kg": 99.5},
    {"target": "METEOSAT 10",       "xe_req_kg": 29.0},
    {"target": "GALAXY 14",         "xe_req_kg": 28.7},
    {"target": "AMC-3",             "xe_req_kg": 10.9},
    {"target": "ECHOSTAR 1",        "xe_req_kg": 17.8},
]

# Tailored tug loading strategy:
# two tugs sized for the worst case, then one for 2nd, 3rd, and 4th worst case.
# This reduces the unused xenon relative to 5 x 99.5 kg.
TUG_XE_LOADS_KG = np.array([99.5, 99.5, 29.0, 28.7, 28.7], dtype=float)
TUG_WET_MASSES_KG = TUG_DRY_MASS_KG + TUG_XE_LOADS_KG

# Release assignment strategy.
# "cg_managed" releases the heaviest capable tug for each target. This usually
# reduces the final one-tug CG offset.
# "smallest_sufficient" releases the smallest capable tug, conserving larger
# tugs but often worsening late-campaign CG.
RELEASE_ASSIGNMENT_STRATEGY = "cg_managed"

# ----------------------------- Burn event sequence ---------------------------
# Use the optimiser delta-V sequence, but recompute chemical propellant per burn
# with the rocket equation because tailored tug loading changes total mass.
MISSION_EVENTS = [
    {"day": 0.0,   "event": "Depart Recycling Hub",       "type": "coast",   "dv": None,  "target": None},
    {"day": 0.0,   "event": "Burn 1 dep + plane chg",     "type": "burn",    "dv": 368.4, "target": None},
    {"day": 0.5,   "event": "Burn 2 circularise",         "type": "burn",    "dv": 12.9,  "target": None},
    {"day": 0.5,   "event": "Burn 3 phasing",             "type": "burn",    "dv": 28.9,  "target": None},
    {"day": 18.2,  "event": "Rendezvous: COMSATBW-1",     "type": "coast",   "dv": None,  "target": "COMSATBW-1"},

    {"day": 28.2,  "event": "Tug release: COMSATBW-1",    "type": "release", "dv": None,  "target": "COMSATBW-1"},
    {"day": 28.2,  "event": "Burn 1 dep + plane chg",     "type": "burn",    "dv": 244.7, "target": None},
    {"day": 28.7,  "event": "Burn 2 circularise",         "type": "burn",    "dv": 0.0,   "target": None},
    {"day": 28.7,  "event": "Burn 3 phasing",             "type": "burn",    "dv": 28.9,  "target": None},
    {"day": 46.4,  "event": "Rendezvous: METEOSAT 10",    "type": "coast",   "dv": None,  "target": "METEOSAT 10"},

    {"day": 56.4,  "event": "Tug release: METEOSAT 10",   "type": "release", "dv": None,  "target": "METEOSAT 10"},
    {"day": 56.4,  "event": "Burn 1 Hohmann dep",         "type": "burn",    "dv": 7.0,   "target": None},
    {"day": 56.9,  "event": "Burn 2 circ + plane chg",    "type": "burn",    "dv": 66.4,  "target": None},
    {"day": 56.9,  "event": "Burn 3 phasing",             "type": "burn",    "dv": 28.7,  "target": None},
    {"day": 74.9,  "event": "Rendezvous: GALAXY 14",      "type": "coast",   "dv": None,  "target": "GALAXY 14"},

    {"day": 84.9,  "event": "Tug release: GALAXY 14",     "type": "release", "dv": None,  "target": "GALAXY 14"},
    {"day": 84.9,  "event": "Burn 1 dep + plane chg",     "type": "burn",    "dv": 178.6, "target": None},
    {"day": 85.4,  "event": "Burn 2 circularise",         "type": "burn",    "dv": 7.0,   "target": None},
    {"day": 85.4,  "event": "Burn 3 phasing",             "type": "burn",    "dv": 28.9,  "target": None},
    {"day": 103.1, "event": "Rendezvous: AMC-3",          "type": "coast",   "dv": None,  "target": "AMC-3"},

    {"day": 113.1, "event": "Tug release: AMC-3",         "type": "release", "dv": None,  "target": "AMC-3"},
    {"day": 113.1, "event": "Burn 1 Hohmann dep",         "type": "burn",    "dv": 6.9,   "target": None},
    {"day": 113.6, "event": "Burn 2 circ + plane chg",    "type": "burn",    "dv": 27.8,  "target": None},
    {"day": 113.6, "event": "Burn 3 phasing",             "type": "burn",    "dv": 28.7,  "target": None},
    {"day": 131.5, "event": "Rendezvous: ECHOSTAR 1",     "type": "coast",   "dv": None,  "target": "ECHOSTAR 1"},

    {"day": 141.5, "event": "Tug release: ECHOSTAR 1",    "type": "release", "dv": None,  "target": "ECHOSTAR 1"},
    {"day": 141.5, "event": "Burn 1 Hohmann dep",         "type": "burn",    "dv": 6.0,   "target": None},
    {"day": 142.0, "event": "Burn 2 circ + plane chg",    "type": "burn",    "dv": 83.8,  "target": None},
    {"day": 142.0, "event": "Burn 3 phasing",             "type": "burn",    "dv": 28.6,  "target": None},
    {"day": 160.2, "event": "Arrive Recycling Hub",       "type": "coast",   "dv": None,  "target": None},
]

# =============================================================================
# 2. GEOMETRY AND MASS HELPER FUNCTIONS
# =============================================================================

def tank_radius_from_volume(volume_m3):
    """Spherical tank radius from volume."""
    return (3.0 * volume_m3 / (4.0 * np.pi)) ** (1.0 / 3.0)


R_TANK_M = tank_radius_from_volume(V_TANK_M3)


def fill_fraction_from_prop_mass(prop_mass_total_kg):
    """Mean fill fraction of the four chemical tanks from total propellant mass."""
    prop_per_tank = prop_mass_total_kg / N_TANKS
    return np.clip(prop_per_tank / (RHO_LMP103S_KG_M3 * V_TANK_M3), 0.0, 1.0)


def chemical_propellant_for_burn(mass_before_kg, delta_v_m_s):
    """Chemical propellant consumed in one burn from the rocket equation."""
    if delta_v_m_s is None or delta_v_m_s <= 0.0:
        return 0.0
    return mass_before_kg * (1.0 - np.exp(-delta_v_m_s / (ISP_CHEM_S * G0)))


def compute_total_cg(m_dry_kg, cg_dry_m, tank_masses_kg, tank_positions_m,
                     tug_masses_kg=None, tug_positions_m=None):
    """Mass-weighted spacecraft CG."""
    total_mass = float(m_dry_kg)
    moment = float(m_dry_kg) * np.array(cg_dry_m, dtype=float)

    for m, r in zip(tank_masses_kg, tank_positions_m):
        total_mass += float(m)
        moment += float(m) * np.array(r, dtype=float)

    if tug_masses_kg is not None and len(tug_masses_kg) > 0:
        for m, r in zip(tug_masses_kg, tug_positions_m):
            total_mass += float(m)
            moment += float(m) * np.array(r, dtype=float)

    return moment / total_mass, total_mass


# =============================================================================
# 3. SPHERICAL TANK LIQUID CENTROID AND SLOSH ENVELOPE
# =============================================================================

def spherical_segment_volume(h, R):
    """Volume of a spherical segment of height h measured from bottom."""
    return np.pi * h**2 * (R - h / 3.0)


def fill_height_from_fraction(phi, R, n_iter=60):
    """Solve fill height h for a spherical tank at fill fraction phi."""
    phi = np.clip(phi, 1e-6, 1.0 - 1e-6)
    V_sphere = 4.0 / 3.0 * np.pi * R**3
    V_target = phi * V_sphere
    h = 2.0 * R * phi

    for _ in range(n_iter):
        V_h = spherical_segment_volume(h, R)
        dVdh = np.pi * h * (2.0 * R - h)
        err = V_h - V_target
        if abs(err) < 1e-13:
            break
        h -= err / max(abs(dVdh), 1e-15) * np.sign(dVdh)
        h = np.clip(h, 1e-9, 2.0 * R - 1e-9)

    return h


def liquid_centroid_offset_from_tank_centre(phi, R):
    """
    Centroid offset of liquid in a partially filled spherical tank.

    Coordinate is along the effective acceleration direction. The returned
    magnitude is used as a conservative maximum displacement scale for lateral
    slosh. For phi=0.5, |offset| = 3R/8.
    """
    phi = np.clip(phi, 1e-6, 1.0 - 1e-6)
    h = fill_height_from_fraction(phi, R)
    V = spherical_segment_volume(h, R)

    # Integrate z*A(z) dz from bottom, where A(z)=pi*(2Rz-z^2)
    numerator = np.pi * ((2.0 * R * h**3) / 3.0 - h**4 / 4.0)
    z_from_bottom = numerator / V
    z_from_centre = z_from_bottom - R
    return z_from_centre


def slosh_static_cg_shift(prop_mass_total_kg, total_spacecraft_mass_kg,
                          pmd_factor=0.5):
    """
    Conservative lateral slosh CG envelope.

    All tanks are assumed to slosh coherently in the same lateral direction.
    PMD factor reduces the effective liquid centroid displacement/mobile mass.
    """
    phi = fill_fraction_from_prop_mass(prop_mass_total_kg)
    offset = abs(liquid_centroid_offset_from_tank_centre(phi, R_TANK_M))
    mobile_prop_mass = pmd_factor * prop_mass_total_kg
    delta_cg_m = mobile_prop_mass * offset / total_spacecraft_mass_kg
    return delta_cg_m, phi, offset


def tank_masses_with_pair_imbalance(prop_mass_total_kg, eps):
    """
    Operational unequal-drain case.

    Two tanks on +y side are high by eps and two tanks on -y side are low by eps.
    Total propellant mass is renormalised to remain constant.
    """
    base = prop_mass_total_kg / N_TANKS
    factors = np.array([1.0 + eps, 1.0 + eps, 1.0 - eps, 1.0 - eps])
    masses = base * factors
    masses *= prop_mass_total_kg / masses.sum()
    return masses


def tank_imbalance_cg_shift(prop_mass_total_kg, attached_tug_masses,
                            attached_tug_positions, eps):
    """Lateral CG shift due to paired unequal tank draining."""
    nominal_tank_masses = np.full(N_TANKS, prop_mass_total_kg / N_TANKS)
    imbalanced_tank_masses = tank_masses_with_pair_imbalance(prop_mass_total_kg, eps)

    cg_nom, _ = compute_total_cg(M_DRY_KG, CG_DRY_M,
                                 nominal_tank_masses, TANK_POSITIONS_M,
                                 attached_tug_masses, attached_tug_positions)
    cg_imb, total_mass = compute_total_cg(M_DRY_KG, CG_DRY_M,
                                           imbalanced_tank_masses, TANK_POSITIONS_M,
                                           attached_tug_masses, attached_tug_positions)
    dxy = np.linalg.norm((cg_imb - cg_nom)[:2])
    return dxy, total_mass


# =============================================================================
# 4. DYNAMIC FIRST-MODE SLOSH MODEL
# =============================================================================

def slosh_frequency(phi, R, omega_ref=OMEGA_REF_RAD_S):
    """Simple fill-dependent engineering fit normalised at phi=0.5."""
    phi = np.clip(phi, 0.02, 0.98)
    h = fill_height_from_fraction(phi, R)
    h_ref = fill_height_from_fraction(0.5, R)
    return omega_ref * np.sqrt(h / h_ref)


def slosh_ode(t, y, zeta, omega, a_func):
    """x'' + 2*zeta*omega*x' + omega^2*x = -a_lateral(t)."""
    x, xd = y
    xdd = -2.0 * zeta * omega * xd - omega**2 * x - a_func(t)
    return [xd, xdd]


def dynamic_slosh_response(phi, prop_mass_total_kg, total_mass_kg,
                           pmd_factor=0.5, excitation="pulse",
                           t_end=80.0):
    """Integrate 1-DOF dynamic lateral slosh response."""
    omega = slosh_frequency(phi, R_TANK_M)
    m_slosh = ALPHA_FUNDAMENTAL * pmd_factor * prop_mass_total_kg

    if excitation == "step":
        a0 = 1e-3
        def a_func(t):
            return a0
    elif excitation == "pulse":
        a0 = 5e-3
        pulse_duration = 0.5
        def a_func(t):
            return a0 if t <= pulse_duration else 0.0
    elif excitation == "sine":
        a0 = 2e-3
        f = 0.5
        def a_func(t):
            return a0 * np.sin(2.0 * np.pi * f * t)
    else:
        raise ValueError("excitation must be 'step', 'pulse', or 'sine'")

    t_eval = np.linspace(0.0, t_end, 2000)
    sol = solve_ivp(
        slosh_ode,
        [0.0, t_end],
        [0.0, 0.0],
        args=(ZETA_SLOSH, omega, a_func),
        t_eval=t_eval,
        rtol=1e-8,
        atol=1e-10,
    )

    x = sol.y[0]
    delta_cg = m_slosh * x / total_mass_kg
    torque = np.abs(delta_cg) * MAIN_THRUST_N

    return {
        "t_s": sol.t,
        "x_slosh_m": x,
        "delta_cg_m": delta_cg,
        "torque_Nm": torque,
        "omega_rad_s": omega,
        "m_slosh_kg": m_slosh,
    }


def mode_criticality(phi, prop_mass_total_kg, total_mass_kg, pmd_factor=0.5):
    """Rank simplified slosh modes by quasi-static torque contribution."""
    offset = abs(liquid_centroid_offset_from_tank_centre(phi, R_TANK_M))
    rows = []
    for mode, factor in MODE_PARTICIPATION.items():
        m_eff = factor * ALPHA_FUNDAMENTAL * pmd_factor * prop_mass_total_kg
        dcg = m_eff * offset / total_mass_kg
        tau = dcg * MAIN_THRUST_N
        rows.append({
            "mode": mode,
            "m_eff_kg": m_eff,
            "dcg_mm": dcg * 1e3,
            "torque_Nm": tau,
        })
    rows = sorted(rows, key=lambda r: r["torque_Nm"], reverse=True)
    return rows


# =============================================================================
# 5. TUG LOADING, RAIL BALANCING AND RELEASE MODEL
# =============================================================================

def tug_positions_from_angles(angles_rad, radius=R_TUG_RING_M, z=Z_TUG_RING_M):
    return np.array([
        [radius * np.cos(th), radius * np.sin(th), z]
        for th in angles_rad
    ], dtype=float)


def tug_lateral_mass_moment(masses_kg, angles_rad, radius=R_TUG_RING_M):
    masses_kg = np.array(masses_kg, dtype=float)
    angles_rad = np.array(angles_rad, dtype=float)
    return radius * np.array([
        np.sum(masses_kg * np.cos(angles_rad)),
        np.sum(masses_kg * np.sin(angles_rad)),
    ])


def optimise_tug_angles(masses_kg, n_starts=12):
    """
    Optimise tug angular locations on the ring to minimise lateral mass moment.
    First angle is fixed to zero to remove rotational degeneracy.
    """
    masses_kg = np.array(masses_kg, dtype=float)
    n = len(masses_kg)

    if n == 0:
        return np.array([]), 0.0
    if n == 1:
        return np.array([0.0]), masses_kg[0] * R_TUG_RING_M

    def objective(theta_free):
        theta = np.concatenate(([0.0], theta_free))
        moment = tug_lateral_mass_moment(masses_kg, theta)
        return np.dot(moment, moment)

    best = None
    rng = np.random.default_rng(7)

    # Try one evenly-spaced initial guess plus random starts.
    starts = [np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)[1:]]
    starts += [rng.uniform(0.0, 2.0 * np.pi, n - 1) for _ in range(n_starts)]

    for x0 in starts:
        sol = minimize(
            objective,
            x0,
            method="Nelder-Mead",
            options={"maxiter": 5000, "xatol": 1e-12, "fatol": 1e-12},
        )
        if best is None or sol.fun < best.fun:
            best = sol

    angles = np.mod(np.concatenate(([0.0], best.x)), 2.0 * np.pi)
    residual_moment_kgm = np.sqrt(best.fun)
    return angles, residual_moment_kgm


def lateral_cg_from_tugs(tug_masses_kg, tug_positions_m, total_mass_kg):
    """Lateral CG contribution from attached tugs only."""
    if tug_masses_kg is None or len(tug_masses_kg) == 0:
        return 0.0
    masses = np.array(tug_masses_kg, dtype=float)
    positions = np.array(tug_positions_m, dtype=float)
    moment_xy = np.sum(masses[:, None] * positions[:, :2], axis=0)
    return np.linalg.norm(moment_xy / total_mass_kg)


def assign_tugs_to_targets(tug_xe_loads_kg, target_sequence,
                           strategy=RELEASE_ASSIGNMENT_STRATEGY):
    """
    Assign tug IDs to target sequence.

    cg_managed: release the heaviest capable tug first, reducing late-campaign
    one-tug CG offsets.

    smallest_sufficient: release the smallest capable tug for each target.
    """
    remaining = list(range(len(tug_xe_loads_kg)))
    assignment = []

    for tgt in target_sequence:
        req = tgt["xe_req_kg"]
        feasible = [i for i in remaining if tug_xe_loads_kg[i] + 1e-9 >= req]
        if not feasible:
            raise RuntimeError(f"No remaining tug can satisfy {tgt['target']} ({req:.1f} kg Xe)")

        if strategy == "cg_managed":
            chosen = max(feasible, key=lambda i: tug_xe_loads_kg[i])
        elif strategy == "smallest_sufficient":
            chosen = min(feasible, key=lambda i: tug_xe_loads_kg[i])
        else:
            raise ValueError("Unknown tug assignment strategy")

        remaining.remove(chosen)
        assignment.append({
            "target": tgt["target"],
            "xe_req_kg": req,
            "tug_id": chosen,
            "tug_xe_load_kg": tug_xe_loads_kg[chosen],
        })

    return assignment


_TUG_BALANCE_CACHE = {}

def current_balanced_tug_state(attached_ids, tug_wet_masses_kg):
    """Return masses, optimised positions, angles and residual rail imbalance.

    Results are cached because the same attached-tug sets are evaluated many
    times during the mission timeline.
    """
    if len(attached_ids) == 0:
        return np.array([]), np.empty((0, 3)), np.array([]), 0.0

    key = tuple(attached_ids)
    masses = np.array([tug_wet_masses_kg[i] for i in attached_ids])

    if key in _TUG_BALANCE_CACHE:
        angles, residual = _TUG_BALANCE_CACHE[key]
    else:
        angles, residual = optimise_tug_angles(masses)
        _TUG_BALANCE_CACHE[key] = (angles, residual)

    positions = tug_positions_from_angles(angles)
    return masses, positions, angles, residual


# =============================================================================
# 6. MISSION SIMULATION
# =============================================================================

def evaluate_worst_cg_components(prop_remaining_kg, attached_masses_kg,
                                  attached_positions_m,
                                  total_mass_kg,
                                  pmd_case=DEFAULT_PMD_CASE,
                                  tank_imbalance_eps=DEFAULT_TANK_IMBALANCE_FOR_COMBINED):
    """
    Compute simple conservative worst-case CG components and torque.

    Components are combined by magnitude for a worst-direction envelope:
        tug CG + tank imbalance CG + static slosh CG.
    """
    pmd_factor = PMD_CASES[pmd_case]

    tug_cg = lateral_cg_from_tugs(attached_masses_kg, attached_positions_m, total_mass_kg)

    tank_cg, _ = tank_imbalance_cg_shift(
        prop_remaining_kg,
        attached_masses_kg,
        attached_positions_m,
        tank_imbalance_eps,
    )

    slosh_cg, phi, slosh_offset = slosh_static_cg_shift(
        prop_remaining_kg,
        total_mass_kg,
        pmd_factor=pmd_factor,
    )

    combined_cg = tug_cg + tank_cg + slosh_cg
    combined_tau = combined_cg * MAIN_THRUST_N

    return {
        "phi": phi,
        "tug_cg_m": tug_cg,
        "tank_imb_cg_m": tank_cg,
        "slosh_cg_m": slosh_cg,
        "slosh_offset_m": slosh_offset,
        "combined_cg_m": combined_cg,
        "combined_tau_Nm": combined_tau,
    }


def simulate_campaign():
    """Simulate burn-by-burn mass, tug release/rebalance and CG envelope."""
    assignment = assign_tugs_to_targets(TUG_XE_LOADS_KG, TARGET_SEQUENCE)
    target_to_tug = {a["target"]: a["tug_id"] for a in assignment}

    attached = list(range(len(TUG_WET_MASSES_KG)))
    prop_remaining = M_CHEM_INITIAL_KG
    rows = []

    # Initial balanced positions.
    attached_masses, attached_positions, angles, residual = current_balanced_tug_state(
        attached, TUG_WET_MASSES_KG
    )

    for ev in MISSION_EVENTS:
        # Always start each event from a balanced rail state if tugs remain.
        attached_masses, attached_positions, angles, residual_before = current_balanced_tug_state(
            attached, TUG_WET_MASSES_KG
        )
        total_mass_before = M_DRY_KG + prop_remaining + np.sum(attached_masses)

        row = {
            "day": ev["day"],
            "event": ev["event"],
            "type": ev["type"],
            "dv_m_s": ev["dv"],
            "target": ev["target"],
            "attached_before": attached.copy(),
            "prop_before_kg": prop_remaining,
            "total_mass_before_kg": total_mass_before,
            "released_tug_id": None,
            "released_mass_kg": 0.0,
            "prop_used_kg": 0.0,
            "cg_release_jump_mm": 0.0,
            "tau_release_jump_Nm": 0.0,
            "rail_residual_before_kgm": residual_before,
            "rail_residual_after_kgm": residual_before,
        }

        if ev["type"] == "burn":
            prop_used = chemical_propellant_for_burn(total_mass_before, ev["dv"])
            prop_used = min(prop_used, prop_remaining)
            prop_remaining -= prop_used
            row["prop_used_kg"] = prop_used

        elif ev["type"] == "release":
            release_id = target_to_tug[ev["target"]]
            if release_id not in attached:
                raise RuntimeError(f"Tug {release_id+1} for {ev['target']} already released")

            # Immediate post-release state: remaining tugs still at old positions.
            local_idx = attached.index(release_id)
            released_mass = TUG_WET_MASSES_KG[release_id]
            remaining_attached = attached.copy()
            remaining_attached.remove(release_id)

            remaining_masses_old = np.delete(attached_masses, local_idx)
            remaining_positions_old = np.delete(attached_positions, local_idx, axis=0)
            total_mass_after_release = M_DRY_KG + prop_remaining + np.sum(remaining_masses_old)

            cg_jump = lateral_cg_from_tugs(
                remaining_masses_old,
                remaining_positions_old,
                total_mass_after_release,
            )

            # Rebalance rail after release.
            attached = remaining_attached
            new_masses, new_positions, new_angles, residual_after = current_balanced_tug_state(
                attached, TUG_WET_MASSES_KG
            )

            row["released_tug_id"] = release_id
            row["released_mass_kg"] = released_mass
            row["cg_release_jump_mm"] = cg_jump * 1e3
            row["tau_release_jump_Nm"] = cg_jump * MAIN_THRUST_N
            row["rail_residual_after_kgm"] = residual_after

        elif ev["type"] == "coast":
            pass
        else:
            raise ValueError(f"Unknown event type {ev['type']}")

        # Current balanced post-event state.
        attached_masses_after, attached_positions_after, _, residual_after_now = current_balanced_tug_state(
            attached, TUG_WET_MASSES_KG
        )
        total_mass_after = M_DRY_KG + prop_remaining + np.sum(attached_masses_after)

        components = evaluate_worst_cg_components(
            prop_remaining,
            attached_masses_after,
            attached_positions_after,
            total_mass_after,
        )

        row.update({
            "attached_after": attached.copy(),
            "prop_after_kg": prop_remaining,
            "total_mass_after_kg": total_mass_after,
            "fill_fraction": components["phi"],
            "tug_cg_balanced_mm": components["tug_cg_m"] * 1e3,
            "tank_imb_cg_mm": components["tank_imb_cg_m"] * 1e3,
            "slosh_cg_mm": components["slosh_cg_m"] * 1e3,
            "combined_cg_mm": components["combined_cg_m"] * 1e3,
            "combined_tau_Nm": components["combined_tau_Nm"],
            "rail_residual_after_kgm": residual_after_now,
        })

        rows.append(row)

    return rows, assignment


# =============================================================================
# 7. PARAMETRIC SWEEPS FOR ENVELOPE PLOTS
# =============================================================================

def static_envelope_sweep():
    """Sweep fill fraction for tank imbalance and slosh envelope plots."""
    phi_vec = np.linspace(0.02, 0.99, 250)
    prop_vec = phi_vec * N_TANKS * RHO_LMP103S_KG_M3 * V_TANK_M3

    results = {}
    for label, eps in TANK_IMBALANCE_CASES.items():
        dcg = []
        tau = []
        for prop in prop_vec:
            # No tugs in this generic tank-only sweep.
            d, total_mass = tank_imbalance_cg_shift(prop, [], np.empty((0, 3)), eps)
            dcg.append(d)
            tau.append(d * MAIN_THRUST_N)
        results[label] = {"phi": phi_vec, "dcg_m": np.array(dcg), "tau_Nm": np.array(tau)}

    slosh = {}
    for label, pmd in PMD_CASES.items():
        dcg = []
        tau = []
        for prop in prop_vec:
            total_mass = M_DRY_KG + prop
            d, _, _ = slosh_static_cg_shift(prop, total_mass, pmd_factor=pmd)
            dcg.append(d)
            tau.append(d * MAIN_THRUST_N)
        slosh[label] = {"phi": phi_vec, "dcg_m": np.array(dcg), "tau_Nm": np.array(tau)}

    return results, slosh


# =============================================================================
# 8. PRINTING FUNCTIONS
# =============================================================================

def print_tug_loading_summary(assignment):
    uniform_xe = 5 * max([t["xe_req_kg"] for t in TARGET_SEQUENCE])
    tailored_xe = np.sum(TUG_XE_LOADS_KG)
    uniform_tug_mass = 5 * (TUG_DRY_MASS_KG + max([t["xe_req_kg"] for t in TARGET_SEQUENCE]))
    tailored_tug_mass = np.sum(TUG_WET_MASSES_KG)

    print("\n" + "=" * 96)
    print("TUG XENON LOADING SUMMARY")
    print("=" * 96)
    print(f"Tug dry mass                      : {TUG_DRY_MASS_KG:.1f} kg")
    print(f"Uniform Xe loading case            : {uniform_xe:.1f} kg Xe total")
    print(f"Tailored Xe loading case           : {tailored_xe:.1f} kg Xe total")
    print(f"Unused Xe reduction                : {uniform_xe - tailored_xe:.1f} kg")
    print(f"Uniform total tug mass             : {uniform_tug_mass:.1f} kg")
    print(f"Tailored total tug mass            : {tailored_tug_mass:.1f} kg")
    print(f"Mothership carried-mass saving     : {uniform_tug_mass - tailored_tug_mass:.1f} kg")
    print("\nTailored tug masses:")
    print(f"{'Tug':>4} {'Xe load [kg]':>14} {'Wet mass [kg]':>15}")
    print("-" * 40)
    for i, (xe, wet) in enumerate(zip(TUG_XE_LOADS_KG, TUG_WET_MASSES_KG), start=1):
        print(f"{i:4d} {xe:14.1f} {wet:15.1f}")

    print("\nRelease assignment:")
    print(f"Strategy: {RELEASE_ASSIGNMENT_STRATEGY}")
    print(f"{'Target':<18} {'Req Xe [kg]':>12} {'Tug':>6} {'Tug Xe [kg]':>14}")
    print("-" * 56)
    for a in assignment:
        print(f"{a['target']:<18} {a['xe_req_kg']:12.1f} {a['tug_id']+1:6d} {a['tug_xe_load_kg']:14.1f}")
    print("=" * 96)


def print_campaign_summary(rows):
    print("\n" + "=" * 150)
    print("MISSION CG AND TORQUE SUMMARY")
    print("=" * 150)
    print(
        f"{'Day':>6} {'Event':<30} {'Type':<8} {'dV':>7} "
        f"{'Prop used':>10} {'Prop rem':>10} {'Mass aft':>10} {'phi':>6} "
        f"{'Rel jump':>10} {'Tug CG':>9} {'Tank CG':>9} {'Slosh CG':>10} {'Comb CG':>9} {'Comb tau':>10}"
    )
    print("-" * 150)
    for r in rows:
        dv = "-" if r["dv_m_s"] is None else f"{r['dv_m_s']:.1f}"
        print(
            f"{r['day']:6.1f} {r['event']:<30.30} {r['type']:<8} {dv:>7} "
            f"{r['prop_used_kg']:10.1f} {r['prop_after_kg']:10.1f} {r['total_mass_after_kg']:10.1f} {r['fill_fraction']:6.3f} "
            f"{r['cg_release_jump_mm']:10.1f} {r['tug_cg_balanced_mm']:9.1f} {r['tank_imb_cg_mm']:9.1f} "
            f"{r['slosh_cg_mm']:10.1f} {r['combined_cg_mm']:9.1f} {r['combined_tau_Nm']:10.2f}"
        )
    print("=" * 150)

    worst_combined = max(rows, key=lambda r: r["combined_tau_Nm"])
    worst_release = max(rows, key=lambda r: r["tau_release_jump_Nm"])

    print("\nWORST-CASE RESULTS")
    print("-" * 72)
    print(f"Worst balanced/combined envelope event : Day {worst_combined['day']:.1f} – {worst_combined['event']}")
    print(f"  Combined CG shift                    : {worst_combined['combined_cg_mm']:.1f} mm")
    print(f"  Main-burn torque envelope            : {worst_combined['combined_tau_Nm']:.2f} N m")
    print(f"Worst immediate release jump event     : Day {worst_release['day']:.1f} – {worst_release['event']}")
    print(f"  Immediate tug-release CG jump        : {worst_release['cg_release_jump_mm']:.1f} mm")
    print(f"  Immediate tug-release torque         : {worst_release['tau_release_jump_Nm']:.2f} N m")
    print("-" * 72)


def print_mode_summary(rows):
    # Use the worst combined event for mode ranking.
    worst = max(rows, key=lambda r: r["combined_tau_Nm"])
    pmd = PMD_CASES[DEFAULT_PMD_CASE]
    mode_rows = mode_criticality(
        worst["fill_fraction"],
        worst["prop_after_kg"],
        worst["total_mass_after_kg"],
        pmd_factor=pmd,
    )

    print("\n" + "=" * 90)
    print("SLOSH MODE CRITICALITY AT WORST COMBINED EVENT")
    print("=" * 90)
    print(f"Event: Day {worst['day']:.1f} – {worst['event']}")
    print(f"Fill fraction: {worst['fill_fraction']:.3f}")
    print(f"PMD case: {DEFAULT_PMD_CASE}")
    print(f"{'Mode':<24} {'m_eff [kg]':>12} {'dCG [mm]':>12} {'Torque [N m]':>14}")
    print("-" * 70)
    for r in mode_rows:
        print(f"{r['mode']:<24} {r['m_eff_kg']:12.1f} {r['dcg_mm']:12.2f} {r['torque_Nm']:14.3f}")
    print("=" * 90)


# =============================================================================
# 9. PLOTTING FUNCTIONS
# =============================================================================

def make_plots(rows, tank_sweep, slosh_sweep):
    plt.rcParams.update({
        "font.size": 9,
        "axes.grid": True,
        "grid.alpha": 0.35,
        "figure.dpi": 130,
    })

    burn_rows = [r for r in rows if r["type"] == "burn"]
    release_rows = [r for r in rows if r["type"] == "release"]

    # ------------------------------------------------------------------
    # Figure 1: mission timeline plots
    # ------------------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(13, 8))
    fig.suptitle("REAVER Mothership – Mission CG and Propellant Timeline", fontweight="bold")

    burn_idx = np.arange(1, len(burn_rows) + 1)
    prop_used = np.array([r["prop_used_kg"] for r in burn_rows])
    axs[0, 0].bar(burn_idx, prop_used)
    axs[0, 0].set_xlabel("Main burn number")
    axs[0, 0].set_ylabel("Chemical propellant used [kg]")
    axs[0, 0].set_title("Propellant Consumed per Main Burn")

    days = np.array([r["day"] for r in rows])
    prop_remaining = np.array([r["prop_after_kg"] for r in rows])
    axs[0, 1].step(days, prop_remaining, where="post")
    axs[0, 1].set_xlabel("Mission day")
    axs[0, 1].set_ylabel("Remaining LMP-103S [kg]")
    axs[0, 1].set_title("Chemical Propellant Remaining")

    phi = np.array([r["fill_fraction"] for r in rows])
    axs[1, 0].plot(days, phi, marker="o")
    axs[1, 0].set_xlabel("Mission day")
    axs[1, 0].set_ylabel("Mean tank fill fraction [-]")
    axs[1, 0].set_title("Tank Fill Fraction Through Campaign")

    combined_tau = np.array([r["combined_tau_Nm"] for r in rows])
    release_tau = np.array([r["tau_release_jump_Nm"] for r in rows])
    axs[1, 1].plot(days, combined_tau, marker="o", label="Balanced combined envelope")
    axs[1, 1].plot(days, release_tau, marker="x", label="Immediate release jump")
    axs[1, 1].set_xlabel("Mission day")
    axs[1, 1].set_ylabel("Torque [N m]")
    axs[1, 1].set_title("Worst-Case Torque Along Campaign")
    axs[1, 1].legend()

    plt.tight_layout()
    plt.savefig("reaver_mission_cg_timeline.png", bbox_inches="tight")

    # ------------------------------------------------------------------
    # Figure 2: tug release and rail rebalancing
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(11, 5))
    labels = [f"{r['target']}\nTug {r['released_tug_id']+1}" for r in release_rows]
    x = np.arange(len(release_rows))
    jump = np.array([r["cg_release_jump_mm"] for r in release_rows])
    balanced = np.array([r["tug_cg_balanced_mm"] for r in release_rows])

    width = 0.38
    ax.bar(x - width/2, jump, width, label="Immediately after release")
    ax.bar(x + width/2, balanced, width, label="After rail rebalancing")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Tug-induced lateral CG shift [mm]")
    ax.set_title("Tug Release CG Jump and Rail Rebalancing")
    ax.legend()
    plt.tight_layout()
    plt.savefig("reaver_tug_release_rebalancing.png", bbox_inches="tight")

    # ------------------------------------------------------------------
    # Figure 3: generic tank/slosh fill-fraction envelopes
    # ------------------------------------------------------------------
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    for label, d in tank_sweep.items():
        axs[0].plot(d["phi"], d["dcg_m"] * 1e3, label=label)
    axs[0].set_xlabel("Fill fraction [-]")
    axs[0].set_ylabel("Tank imbalance CG shift [mm]")
    axs[0].set_title("Chemical Tank Drain-Imbalance Envelope")
    axs[0].legend(fontsize=8)

    for label, d in slosh_sweep.items():
        axs[1].plot(d["phi"], d["dcg_m"] * 1e3, label=label)
    axs[1].set_xlabel("Fill fraction [-]")
    axs[1].set_ylabel("Static slosh CG envelope [mm]")
    axs[1].set_title("Conservative Slosh CG Envelope")
    axs[1].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig("reaver_tank_slosh_envelopes.png", bbox_inches="tight")

    # ------------------------------------------------------------------
    # Figure 4: dynamic slosh at worst combined event
    # ------------------------------------------------------------------
    worst = max(rows, key=lambda r: r["combined_tau_Nm"])
    dyn_cases = {
        "step": dynamic_slosh_response(worst["fill_fraction"], worst["prop_after_kg"],
                                        worst["total_mass_after_kg"],
                                        PMD_CASES[DEFAULT_PMD_CASE], "step"),
        "pulse": dynamic_slosh_response(worst["fill_fraction"], worst["prop_after_kg"],
                                         worst["total_mass_after_kg"],
                                         PMD_CASES[DEFAULT_PMD_CASE], "pulse"),
        "sine": dynamic_slosh_response(worst["fill_fraction"], worst["prop_after_kg"],
                                        worst["total_mass_after_kg"],
                                        PMD_CASES[DEFAULT_PMD_CASE], "sine"),
    }

    fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for name, d in dyn_cases.items():
        axs[0].plot(d["t_s"], d["x_slosh_m"] * 1e3, label=name)
        axs[1].plot(d["t_s"], d["delta_cg_m"] * 1e3, label=name)
        axs[2].plot(d["t_s"], d["torque_Nm"], label=name)
    axs[0].set_ylabel("Slosh mass displacement [mm]")
    axs[1].set_ylabel("Dynamic CG shift [mm]")
    axs[2].set_ylabel("Torque [N m]")
    axs[2].set_xlabel("Time [s]")
    axs[0].set_title("Dynamic First-Mode Slosh Response at Worst Event")
    for ax in axs:
        ax.legend()
    plt.tight_layout()
    plt.savefig("reaver_dynamic_slosh_response.png", bbox_inches="tight")

    print("\nPlots saved:")
    print("  reaver_mission_cg_timeline.png")
    print("  reaver_tug_release_rebalancing.png")
    print("  reaver_tank_slosh_envelopes.png")
    print("  reaver_dynamic_slosh_response.png")


# =============================================================================
# 10. REPORT-READY TEXT OUTPUT
# =============================================================================

def print_report_paragraph(rows):
    worst_combined = max(rows, key=lambda r: r["combined_tau_Nm"])
    worst_release = max(rows, key=lambda r: r["tau_release_jump_Nm"])

    print("\n" + "=" * 96)
    print("REPORT-READY SUMMARY PARAGRAPH")
    print("=" * 96)
    print(f"""
The mothership CG envelope was evaluated by combining chemical propellant
slosh, unequal tank draining, tailored tug loading and sequential tug release.
The four LMP-103S tanks were modelled as spherical PMD tanks placed at the
corners of a 1 m square around the spacecraft centreline. Equal tank draining
therefore cancels first-order lateral CG motion, while paired drain imbalance
was retained as an operational worst-case envelope. Sloshing was modelled with
a conservative spherical-tank liquid-centroid formulation, reduced by a PMD
mobile-mass factor, and with a first-mode spring-mass-damper model for dynamic
excitation.

The tailored tug loading reduces the total carried xenon from {5*99.5:.1f} kg
to {np.sum(TUG_XE_LOADS_KG):.1f} kg, saving approximately
{5*(TUG_DRY_MASS_KG+99.5)-np.sum(TUG_WET_MASSES_KG):.1f} kg of carried tug mass
relative to the uniform worst-case loading. However, the tugs no longer have
equal wet mass, so equal angular spacing on the 1 m radius tug rail would not
cancel the lateral mass moment. The rail is therefore modelled as a CG management
mechanism: before and after each release, the remaining tug angles are optimised
to minimise the lateral mass moment. The largest balanced combined CG envelope
occurs at day {worst_combined['day']:.1f} during '{worst_combined['event']}',
with a conservative CG shift of {worst_combined['combined_cg_mm']:.1f} mm,
corresponding to a main-burn torque of {worst_combined['combined_tau_Nm']:.2f} N m
for the 88 N thruster cluster. The largest transient release jump occurs during
'{worst_release['event']}', producing an immediate tug-induced CG jump of
{worst_release['cg_release_jump_mm']:.1f} mm before rail rebalancing, equivalent
to {worst_release['tau_release_jump_Nm']:.2f} N m. This shows that the dominant
CG driver is the sequential release and rebalancing of unequal-mass tugs, while
chemical slosh remains a secondary but relevant disturbance contribution for the
AOCS torque margin.
""")
    print("=" * 96)


# =============================================================================
# 11. MAIN EXECUTION
# =============================================================================

def main():
    rows, assignment = simulate_campaign()
    tank_sweep, slosh_sweep = static_envelope_sweep()

    print("=" * 96)
    print("REAVER MOTHERSHIP – SLOSH, TUG RELEASE, CG SHIFT AND TORQUE ANALYSIS")
    print("=" * 96)
    print(f"Dry mass                         : {M_DRY_KG:.1f} kg")
    print(f"Initial chemical propellant       : {M_CHEM_INITIAL_KG:.1f} kg")
    print(f"Tank volume per tank              : {V_TANK_M3:.4f} m^3")
    print(f"Tank radius                       : {R_TANK_M:.3f} m")
    print(f"Initial fill fraction             : {fill_fraction_from_prop_mass(M_CHEM_INITIAL_KG):.3f}")
    print(f"Main thrust                       : {MAIN_THRUST_N:.1f} N")
    print(f"Default tank imbalance            : {DEFAULT_TANK_IMBALANCE_FOR_COMBINED*100:.1f}%")
    print(f"Default PMD case                  : {DEFAULT_PMD_CASE}")

    print_tug_loading_summary(assignment)
    print_campaign_summary(rows)
    print_mode_summary(rows)
    print_report_paragraph(rows)
    make_plots(rows, tank_sweep, slosh_sweep)
    plt.show()


if __name__ == "__main__":
    main()
