import numpy as np
from scipy.optimize import fsolve
import itertools
import random
import warnings

warnings.filterwarnings('ignore', 'The iteration is not making good progress')

# ==========================================
# EDITABLE VARIABLES (Original Units: km, deg)
# ==========================================
NUM_SIMULATIONS = 100
n_targets = 5

ALT_MIN = 35700  # [km]
ALT_MAX = 36500  # [km]
INC_MIN = 0      # [deg]
INC_MAX = 20     # [deg]

RAAN_START = 240 # [deg] #this is from statistical analysis by Jules
RAAN_END = 90    # [deg]

# SIMILARITY THRESHOLDS: Regenerate if a new orbit is closer than ALL these values to an existing one
MIN_ALT_DIFF = 10.0   # [km]
MIN_INC_DIFF = 0.5    # [deg]
MIN_RAAN_DIFF = 0.5   # [deg]
# ==========================================

# Earth parameters
mu_Earth = 398600.4415  # [km^3 / s^2]
r_Earth = 6378.137      # [km]

def orbit_v(r):
    return np.sqrt(mu_Earth/r)

def deltaV(orbit1, orbit2):
    r_1 = orbit1[0]
    r_2 = orbit2[0]
    i_1 = orbit1[1]
    i_2 = orbit2[1]
    RAAN_1 = orbit1[2]
    RAAN_2 = orbit2[2]

    v1 = orbit_v(r_1)
    v2 = orbit_v(r_2)
    vt1 = np.sqrt(mu_Earth * (2/r_1 - 2/(r_1 + r_2)))
    vt2 = np.sqrt(mu_Earth * (2/r_2 - 2/(r_1 + r_2)))

    cos_theta = np.cos(i_1) * np.cos(i_2) + np.sin(i_1) * np.sin(i_2) * np.cos(RAAN_2 - RAAN_1) 
    cos_theta = np.clip(cos_theta, -1.0, 1.0) 
    theta = np.arccos(cos_theta)

    if theta == 0.0:  
        return abs(v1 - v2)

    def dv_total(t1):
        dV1 = np.sqrt(v1**2 + vt1**2 - 2*v1*vt1*np.cos(t1))
        dV2 = np.sqrt(v2**2 + vt2**2 - 2*v2*vt2*np.cos(theta - t1))
        return dV1 + dV2

    def f(theta1):
        denom1 = np.sqrt(v1**2 + vt1**2 - 2*v1*vt1*np.cos(theta1))
        denom2 = np.sqrt(v2**2 + vt2**2 - 2*v2*vt2*np.cos(theta - theta1))
        term1 = 0 if denom1 == 0 else (v1 * vt1 * np.sin(theta1) / denom1)
        term2 = 0 if denom2 == 0 else (v2 * vt2 * np.sin(theta - theta1) / denom2)
        return term1 - term2

    # The minimum for near-circular orbits is virtually always at 0 or theta
    candidates = [0.0, theta]
    
    # Check the fsolve root in case there is a valid interior minimum
    theta1_fsolve = fsolve(f, theta/2)[0]
    if not np.isnan(theta1_fsolve) and 0 <= theta1_fsolve <= theta:
        candidates.append(theta1_fsolve)

    # Select the physically accurate minimum delta V
    best_theta1 = min(candidates, key=dv_total)
    return dv_total(best_theta1)


'''ORBITAL PARAMETERS'''
# recycling hub orbital parameters
h_0 = 37586  # [km]
i_0 = 7      # [deg]
RAAN_0 = 0   # [deg]
rh_orbit = (h_0 + r_Earth, i_0 * np.pi/180, RAAN_0 * np.pi/180)

debris_list = list(itertools.permutations(range(1, n_targets + 1))) 

# Tracking variables for worst overall route
worst_case_total_dv = -1
worst_case_orbits = []
worst_case_opt_sequence = []
worst_case_deltav_list = []
worst_case_tug_dvs = []

# Tracking variable for the absolute worst tug from ANY generated debris
absolute_worst_tug_dv = -1
absolute_worst_tug_debris_orbit = None

for sim in range(NUM_SIMULATIONS):
    orbits = [rh_orbit] 
    for _ in range(n_targets):
        valid_orbit = False
        while not valid_orbit:
            h = random.uniform(ALT_MIN, ALT_MAX)
            i = random.uniform(INC_MIN, INC_MAX)
            
            if RAAN_START <= RAAN_END:
                raan = random.uniform(RAAN_START, RAAN_END)
            else:
                span = (360 - RAAN_START) + RAAN_END
                raan = (RAAN_START + random.uniform(0, span)) % 360
                
            new_orbit = (h + r_Earth, i * np.pi/180, raan * np.pi/180)
            
            # Check against all existing orbits in this simulation
            too_similar = False
            for existing_orbit in orbits:
                alt_diff = abs(new_orbit[0] - existing_orbit[0])
                inc_diff = abs(new_orbit[1] - existing_orbit[1]) * 180/np.pi
                
                # Handle RAAN wrap-around at 360 degrees
                r1 = new_orbit[2] * 180/np.pi
                r2 = existing_orbit[2] * 180/np.pi
                raan_diff = abs(r1 - r2)
                raan_diff = min(raan_diff, 360 - raan_diff)
                
                if alt_diff < MIN_ALT_DIFF and inc_diff < MIN_INC_DIFF and raan_diff < MIN_RAAN_DIFF:
                    too_similar = True
                    break  # Break out of the checking loop, orbit is rejected
            
            if not too_similar:
                valid_orbit = True
                orbits.append(new_orbit)

    # Evaluate all debris in this simulation against the RH for the worst tug case
    for d_idx in range(1, n_targets + 1):
        tug_dv = float(deltaV(orbits[d_idx], rh_orbit))
        if tug_dv > absolute_worst_tug_dv:
            absolute_worst_tug_dv = tug_dv
            absolute_worst_tug_debris_orbit = orbits[d_idx]

    # Evaluate all route permutations
    orbits_list = []
    for debris_perm in debris_list:
        orbit_perm = [orbits[0]]  
        for idx in debris_perm:
            orbit_perm.append(orbits[idx])  
        orbit_perm.append(orbits[0])  
        orbits_list.append(orbit_perm)

    deltav_list = []
    for path in orbits_list:  
        dv_tot = []
        for idx in range(len(path) - 1):
            dv_tot.append(float(deltaV(path[idx], path[idx+1])))  
        deltav_list.append(dv_tot)  

    tot_dv = [sum(dv_path) for dv_path in deltav_list]

    # Find optimal sequence for this iteration
    optimal_path_idx = tot_dv.index(min(tot_dv))
    min_dv_for_sim = tot_dv[optimal_path_idx]
    opt_sequence = debris_list[optimal_path_idx]
    opt_dvs = deltav_list[optimal_path_idx]

    # Tug return delta V for the debris stops in the optimal route
    tug_dvs = []
    for k in range(n_targets - 1):
        debris_index = opt_sequence[k]
        tug_dv = float(deltaV(orbits[debris_index], rh_orbit))
        tug_dvs.append(tug_dv)

    # Check against worst overall route seen across simulations
    if min_dv_for_sim > worst_case_total_dv:
        worst_case_total_dv = min_dv_for_sim
        worst_case_orbits = orbits.copy()
        worst_case_opt_sequence = opt_sequence
        worst_case_deltav_list = opt_dvs
        worst_case_tug_dvs = tug_dvs

# Build readable route string
visit_labels = ["RH"] + [f"D{d}" for d in worst_case_opt_sequence] + ["RH"]
visit_str = "  ->  ".join(visit_labels)

print("=" * 65)
print(" WORST CASE OVERALL ROUTE FOUND IN SIMULATIONS")
print("=" * 65)
print(f" Route:     {visit_str}")
print(f" Total dv:  {worst_case_total_dv:.4f} km/s")
print()
print(" dv per manoeuvre leg:")
for k, (frm, to, dv) in enumerate(zip(visit_labels[:-1], visit_labels[1:], worst_case_deltav_list), start=1):
    print(f"   Leg {k} ({frm} -> {to}): {dv:.4f} km/s")

print()
print("=" * 65)
print(" TUG RETURN DELTA V TO RH (For this route)")
print("=" * 65)
for k in range(len(worst_case_tug_dvs)):
    debris_label = visit_labels[k + 1]
    print(f"   Tug {k+1} ({debris_label} -> RH): {worst_case_tug_dvs[k]:.4f} km/s")

print()
print("=" * 65)
print(" ORBITAL PARAMETERS FOR THIS ROUTE")
print("=" * 65)
print(f" RH: Alt = {worst_case_orbits[0][0] - r_Earth:.0f} km, Inc = {worst_case_orbits[0][1]*180/np.pi:.2f} deg, RAAN = {worst_case_orbits[0][2]*180/np.pi:.2f} deg")
for k in range(1, n_targets + 1):
    alt = worst_case_orbits[k][0] - r_Earth
    inc = worst_case_orbits[k][1] * 180/np.pi
    raan = worst_case_orbits[k][2] * 180/np.pi
    print(f" D{k}: Alt = {alt:.0f} km, Inc = {inc:.2f} deg, RAAN = {raan:.2f} deg")

print()
print("=" * 65)
print(" ABSOLUTE WORST TUG DELTA V SEEN IN ALL SIMULATIONS")
print("=" * 65)
print(f" Delta V: {absolute_worst_tug_dv:.4f} km/s")
print(f" Debris Orbit: Alt = {absolute_worst_tug_debris_orbit[0] - r_Earth:.0f} km, "
      f"Inc = {absolute_worst_tug_debris_orbit[1]*180/np.pi:.2f} deg, "
      f"RAAN = {absolute_worst_tug_debris_orbit[2]*180/np.pi:.2f} deg")
print("=" * 65)