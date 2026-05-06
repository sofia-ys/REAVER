import numpy as np

'''HOW TO USE:
inputs for each orbit: 
    altitude, inclination, RAAN
1)  create a variable with the orbit input stored
    orbit1 = ()

1)  create an instance of the class:
    deltaV = DeltaV(r_initial, r_final, i_initial, i_final)
2)  run the combined deltav manouevre function (does altitude and inclinatino change)
    delta_v = deltaV.combined()
3)  yay you now have a maybe wrong value for delta v 
'''

# class DeltaV:
#     def __init__(self, r_initial, r_final, i_initial, i_final):
#         # Earth parameters
#         self.mu_Earth = 398600.4415  # [km^2 / s^2]
#         self.r_Earth = 6378.137  # [km]
#         self.v_initial = np.sqrt(self.mu_Earth/r_initial)
#         self.v_final = np.sqrt(self.mu_Earth/r_final)
#         self.delta_i = i_final - i_initial
    
#     # inclination and altitude change
#     def combined(self):
#         deltaV = self.v_initial**2 + self.v_final**2 - 2 * self.v_initial * self.v_final * np.cos(self.delta_i)
#         return deltaV
    
#     def phase(self):
#         return "hello"
    
# # recycling hub parameters
# r_rh = 42164  # GEO radius [km]
# i_rh = 7 * np.pi/180  # inclination [rad]

mu_Earth = 398600.4415  # [km^2 / s^2]
r_Earth = 6378.137  # [km]

def orbit_v(r):
    return np.sqrt(mu_Earth/r)

from scipy.optimize import fsolve

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
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # jusssst in case we get bad answer
    theta = np.acos(cos_theta)

    # numerically solving for theta1 with scipy thing
    def f(theta1):
        return (
            v1 * vt1 * np.sin(theta1)
            / np.sqrt(v1**2 + vt1**2 - 2*v1*vt1*np.cos(theta1))
            -
            v2 * vt2 * np.sin(theta - theta1)
            / np.sqrt(v2**2 + vt2**2 - 2*v2*vt2*np.cos(theta - theta1))
        )

    theta1 = fsolve(f, theta/2)[0]

    dV1 = np.sqrt(v1**2 + vt1**2 - 2*v1*vt1*np.cos(theta1))
    dV2 = np.sqrt(v2**2 + vt2**2 - 2*v2*vt2*np.cos(theta - theta1))
    return dV1 + dV2

import itertools

'''delta_v optimising --> assuming delta-v function'''

# recycling hub orbital parameters
h_0 = 37586  # [km]
i_0 = 7  # [deg]
RAAN_0 = 0  # [deg]
rh_orbit = (h_0 + r_Earth, i_0 * np.pi/180, RAAN_0 * np.pi/180)

# debris 1 orbital parameters
h_1 = 3570  # [km]
i_1 = 20  # [deg]
RAAN_1 = 0  # [deg]
d1_orbit = (h_1 + r_Earth, i_1 * np.pi/180, RAAN_1 * np.pi/180)

# debris 2 orbital parameters
h_2 = 3800  # [km]
i_2 = 10  # [deg]
RAAN_2 = 0  # [deg]
d2_orbit = (h_2 + r_Earth, i_2 * np.pi/180, RAAN_2 * np.pi/180)

# debris 3 orbital parameters
h_3 = 3630  # [km]
i_3 = 10  # [deg]
RAAN_3 = 10  # [deg]
d3_orbit = (h_3 + r_Earth, i_3 * np.pi/180, RAAN_3 * np.pi/180)

orbits = [rh_orbit, d1_orbit, d2_orbit, d3_orbit]

n_targets = 3
debris_list = list(itertools.permutations(range(1, n_targets + 1)))

# all the orbital data in all the possible perms
orbits_list = []
for debris_perm in debris_list:
    orbit_perm = [orbits[0]]  # start is always RH
    for i in debris_perm:
        orbit_perm.append(orbits[i])  # switching up the different debris paths
    orbit_perm.append(orbits[0])  # end is always RH
    orbits_list.append(orbit_perm)

deltav_list = []
for path in orbits_list:  # for each possible path in the orbits list: rh -- debrisx -- debrisx -- ... -- rh
    dv_tot = 0
    for i in range(len(path) - 1):
        dv_tot += deltaV(path[i], path[i+1])

    deltav_list.append(dv_tot)

optimal_path = deltav_list.index(min(deltav_list))
print(debris_list[optimal_path], deltav_list[optimal_path])