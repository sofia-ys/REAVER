import numpy as np
from scipy.optimize import fsolve
import itertools

'''HOW TO USE:
inputs for each orbit: 
    altitude, inclination, RAAN
1)  copy this format for the debris and paste it around the line 76 etc
    # debris x orbital parameters
    h_x = 35700  # [km]
    i_x = 20  # [deg]
    RAAN_x = 0  # [deg]
    dx_orbit = (h_x + r_Earth, i_x * np.pi/180, RAAN_x * np.pi/180)
2)  fill in the number of targets you wrote orbital parameters for around line 93
    n_targets = x
3)  run the code
4)  yay you now have a maybe wrong value for delta v and the optimal path to get there
    (3, 1, 2) 1.4242006554439008 --> visit debris 3, 1, then 2 for a deltav total of 1.42 km/s
'''

# Earth parameters
mu_Earth = 398600.4415  # [km^2 / s^2]
r_Earth = 6378.137  # [km]

# calculating orbital velocity given an orbital radius
def orbit_v(r):
    return np.sqrt(mu_Earth/r)

# calculates deltaV of a manouevre (given the radius, inclination, and RAAN of two different orbits)
def deltaV(orbit1, orbit2):
    # extracting the orbital parameters from the inputs
    r_1 = orbit1[0]
    r_2 = orbit2[0]
    i_1 = orbit1[1]
    i_2 = orbit2[1]
    RAAN_1 = orbit1[2]
    RAAN_2 = orbit2[2]

    # velocity things
    v1 = orbit_v(r_1)
    v2 = orbit_v(r_2)
    vt1 = np.sqrt(mu_Earth * (2/r_1 - 2/(r_1 + r_2)))
    vt2 = np.sqrt(mu_Earth * (2/r_2 - 2/(r_1 + r_2)))

    # angle stuff
    cos_theta = np.cos(i_1) * np.cos(i_2) + np.sin(i_1) * np.sin(i_2) * np.cos(RAAN_2 - RAAN_1) 
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # jusssst in case we get bad answer
    theta = np.acos(cos_theta)

    # numerically solving for theta1 with scipy
    def f(theta1):
        return (
            v1 * vt1 * np.sin(theta1)
            / np.sqrt(v1**2 + vt1**2 - 2*v1*vt1*np.cos(theta1))
            -
            v2 * vt2 * np.sin(theta - theta1)
            / np.sqrt(v2**2 + vt2**2 - 2*v2*vt2*np.cos(theta - theta1))
        )

    theta1 = fsolve(f, theta/2)[0]

    # delta V of the two parts
    dV1 = np.sqrt(v1**2 + vt1**2 - 2*v1*vt1*np.cos(theta1))
    dV2 = np.sqrt(v2**2 + vt2**2 - 2*v2*vt2*np.cos(theta - theta1))
    return dV1 + dV2

'''ORBITAL PARAMETERS'''
# recycling hub orbital parameters
h_0 = 37586  # [km]
i_0 = 7  # [deg]
RAAN_0 = 0  # [deg]
rh_orbit = (h_0 + r_Earth, i_0 * np.pi/180, RAAN_0 * np.pi/180)

# debris 1 orbital parameters
h_1 = 35700  # [km]
i_1 = 20  # [deg]
RAAN_1 = 0  # [deg]
d1_orbit = (h_1 + r_Earth, i_1 * np.pi/180, RAAN_1 * np.pi/180)

# debris 2 orbital parameters
h_2 = 38000  # [km]
i_2 = 10  # [deg]
RAAN_2 = 0  # [deg]
d2_orbit = (h_2 + r_Earth, i_2 * np.pi/180, RAAN_2 * np.pi/180)

# debris 3 orbital parameters
h_3 = 36300  # [km]
i_3 = 10  # [deg]
RAAN_3 = 10  # [deg]
d3_orbit = (h_3 + r_Earth, i_3 * np.pi/180, RAAN_3 * np.pi/180)

# debris 4 orbital parameters
h_4 = 36300  # [km]
i_4 = 2  # [deg]
RAAN_4 = 10  # [deg]
d4_orbit = (h_4 + r_Earth, i_4 * np.pi/180, RAAN_4 * np.pi/180)

# debris 5 orbital parameters
h_5 = 36300  # [km]
i_5 = 20  # [deg]
RAAN_5 = 10  # [deg]
d5_orbit = (h_5 + r_Earth, i_5 * np.pi/180, RAAN_5 * np.pi/180)

orbits = [rh_orbit, d1_orbit, d2_orbit, d3_orbit, d4_orbit, d5_orbit]  # list of tuples we can iterate over

n_targets = 5  # CHANGE THIS IF ADDING MORE DEBRIS
debris_list = list(itertools.permutations(range(1, n_targets + 1)))  # all permutations of visiting debris 

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
    dv_tot = []
    for i in range(len(path) - 1):
        dv_tot.append(float(deltaV(path[i], path[i+1])))  # dv for each manuevre taking this certain path

    deltav_list.append(dv_tot)  # now we add it to the list so we can see what total dv each path requires

tot_dv = []
for dv_path in deltav_list:
    tot_dv.append(sum(dv_path))

# RUN HERE
optimal_path = tot_dv.index(min(tot_dv))  # order of the debris visiting that is best
print(debris_list[optimal_path], tot_dv[optimal_path])  # printing that and also how much total delta v it is to visit all debris and go back and forth from rh
print(deltav_list[optimal_path])  # dv of each manuevre taking the optimal path