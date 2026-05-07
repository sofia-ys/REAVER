import numpy as np
from scipy.optimize import fsolve

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
h_1 = 36102  # [km]
i_1 = 2.57  # [deg]
RAAN_1 = 359.47  # [deg]
d1_orbit = (h_1 + r_Earth, i_1 * np.pi/180, RAAN_1 * np.pi/180)

# debris 2 orbital parameters
h_2 = 36015  # [km]
i_2 = 0.48  # [deg]
RAAN_2 = 341.91  # [deg]
d2_orbit = (h_2 + r_Earth, i_2 * np.pi/180, RAAN_2 * np.pi/180)

# debris 3 orbital parameters
h_3 = 36076  # [km]
i_3 = 19.69  # [deg]
RAAN_3 = 7.22  # [deg]
d3_orbit = (h_3 + r_Earth, i_3 * np.pi/180, RAAN_3 * np.pi/180)

# debris 4 orbital parameters
h_4 = 35888  # [km]
i_4 = 18.03  # [deg]
RAAN_4 = 15.44  # [deg]
d4_orbit = (h_4 + r_Earth, i_4 * np.pi/180, RAAN_4 * np.pi/180)

# debris 5 orbital parameters
h_5 = 35775  # [km]
i_5 = 19.42  # [deg]
RAAN_5 = 330.11  # [deg]
d5_orbit = (h_5 + r_Earth, i_5 * np.pi/180, RAAN_5 * np.pi/180)

orbits = [rh_orbit, d1_orbit, d2_orbit, d3_orbit, d4_orbit, d5_orbit]  # list of tuples we can iterate over

# how many debris do we visit???
n_targets = 5

orbits_list = [orbits[0]]
for i in range(1, n_targets + 1):
    orbits_list.append(orbits[i])  # go to debris
    orbits_list.append(orbits[0])  # return to RH

dv_tot = []
for i in range(len(orbits_list) - 1):
    dv_tot.append(float(deltaV(orbits_list[i], orbits_list[i+1])))  # dv for each manuevre taking this certain path

# RUN HERE
print("------------------RH -> 1 DEBRIS -> RH -> 1 DEBRIS -> RH-----------------")
print("total dv: ", sum(dv_tot))  # printing that and also how much total delta v it is to visit all debris and go back and forth from rh
print(dv_tot)  # dv of each manuevre taking the optimal path


'''NUMBER AND MASS OF DRONES'''

# PROPELLANT TYPES
Isp = 220 # [s] LMP-103S [https://ntrs.nasa.gov/api/citations/20140002595/downloads/20140002595.pdf]
# Isp = 261 # [s] ASCENT [https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=5280&context=smallsat]


def propellant_m(dv, Isp, m_wet):
    return m_wet * (1 - np.e**(-dv/(Isp * 9.80665)))

# for the way there we need to transport just drones 
# iteration loop to find prop mass?

# "while a 12-thruster configuration offers an optimal balance between control authority and thrust efficiency."
# [https://arxiv.org/pdf/2601.11802v1]
# because travelling so much distance, HUGE swarm doesn't make sense (since s/c won't be able to carry all propellant)
n_drones = 12
m_drone = 120  # preliminary estimate 
m_dry = m_drone * n_drones
m_debris = 2000
dv = sum(dv_tot)/10 * 1000  # average dv for each way
print("-------------------------")

# iterate to find propellant mass -- rocket equation
def prop_mass_iteration(m_dry):
    m_wet = m_dry + m_debris
    for j in range(10):
        m_wet = m_dry + m_debris + propellant_m(dv, Isp, m_wet)
    m_return_prop = m_wet - (m_dry + m_debris)

    m_wet = m_dry + m_return_prop
    for i in range(10):
        m_wet = m_dry + m_return_prop + propellant_m(dv, Isp, m_wet)
    m_there_prop = m_wet - (m_dry + m_return_prop)

    m_prop = m_there_prop + m_return_prop
    return m_prop

# initial propellant mass estimation using estimated dry mass 
m_prop = prop_mass_iteration(m_dry)
print(m_prop)

# # monopropellant thruster dry mass estimation
# m_dry_prop_sys = 0.178 * m_prop/12 + 7.69
# print(m_dry_prop_sys)

# for TWO targets 
# dv_two_target = [0.6880360287114226, 0.5325970850852726, 0.15780588989264369]

# # iterating for multiple targets 
# def prop_mass_multiTarget_iteration(m_dry):
    
#     # return with both debris 
#     dv = dv_two_target[-1]
#     m_wet = m_dry + 2 * m_debris
#     for j in range(10):
#         m_wet = m_dry + 2 * m_debris + propellant_m(dv, Isp, m_wet)
#     m_return_prop = m_wet - (m_dry + 2 * m_debris)
#     print(m_return_prop)

#     dv = dv_two_target[-2]
#     m_wet = m_dry + m_debris + m_return_prop
#     for k in range(10):
#             m_wet = m_dry + m_return_prop + m_debris + propellant_m(dv, Isp, m_wet)
#     m_middle_prop = m_wet - (m_dry + m_debris + m_return_prop)
#     print(m_middle_prop)

#     dv = dv_two_target[-3]
#     m_wet = m_dry + m_return_prop + m_middle_prop
#     for i in range(10):
#         m_wet = m_dry + m_return_prop + m_middle_prop + propellant_m(dv, Isp, m_wet)
#     m_there_prop = m_wet - (m_dry + m_return_prop + m_middle_prop)
#     print(m_there_prop)

#     m_prop = m_there_prop + m_return_prop + m_middle_prop
#     return m_prop

# m_prop = prop_mass_iteration(m_dry)
# print(m_prop)




def prop_mass_multiTarget_iteration(m_dry, n_targets, dv_list):
    
    m_prop = 0
    for i in range(n_targets):
        dv = 1000 * dv_list[-(i+1)]  # we start with the last manouvre
        m_wet = m_dry + (n_targets - i) * m_debris  
        for j in range(10):
            m_wet = m_dry + (n_targets - i) * m_debris + propellant_m(dv, Isp, m_wet)
        m_prop += m_wet - (m_dry + (n_targets - i) * m_debris)

    return m_prop


# m_prop = prop_mass_multiTarget_iteration(m_dry, 2, [0.6880360287114226, 0.5325970850852726, 0.15780588989264369])

m_prop = prop_mass_multiTarget_iteration(m_dry, 1, [0.5196472774421069, 0.5196472774421069])
print(m_prop)