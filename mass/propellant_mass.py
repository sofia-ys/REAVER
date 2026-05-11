import numpy as np

'''HOW TO USE:
this tool calculates how much propellant mass is needed for a capture sequence before coming back to the RH
1)  adjust the debris mass (line 19) -- this is the masses of the debris IN THE ORDER YOU'RE CAPTURING THEM (same as the delta-v manouevre order)
    m_debris = [2000, 1200]
2)  adjust the Isp for your type of propulsion (line 20)
    Isp = 220
3)  adjust the dry mass of your spacecraft that is travelling (line 21))
    m_dry = 1440
4)  set the number of targets you are capturing IN A ROW -- between 1 and 5 (line 52)
    n_targets = 1
5)  copy and paste the list of dv manouevres from the dv tool! (line 53)
    dv_list = [0.5196472774421069, 0.5196472774421069]
now you have the propellant mass required for your set of dv manouevres! this is iterative so it converges to a propellant mass (:
    '''

# change these parameters as needed
m_debris_list = [2000, 1800]  # [kg]
Isp = 220  # [s] LMP-103S [https://ntrs.nasa.gov/api/citations/20140002595/downloads/20140002595.pdf]
m_dry = 1440  # [kg] just the spacecraft WITHOUT any propellant (or propellant tanks etc since that's estimated)

# rocket mass equation: dv [km/s], Isp [s], m_final (at the end of the manouevre) [kg] 
# when doing multi-target we need to carry the fuel for the next manouevre for us each time so the m_final varies between manouevres
def propellant_m(dv, Isp, m_final):
    return m_final * (np.e**(1000*dv/(Isp * 9.80665)) - 1)

# increase in dry mass from propellant mass -> RCS mass (Reaction Control System)
# specifically for monopropellant beacuse it's the only non toxic fuel ADSEE p. 166
# RSE = 8.1%, propellant mass in range 30-300kg (extrapolating for this project)
def rcs_m(m_prop):
    return 0.178 * m_prop + 7.69

# dv_list is the delta-v required for the sequence of manouevres
# for 1 target: dv_list has the manoeuvres RH --> debris --> RH e.g. dv_list = [0.5196472774421069, 0.5196472774421069]
# for multiple targets: dv_list has the manoeuvres RH --> debris 1 --> debris 2 --> ... --> RH e.g. dv_list = [0.6880360287114226, 0.5325970850852726, 0.15780588989264369] (for two debris)
# dv_list is in [km/s]
# m_dry is the mass of the spacecraft without ANY propellant [kg]
# n_targets is how many debris are being captured
def prop_mass_multiTarget(m_dry, n_targets, dv_list, m_debris_list):
    
    m_prop = 0
    m_prop_list = []
    for i in range(n_targets + 1):
        m_final = m_dry + sum(m_debris_list[:(n_targets - i)]) + m_prop 
        m_prop += propellant_m(dv_list[-(i+1)], Isp, m_final)
        m_prop_list.append(m_prop)

    return m_prop, m_prop_list[::-1]

# here adjust your parameters
n_targets = 2
dv_list = [0.6880360287114226, 0.5325970850852726, 0.15780588989264369]

# iterating to include propellant tanks dry mass as a part of the dry mass
m_prop, _ = prop_mass_multiTarget(m_dry, n_targets, dv_list, m_debris_list)  # initialising our propellant mass estimates
m_prop_prev = 0 
while abs(m_prop - m_prop_prev) > 1:  # convergence condition
    m_prop_prev = m_prop  # setting the previous estimate so we can compare
    m_prop, m_prop_list = prop_mass_multiTarget(m_dry + rcs_m(m_prop), n_targets, dv_list, m_debris_list)

# just for fancy printing
manouevre_list = ""
for j in range(1, n_targets + 1):
    manouevre_list += "DEBRIS-" + str(j) + " → "

print(f"========================{n_targets}-DEBRIS TRAJECTORY========================")
print(f"Total propellant mass: {m_prop} [kg]")
print(f"RCS dry mass: {rcs_m(m_prop)} [kg]")
print(f"Manouevre trajectory: RH → {manouevre_list}RH ")
print(f"Propellant mass at the start of each manouevre: {m_prop_list} [kg]")
