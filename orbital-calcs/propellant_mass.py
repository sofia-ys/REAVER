import numpy as np

'''HOW TO USE:
this tool calculates how much propellant mass is needed for a capture sequence before coming back to the RH
1)  adjust the debris mass -- currently it uses one constant debris mass (so if you capture 5 debris they're all 2000kg) around line 19
    m_debris = 2000
2)  adjust the Isp for your type of propulsion around line 20
    Isp = 220
3)  adjust the dry mass of your spacecraft that is travelling  around line 21
    m_dry = 1440
4)  set the number of targets you are capturing IN A ROW (between 1 and 5) around line 47
    n_targets = 1
5)  copy and paste the list of dv manouevres from the dv tool! around line 48
    dv_list = [0.5196472774421069, 0.5196472774421069]
now you have the propellant mass required for your set of dv manouevres! this is iterative so it converges to a propellant mass (:
    '''

# change these parameters as needed
m_debris = 2000  # [kg]
Isp = 220  # [s] LMP-103S [https://ntrs.nasa.gov/api/citations/20140002595/downloads/20140002595.pdf]
m_dry = 1440  # [kg] just the spacecraft WITHOUT any propellant

# rocket mass equation: dv [km/s], Isp [s], m_final (at the end of the manouevre) [kg] 
# when doing multi-target we need to carry the fuel for the next manouevre for us each time so the m_final varies between manouevres
def propellant_m(dv, Isp, m_final):
    return m_final * (np.e**(1000*dv/(Isp * 9.80665)) - 1)

# increase in dry mass from propellant mass 
# RCS mass (Reaction Control System)
# specifically for monopropellant beacuse it's the only non toxic fuel
def rcs_m(m_prop):
    return 0.178 * m_prop + 7.69

# dv_list is the delta-v required for the sequence of manouevres
# for 1 target: dv_list has the manoeuvres RH --> debris --> RH e.g. dv_list = [0.5196472774421069, 0.5196472774421069]
# for multiple targets: dv_list has the manoeuvres RH --> debris 1 --> debris 2 --> ... --> RH e.g. dv_list = [0.6880360287114226, 0.5325970850852726, 0.15780588989264369] (for two debris)
# dv_list is in [km/s]
# m_dry is the mass of the spacecraft without ANY propellant [kg]
# n_targets is how many debris are being captured
def prop_mass_multiTarget_iteration(m_dry, n_targets, dv_list):
    
    m_prop = 0
    for i in range(n_targets + 1):
        m_final = m_dry + (n_targets - i) * m_debris + m_prop + rcs_m(m_prop)
        m_prop += propellant_m(dv_list[-(i+1)], Isp, m_final)

    return m_prop

# here adjust your parameters
n_targets = 2
dv_list = [0.6880360287114226, 0.5325970850852726, 0.15780588989264369]
m_prop = prop_mass_multiTarget_iteration(m_dry, n_targets, dv_list)
print(m_prop)