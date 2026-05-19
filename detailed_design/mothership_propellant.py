import numpy as np

def propellant_m(dv, Isp, m_final):
    return m_final * (np.e**(1000*dv/(Isp * 9.80665)) - 1)

def rcs_mass(m_prop):
    return 0.178 * m_prop + 7.69

def sequence_prop_mass(dv_list, m_ms_dry, tug_m_wet, m_debris, ms_isp):
    m_prop = 0
    for i in range(len(dv_list)):  # number of manouevres
        m_final = m_ms_dry + (i - 1) * tug_m_wet + m_prop
        if i == 0:
            m_final += m_debris[-1] + tug_m_wet  # for the last manouevre, the mothership also carries the last debris
        m_prop += propellant_m(dv_list[-(i+1)], ms_isp, m_final)
    return m_prop

def mothership_m_prop(dv_list, m_ms_dry, tug_m_wet, m_debris, ms_isp):
    m_prop = sequence_prop_mass(dv_list, m_ms_dry, tug_m_wet, m_debris, ms_isp)  # initialising our propellant mass estimates
    m_prop_prev = 0 
    while abs(m_prop - m_prop_prev) > 1:  # convergence condition
        m_prop_prev = m_prop  # setting the previous estimate so we can compare
        m_prop = sequence_prop_mass(dv_list, m_ms_dry + rcs_mass(m_prop), tug_m_wet, m_debris, ms_isp)
    m_rcs = rcs_mass(m_prop)
    return m_prop, m_rcs

def tug_m_prop(tug_m_dry, m_debris, dv_list_tug, tug_isp):
    tugs_m_prop = []
    for i in range(4):  # 4 tugs
        m_final = tug_m_dry + m_debris[i]
        tug_m_prop = propellant_m(dv_list_tug[i], tug_isp, m_final)
        tugs_m_prop.append(tug_m_prop)
    return tugs_m_prop

dv_list_ms = [0.270169420708688, 0.4266240324320108, 0.5332191202826168, 0.18398021120972913, 0.5325970850852726, 0.15780588989264369]  # [km/s] list of dv for manouevres of the mothership: RH -> D1 -> D2 -> D3 -> D4 -> D5 -> RH
dv_list_tug = [0.36, 0.36, 0.36, 0.36]  # [km/s] dv for each tug going from its debris back to the RH
ms_m_dry = 1900  # [kg] dry mass of the mothership 
tug_m_dry = 200  # [kg] wet mass of the tug (it carries its own propellant)
m_debris = [1800, 1200, 2000, 1000, 1500]  # [kg] mass of the collected debris
ms_isp = 250  # [s] isp for the mothership (chemical propellant)
tug_isp = 4200  # [s] isp for the tug (electric propellant)

tugs_m_prop = tug_m_prop(tug_m_dry, m_debris, dv_list_tug, tug_isp)
tug_m_wet = max(tugs_m_prop) + tug_m_dry  # for all tugs same prop, otherwise [tug_m_dry + prop for prop in tugs_m_prop] if each tug will be sized same, but some underfilled 
ms_m_prop, ms_m_rcs = mothership_m_prop(dv_list_ms, ms_m_dry, tug_m_wet, m_debris, ms_isp)

print(f"========================MOTHERSHIP PROPELLANT========================")
print(f"Total propellant mass: {ms_m_prop} [kg]")
print(f"RCS dry mass: {ms_m_rcs} [kg]")