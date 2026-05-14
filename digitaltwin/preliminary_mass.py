from spacecraft import SpaceCraft

n_targets = 3
dv_list = [0.1819582576744156, 0.5518394227915688, 0.5325970850852726, 0.15780588989264369]
m_debris_list = [1200, 2000, 1800]
Isp = 220
sc_type = "complex"
n_redundancy = 2

sc = SpaceCraft(n_targets, dv_list, m_debris_list, Isp, sc_type, n_redundancy)
m_dry = sc.mass()

print("------------------------------------------")
print(m_dry)


# import numpy as np
# from propellant_mass import prop_mass_multiTarget

# # certain subsystems (or components of them) are dry mass independent, while others are dependent
# # IDEA: iterative loop that converges to a dry mass
# # initial fixed dry mass made up of the DRY MASS INDEPENDENT components --> plug in dry mass --> recalculate dry mass with all DRY MASS DEPENDENT components --> iterate

# # DRY MASS INDEPENDENT PARAMETERS
# m_dry = 1800  # [kg] initial guess base on NASA NM1 
# n_targets = 1  # targets debris consecutively 
# dv_list = [0.66, 0.66]  # [km/s] single retrieval 0.66 per debris, multi retrieval 1.84 per debris
# m_debris_list = [1500]  # [kg] mass of debris collected
# m_aocs = 141.1  # [kg] for single s/c missions

# # DRY MASS DEPENDENT SUBSYSTEMS  (propulsion, )


# for i in range(10):
#     _, _, m_rcs = prop_mass_multiTarget(m_dry, n_targets, dv_list, m_debris_list)
#     m_eps = 0.294 * m_dry  # based on percentages
    


#     m_dry = m_rcs + m_aocs + m_eps

#     print(m_dry)


# # ''' dry mass from payload mass
# #     earth orbiting s/c [ADSEE p. 256] '''
# # # payload mass range: 20-550 kg
# # def m_dry_brown(m_payload):
# #     # a value: range 3-7, avg 4.8
# #     return 4.8 * m_payload 

# # # payload mass range: 50-950 kg
# # def m_dry_zandbergen(m_payload):
# #     # statisitics: R2 = 0.9436, RSE = 14.6%
# #     return 2.058 * m_payload + 342.8

# # ''' mass fractions
# #     ADSEE p. 268'''

# # # baseline report
# # mass_fractions_br = {
# #     "Payload": 0.3,
# #     "Structures": 0.25,
# #     "TCS": 0.03,
# #     "EPS": 0.2,
# #     "AOCS": 0.1,
# #     "Propulsion": 0.08,
# #     "TTC": 0.02,
# #     "CDH": 0.02
# # }

# # # SMAD, SSE (large GEO telecommunications satellites)
# # mass_fractions_smad = {
# #     "Payload": 0.287,
# #     "Structures": 0.209,
# #     "TCS": 0.043,
# #     "EPS": 0.292,
# #     "AOCS": 0.061,
# #     "Propulsion": 0.049,
# #     "TTC": 0.044
# # }

# # # MediaGlobe (large GEO telecommunications satellites)
# # mass_fractions_smad = {
# #     "Payload": 0.284,
# #     "Structures": 0.171,
# #     "TCS": 0.05,
# #     "EPS": 0.295,
# #     "AOCS": 0.069,
# #     "Propulsion": 0.087,
# #     "TTC": 0.031
# # }

# # '''past mission data'''
# # # NASA nm study
# # # AOCS: 141.1 [kg]
# # # capture: 157 [kg]
# # # dry mass (30% contingency): 2352 [kg]
# # # wet mass: 3694 [kg]

# # # FIXED: propellant tanks, ADCS, 
# # # structure around prop tanks, battery, solar array 
# # # capture mechanism varies 

# # # m_prop, _, m_rcs = prop_mass_multiTarget(m_dry, n_targets, dv_list, m_debris_list)