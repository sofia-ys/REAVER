from spacecraft import SpaceCraft

'''ISC CALCULATIONS'''
n_targets = 1
dv_list = [0.33, 0.33] # [0.306, 0.306, 0.306, 0.306, 0.306, 0.306]
m_debris_list = [1800]  # [2000, 1800, 1500, 1800, 1000]
Isp = 4200  # electric: 4200, chemical: 270
sc_type = "complex"
n_redundancy = 1
capture_type = "gecko_adhesive"  # choose from: "robotic_manipulator", "gecko_adhesive", "clamp", "payload_adapter", "mev_mechanism"
mission_type = "ISC"  # choose from: "STR", "MTR", "PSW", "M&T", "ISC"

sc = SpaceCraft(n_targets, dv_list, m_debris_list, Isp, sc_type, n_redundancy, capture_type, mission_type)
m_dry = sc.mass()
percentages = sc.mass_breakdown()

print("------------------------------------------")
print(m_dry)
print(percentages)

'''MTR CALCULATIONS'''
# n_targets = 5
# dv_list = [0.306, 0.306, 0.306, 0.306, 0.306, 0.306] # [0.306, 0.306, 0.306, 0.306, 0.306, 0.306]
# m_debris_list = [2000, 1800, 1500, 1800, 1000]
# Isp = 4200  # electric: 4200, chemical: 270
# sc_type = "complex"
# n_redundancy = 1
# capture_type = "robotic_manipulator"  # choose from: "robotic_manipulator", "gecko_adhesive", "clamp", "payload_adapter", "mev_mechanism"
# mission_type = "MTR"  # choose from: "STR", "MTR", "PSW", "M&T", "ISC"

# sc = SpaceCraft(n_targets, dv_list, m_debris_list, Isp, sc_type, n_redundancy, capture_type, mission_type)
# m_dry = sc.mass()
# percentages = sc.mass_breakdown()

# print("------------------------------------------")
# print(m_dry)
# print(percentages)


'''STR CALCULATIONS'''
# n_targets = 1
# dv_list = [0.33, 0.33] # [0.306, 0.306, 0.306, 0.306, 0.306, 0.306]
# m_debris_list = [1800]  # [2000, 1800, 1500, 1800, 1000]
# Isp = 270  # electric: 4200, chemical: 270
# sc_type = "complex"
# n_redundancy = 1
# capture_type = "robotic_manipulator"  # choose from: "robotic_manipulator", "gecko_adhesive", "clamp", "payload_adapter", "mev_mechanism"
# mission_type = "STR"  # choose from: "STR", "MTR", "PSW", "M&T", "ISC"

# sc = SpaceCraft(n_targets, dv_list, m_debris_list, Isp, sc_type, n_redundancy, capture_type, mission_type)
# m_dry = sc.mass()
# percentages = sc.mass_breakdown()

# print("------------------------------------------")
# print(m_dry)
# print(percentages)
