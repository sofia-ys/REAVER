from spacecraft import SpaceCraft

# chemical propulsion per debris delta-v (RH -> debris -> RH)
worst_dv_per_debris = [0.76, 0.76]
avg_dv_per_debris = [0.33, 0.33]
opt_dv_per_debris = [0.025, 0.025]

# chemical propulsion grand tour delta-v (RH -> debris1 -> debris2 -> ... -> debris5 -> RH)
worst_dv_grand_tour = [0.477, 0.477, 0.477, 0.477, 0.477, 0.477]
avg_dv_grand_tour = [0.306, 0.306, 0.306, 0.306, 0.306, 0.306]
opt_dv_grand_tour = [0.123, 0.123, 0.123, 0.123, 0.123, 0.123]

# electric propulsion one way delta-v (debris -> RH)
worst_dv_one_way_elec = [1.22]
avg_dv_one_way_elec = [0.52]
opt_dv_one_way_elec = [0.05]

# electric propulsion per debris delta-v (RH -> debris -> RH)
worst_dv_per_debris = [1.185, 1.185]
avg_dv_per_debris = [0.515, 0.515]
opt_dv_per_debris = [0.04, 0.04]

# isp values
isp_chem = 253
isp_elec_low = 1800
isp_elec_high = 4200


'''CALCULATIONS'''
n_targets = 1
dv_list = worst_dv_per_debris
m_debris_list = [2000]*5
Isp = isp_elec_high
sc_type = "complex"
n_redundancy = 1
capture_type = "robotic_manipulator"  # choose from: "robotic_manipulator", "gecko_adhesive", "clamp", "payload_adapter", "mev_mechanism"
mission_type = "ISC"  # choose from: "STR", "MTR", "PSW", "M&T", "ISC"

sc = SpaceCraft(n_targets, dv_list, m_debris_list, Isp, sc_type, n_redundancy, capture_type, mission_type)
m_dry = sc.mass()
m_prop = sc.propulsion.m_prop
percentages = sc.mass_breakdown()

print("------------------------------------------")
print(f"DRY MASS: {m_dry}")
print(f"PROPELLANT MASS: {m_prop*5}")
# print(percentages)
