from mission_spacecraft import SpaceCraft

tug = SpaceCraft(sc_type="tug")
tug_wet_mass = tug.m_prop + tug.m_dry
mothership = SpaceCraft(sc_type="ms", m_wet_t=tug_wet_mass)

print("===============================DRY MASS===============================")
print(f"TUG: {tug.m_dry * 1.2}")
print(f"MOTHERSHIP: {mothership.m_dry * 1.2}")

print("============================PROPELLANT MASS============================")
print(f"TUG: {tug.m_prop}")
print(f"MOTHERSHIP: {mothership.m_prop}")