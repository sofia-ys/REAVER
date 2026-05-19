from mission_spacecraft import SpaceCraft

tug = SpaceCraft(sc_type="tug")
tug_wet_mass = tug.m_prop + tug.m_dry
mothership = SpaceCraft(sc_type="ms", m_wet_t=tug_wet_mass)

print("===============================DRY MASS===============================")
print(f"TUG: {tug.m_dry}")
print(f"MOTHERSHIP: {mothership.m_dry}")

print("============================PROPELLANT MASS============================")
print(f"TUG: {tug.m_prop * 5}")
print(f"MOTHERSHIP: {mothership.m_prop}")

print("==============================TOTAL MASS==============================")
print(f"LAUNCH MASS: {4*(tug.m_dry + tug.m_prop) + (mothership.m_dry + mothership.m_prop)}")