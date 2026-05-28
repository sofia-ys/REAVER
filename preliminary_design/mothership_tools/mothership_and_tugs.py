from mission_spacecraft import SpaceCraft

n_tugs = 5  # number of tugs 
tug = SpaceCraft(sc_type="tug")
tug_wet_mass = tug.m_prop + tug.m_dry
mothership = SpaceCraft(sc_type="ms", m_wet_t=tug_wet_mass)

print("===============================DRY MASS===============================")
print(f"EACH TUG: {tug.m_dry}")
print(f"ALL TUGS ({n_tugs}x): {n_tugs * tug.m_dry}")
print(f"MOTHERSHIP: {mothership.m_dry}")

print("============================PROPELLANT MASS===========================")
print(f"EACH TUG: {tug.m_prop}")
print(f"ALL TUGS ({n_tugs}x): {n_tugs * tug.m_prop}")
print(f"MOTHERSHIP: {mothership.m_prop}")

print("==============================TOTAL MASS==============================")
print(f"TOTAL DRY MASS: {n_tugs * tug.m_dry + mothership.m_dry}")
print(f"LAUNCH MASS (DRY + PROPELLANT): {n_tugs * (tug.m_dry + tug.m_prop) + (mothership.m_dry + mothership.m_prop)}")

print("==========================TUG MASS BREAKDOWN==========================")
for item in tug.mass_breakdown():
    print(f"{item}: {tug.mass_breakdown()[item]} [kg]")

print("=======================MOTHERSHIP MASS BREAKDOWN======================")
for item in mothership.mass_breakdown():
    print(f"{item}: {mothership.mass_breakdown()[item]} [kg]")