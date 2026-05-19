from mothership_subsystems import Propulsion

m_dry_ms = 1800
m_dry_t = 250

propulsion = Propulsion(m_dry_ms, m_dry_t)
rcs_masses = propulsion._base_mass_items()
prop_masses = propulsion.propellant_mass()

print("\n================================TUGS================================")
print(f"Propulsion system dry mass: {rcs_masses[0].mass_kg} [kg]")
print(f"Propellant mass: {prop_masses[0].mass_kg} [kg]")
print("\n=============================MOTHERSHIP=============================")
print(f"Propulsion system dry mass: {rcs_masses[1].mass_kg} [kg]")
print(f"Propellant mass: {prop_masses[1].mass_kg} [kg]")