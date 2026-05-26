import math
import pandas as pd

# ============================================================
# Constants
# ============================================================

G0 = 9.80665          # m/s^2
N_TUGS = 4
TARGET_MASS = 2000   # kg, worst-case defunct satellite mass


# ============================================================
# Helper functions
# ============================================================

def ask_float(prompt):
    return float(input(prompt))


def ask_dv_list(prompt, expected_length):
    """
    Allows user to input either:
    - one value, e.g. 1000
    - several values separated by commas, e.g. 500, 400, 300
    """
    raw = input(prompt).strip()
    values = [float(x.strip()) for x in raw.split(",")]

    if len(values) == 1:
        values = values * expected_length

    if len(values) != expected_length:
        raise ValueError(f"Expected {expected_length} values, but got {len(values)}.")

    return values


def propellant_needed(final_mass, delta_v, isp):
    """
    Rocket equation:
    m0 = mf * exp(dv / (Isp*g0))
    propellant = m0 - mf
    """
    mass_ratio = math.exp(delta_v / (isp * G0))
    initial_mass = final_mass * mass_ratio
    propellant = initial_mass - final_mass
    return propellant, initial_mass, mass_ratio


# ============================================================
# User inputs
# ============================================================

print("\n=== REAVER Mothership + Tug Propellant Sizing ===\n")

mother_dry = ask_float("Mothership dry mass excluding tugs [kg]: ")
tug_dry = ask_float("Dry mass of each tug [kg]: ")

mother_isp = ask_float("Mothership propulsion Isp [s]: ")
tug_isp = ask_float("Tug propulsion Isp [s]: ")

print("\nEnter mothership delta-v values in m/s.")
print("There are 6 legs:")
print("1: Base station -> Debris 1")
print("2: Debris 1 -> Debris 2")
print("3: Debris 2 -> Debris 3")
print("4: Debris 3 -> Debris 4")
print("5: Debris 4 -> Debris 5")
print("6: Debris 5 -> Base station with captured debris")

mother_dvs = ask_dv_list(
    "\nEnter 6 mothership delta-v values separated by commas [m/s]: ",
    6
)

print("\nEnter tug return delta-v values in m/s.")
print("These are for each tug carrying one 2000 kg target back to the base station.")
print("You can enter one value to use the same delta-v for all 4 tugs.")

tug_dvs = ask_dv_list(
    "\nEnter 4 tug return delta-v values separated by commas [m/s]: ",
    4
)


# ============================================================
# 1. Calculate tug propellant and wet masses
# ============================================================

tug_rows = []
tug_wet_masses = []

for i, dv in enumerate(tug_dvs, start=1):
    final_mass = TARGET_MASS + tug_dry

    prop, initial_mass_with_target, mass_ratio = propellant_needed(
        final_mass=final_mass,
        delta_v=dv,
        isp=tug_isp
    )

    tug_wet = tug_dry + prop
    tug_wet_masses.append(tug_wet)

    tug_rows.append({
        "Tug": f"Tug {i}",
        "Delta-v [m/s]": dv,
        "Target mass [kg]": TARGET_MASS,
        "Tug dry mass [kg]": tug_dry,
        "Tug propellant [kg]": prop,
        "Tug wet mass before deployment [kg]": tug_wet,
        "Combined mass with target before return [kg]": initial_mass_with_target,
        "Mass ratio": mass_ratio
    })

tug_table = pd.DataFrame(tug_rows)


# ============================================================
# 2. Calculate mothership propellant backwards
# ============================================================
# We work backwards because the mothership mass changes after
# each tug deployment and after capturing the 5th target.

# Final condition after arriving back at base:
# mothership dry + captured 5th debris, with mothership propellant spent.
mass_after_leg_6 = mother_dry + TARGET_MASS

# Leg 6: Debris 5 -> Base station with captured debris
mr6 = math.exp(mother_dvs[5] / (mother_isp * G0))
mass_before_leg_6 = mass_after_leg_6 * mr6

# Before capturing target 5
mass_after_leg_5 = mass_before_leg_6 - TARGET_MASS

# Leg 5: Debris 4 -> Debris 5
mr5 = math.exp(mother_dvs[4] / (mother_isp * G0))
mass_before_leg_5 = mass_after_leg_5 * mr5

# Go backwards through tug deployments
# Forward: after leg 4, deploy tug 4, then leg 5
# Backwards: before leg 5, add tug 4 to get after leg 4
mass_after_leg_4 = mass_before_leg_5 + tug_wet_masses[3]

mr4 = math.exp(mother_dvs[3] / (mother_isp * G0))
mass_before_leg_4 = mass_after_leg_4 * mr4

mass_after_leg_3 = mass_before_leg_4 + tug_wet_masses[2]

mr3 = math.exp(mother_dvs[2] / (mother_isp * G0))
mass_before_leg_3 = mass_after_leg_3 * mr3

mass_after_leg_2 = mass_before_leg_3 + tug_wet_masses[1]

mr2 = math.exp(mother_dvs[1] / (mother_isp * G0))
mass_before_leg_2 = mass_after_leg_2 * mr2

mass_after_leg_1 = mass_before_leg_2 + tug_wet_masses[0]

mr1 = math.exp(mother_dvs[0] / (mother_isp * G0))
initial_stack_mass = mass_after_leg_1 * mr1

total_tug_wet_mass = sum(tug_wet_masses)
mother_propellant = initial_stack_mass - mother_dry - total_tug_wet_mass


# ============================================================
# 3. Forward stage table for mothership
# ============================================================

stage_rows = []

current_mass = initial_stack_mass

for i, dv in enumerate(mother_dvs, start=1):
    mass_before_burn = current_mass
    mass_ratio = math.exp(dv / (mother_isp * G0))
    mass_after_burn = mass_before_burn / mass_ratio
    prop_burned = mass_before_burn - mass_after_burn

    action = ""
    mass_after_action = mass_after_burn

    if i <= 4:
        tug_mass = tug_wet_masses[i - 1]
        mass_after_action = mass_after_burn - tug_mass
        action = f"Deploy Tug {i} ({tug_mass:.1f} kg removed)"

    elif i == 5:
        mass_after_action = mass_after_burn + TARGET_MASS
        action = f"Capture target 5 ({TARGET_MASS:.1f} kg added)"

    elif i == 6:
        action = "Arrive at base station with target 5"

    stage_rows.append({
        "Leg": i,
        "Manoeuvre": [
            "Base station -> Debris 1",
            "Debris 1 -> Debris 2",
            "Debris 2 -> Debris 3",
            "Debris 3 -> Debris 4",
            "Debris 4 -> Debris 5",
            "Debris 5 -> Base station"
        ][i - 1],
        "Delta-v [m/s]": dv,
        "Mass before burn [kg]": mass_before_burn,
        "Mothership propellant burned [kg]": prop_burned,
        "Mass after burn [kg]": mass_after_burn,
        "Action after burn": action,
        "Mass after action [kg]": mass_after_action
    })

    current_mass = mass_after_action

stage_table = pd.DataFrame(stage_rows)


# ============================================================
# 4. Summary
# ============================================================

total_tug_propellant = tug_table["Tug propellant [kg]"].sum()
total_propellant = mother_propellant + total_tug_propellant

summary = pd.DataFrame([{
    "Mothership dry mass [kg]": mother_dry,
    "Tug dry mass each [kg]": tug_dry,
    "Number of tugs": N_TUGS,
    "Target mass each [kg]": TARGET_MASS,
    "Total tug wet mass [kg]": total_tug_wet_mass,
    "Mothership propellant [kg]": mother_propellant,
    "Total tug propellant [kg]": total_tug_propellant,
    "Total propellant [kg]": total_propellant,
    "Initial stack mass at base [kg]": initial_stack_mass,
    "Final mass at base with target 5 [kg]": current_mass
}])


# ============================================================
# 5. Print and save results
# ============================================================

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 180)

print("\n=== Tug Propellant Table ===")
print(tug_table.round(2))

print("\n=== Mothership Stage Table ===")
print(stage_table.round(2))

print("\n=== Mission Summary ===")
print(summary.round(2))

tug_table.to_csv("tug_propellant_table.csv", index=False)
stage_table.to_csv("mothership_stage_table.csv", index=False)
summary.to_csv("mission_mass_summary.csv", index=False)

print("\nFiles saved:")
print("- tug_propellant_table.csv")
print("- mothership_stage_table.csv")
print("- mission_mass_summary.csv")