import math
import itertools
import random
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Constants
# ============================================================

MU_EARTH = 398600.4418   # km^3/s^2
G0 = 9.80665             # m/s^2
TARGET_MASS = 2000       # kg, worst-case defunct satellite
N_TUGS = 4
BASE_A = 42164.0         # km, approximate GEO radius
BASE_INC = 7.0           # deg, base station inclination
BASE_RAAN = 0.0          # deg, assumed base station RAAN


# ============================================================
# User input helpers
# ============================================================

def ask_float(prompt, default=None):
    if default is None:
        return float(input(prompt))

    raw = input(f"{prompt} [{default}]: ").strip()
    if raw == "":
        return float(default)
    return float(raw)


def ask_int(prompt, default=None):
    if default is None:
        return int(input(prompt))

    raw = input(f"{prompt} [{default}]: ").strip()
    if raw == "":
        return int(default)
    return int(raw)


# ============================================================
# Orbital mechanics helper functions
# ============================================================

def deg2rad(angle_deg):
    return angle_deg * math.pi / 180.0


def angular_difference_deg(a, b):
    """
    Returns smallest angular difference between two angles in degrees.
    """
    return abs((a - b + 180.0) % 360.0 - 180.0)


def plane_angle_deg(i1, raan1, i2, raan2):
    """
    Computes true angle between two orbital planes.
    """
    i1 = deg2rad(i1)
    i2 = deg2rad(i2)
    delta_raan = deg2rad(angular_difference_deg(raan1, raan2))

    cos_angle = (
        math.cos(i1) * math.cos(i2)
        + math.sin(i1) * math.sin(i2) * math.cos(delta_raan)
    )

    cos_angle = max(-1.0, min(1.0, cos_angle))
    return math.degrees(math.acos(cos_angle))


def hohmann_delta_v(r1, r2):
    """
    Approximate Hohmann transfer delta-v between circular orbits.
    Inputs in km.
    Output in m/s.
    """
    if abs(r2 - r1) < 1e-6:
        return 0.0

    v1 = math.sqrt(MU_EARTH / r1)
    v2 = math.sqrt(MU_EARTH / r2)

    dv1 = v1 * (math.sqrt(2 * r2 / (r1 + r2)) - 1)
    dv2 = v2 * (1 - math.sqrt(2 * r1 / (r1 + r2)))

    return (abs(dv1) + abs(dv2)) * 1000.0


def plane_change_delta_v(a_avg, plane_angle):
    """
    Approximate pure plane-change delta-v.
    a_avg in km.
    plane_angle in deg.
    Output in m/s.
    """
    v = math.sqrt(MU_EARTH / a_avg)  # km/s
    angle_rad = deg2rad(plane_angle)
    dv = 2 * v * math.sin(angle_rad / 2)
    return dv * 1000.0


def estimate_transfer_dv(obj1, obj2, margin_percent=20.0, extra_dv=0.0):
    """
    Simple approximate transfer delta-v between two GEO-like objects.

    It combines:
    - Hohmann altitude transfer
    - orbital plane change
    - margin
    - fixed extra dv

    You can later replace this function with your accurate optimiser.
    """

    a1 = obj1["a"]
    i1 = obj1["inc"]
    raan1 = obj1["raan"]

    a2 = obj2["a"]
    i2 = obj2["inc"]
    raan2 = obj2["raan"]

    plane_angle = plane_angle_deg(i1, raan1, i2, raan2)

    dv_alt = hohmann_delta_v(a1, a2)
    dv_plane = plane_change_delta_v((a1 + a2) / 2, plane_angle)

    # Simple combination
    dv_raw = math.sqrt(dv_alt**2 + dv_plane**2)

    dv_with_margin = dv_raw * (1 + margin_percent / 100.0) + extra_dv

    return dv_with_margin


# ============================================================
# Rocket equation and mass model
# ============================================================

def propellant_needed(final_mass, delta_v, isp):
    """
    final_mass: mass after burn [kg]
    delta_v: m/s
    isp: s
    """
    mass_ratio = math.exp(delta_v / (isp * G0))
    initial_mass = final_mass * mass_ratio
    propellant = initial_mass - final_mass
    return propellant, initial_mass


def calculate_tug_wet_mass(tug_dry, tug_dv, tug_isp):
    """
    Tug returns one 2000 kg target to base station.
    """
    final_mass = TARGET_MASS + tug_dry
    propellant, initial_mass_with_target = propellant_needed(
        final_mass=final_mass,
        delta_v=tug_dv,
        isp=tug_isp
    )

    tug_wet = tug_dry + propellant

    return tug_wet, propellant, initial_mass_with_target


def calculate_campaign_mass(
    ordered_targets,
    mother_dry,
    tug_dry,
    mother_isp,
    tug_isp,
    margin_percent,
    mother_extra_dv,
    tug_extra_dv,
):
    """
    ordered_targets: list of 5 target dictionaries.

    Mission:
    base -> target 1 -> target 2 -> target 3 -> target 4 -> target 5 -> base

    At targets 1-4:
    - mothership deploys one tug

    At target 5:
    - mothership captures the target itself

    Returns:
    - total initial stack mass
    - mothership propellant
    - tug propellant total
    - detailed stage data
    """

    base = {
        "name": "BASE",
        "norad": "BASE",
        "a": BASE_A,
        "inc": BASE_INC,
        "raan": BASE_RAAN,
    }

    # --------------------------------------------------------
    # Tug return delta-v and wet mass
    # --------------------------------------------------------

    tug_wet_masses = []
    tug_propellants = []
    tug_dvs = []

    for target in ordered_targets[:4]:
        dv_tug = estimate_transfer_dv(
            target,
            base,
            margin_percent=margin_percent,
            extra_dv=tug_extra_dv
        )

        tug_wet, tug_prop, _ = calculate_tug_wet_mass(
            tug_dry=tug_dry,
            tug_dv=dv_tug,
            tug_isp=tug_isp
        )

        tug_dvs.append(dv_tug)
        tug_wet_masses.append(tug_wet)
        tug_propellants.append(tug_prop)

    # --------------------------------------------------------
    # Mothership leg delta-v
    # --------------------------------------------------------

    route = [base] + ordered_targets + [base]

    mother_dvs = []

    for i in range(len(route) - 1):
        dv = estimate_transfer_dv(
            route[i],
            route[i + 1],
            margin_percent=margin_percent,
            extra_dv=mother_extra_dv
        )
        mother_dvs.append(dv)

    # --------------------------------------------------------
    # Mothership mass calculated backwards
    # --------------------------------------------------------

    # Final mass after arriving at base with target 5
    final_mass = mother_dry + TARGET_MASS

    # Work backwards through leg 6
    prop, mass_before_leg_6 = propellant_needed(
        final_mass=final_mass,
        delta_v=mother_dvs[5],
        isp=mother_isp
    )

    # Before capturing target 5
    mass_after_leg_5 = mass_before_leg_6 - TARGET_MASS

    # Work backwards through leg 5
    prop, mass_before_leg_5 = propellant_needed(
        final_mass=mass_after_leg_5,
        delta_v=mother_dvs[4],
        isp=mother_isp
    )

    # Add tug 4 before leg 4
    mass_after_leg_4 = mass_before_leg_5 + tug_wet_masses[3]

    prop, mass_before_leg_4 = propellant_needed(
        final_mass=mass_after_leg_4,
        delta_v=mother_dvs[3],
        isp=mother_isp
    )

    # Add tug 3 before leg 3
    mass_after_leg_3 = mass_before_leg_4 + tug_wet_masses[2]

    prop, mass_before_leg_3 = propellant_needed(
        final_mass=mass_after_leg_3,
        delta_v=mother_dvs[2],
        isp=mother_isp
    )

    # Add tug 2 before leg 2
    mass_after_leg_2 = mass_before_leg_3 + tug_wet_masses[1]

    prop, mass_before_leg_2 = propellant_needed(
        final_mass=mass_after_leg_2,
        delta_v=mother_dvs[1],
        isp=mother_isp
    )

    # Add tug 1 before leg 1
    mass_after_leg_1 = mass_before_leg_2 + tug_wet_masses[0]

    prop, initial_stack_mass = propellant_needed(
        final_mass=mass_after_leg_1,
        delta_v=mother_dvs[0],
        isp=mother_isp
    )

    total_tug_wet_mass = sum(tug_wet_masses)
    mother_propellant = initial_stack_mass - mother_dry - total_tug_wet_mass

    total_tug_propellant = sum(tug_propellants)
    total_propellant = mother_propellant + total_tug_propellant

    return {
        "initial_stack_mass": initial_stack_mass,
        "mother_propellant": mother_propellant,
        "total_tug_propellant": total_tug_propellant,
        "total_propellant": total_propellant,
        "mother_dvs": mother_dvs,
        "tug_dvs": tug_dvs,
        "tug_wet_masses": tug_wet_masses,
        "targets": ordered_targets,
    }


# ============================================================
# Data loading
# ============================================================

def load_targets(file_path):
    df = pd.read_csv(file_path)
    df.columns = df.columns.str.strip()

    required = [
        "OBJECT_NAME",
        "NORAD_CAT_ID",
        "SEMIMAJOR_AXIS",
        "INCLINATION",
        "RA_OF_ASC_NODE",
    ]

    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    for col in ["SEMIMAJOR_AXIS", "INCLINATION", "RA_OF_ASC_NODE"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["SEMIMAJOR_AXIS", "INCLINATION", "RA_OF_ASC_NODE"])

    targets = []

    for _, row in df.iterrows():
        targets.append({
            "name": str(row["OBJECT_NAME"]),
            "norad": str(row["NORAD_CAT_ID"]),
            "a": float(row["SEMIMAJOR_AXIS"]),
            "inc": float(row["INCLINATION"]),
            "raan": float(row["RA_OF_ASC_NODE"]),
        })

    return targets


# ============================================================
# Campaign search
# ============================================================

def evaluate_best_order_for_five_targets(
    five_targets,
    mother_dry,
    tug_dry,
    mother_isp,
    tug_isp,
    margin_percent,
    mother_extra_dv,
    tug_extra_dv,
):
    """
    For a selected group of 5 targets, test all 120 possible orders.
    Return the best and worst order by initial stack mass.
    """

    best_result = None
    worst_result = None

    for order in itertools.permutations(five_targets, 5):
        result = calculate_campaign_mass(
            ordered_targets=list(order),
            mother_dry=mother_dry,
            tug_dry=tug_dry,
            mother_isp=mother_isp,
            tug_isp=tug_isp,
            margin_percent=margin_percent,
            mother_extra_dv=mother_extra_dv,
            tug_extra_dv=tug_extra_dv,
        )

        if best_result is None or result["initial_stack_mass"] < best_result["initial_stack_mass"]:
            best_result = result

        if worst_result is None or result["initial_stack_mass"] > worst_result["initial_stack_mass"]:
            worst_result = result

    return best_result, worst_result


def campaign_to_row(result, label):
    names = [t["name"] for t in result["targets"]]
    norads = [t["norad"] for t in result["targets"]]

    return {
        "case": label,
        "initial_stack_mass_kg": result["initial_stack_mass"],
        "mother_propellant_kg": result["mother_propellant"],
        "total_tug_propellant_kg": result["total_tug_propellant"],
        "total_propellant_kg": result["total_propellant"],
        "mother_total_dv_m_s": sum(result["mother_dvs"]),
        "average_tug_dv_m_s": np.mean(result["tug_dvs"]),
        "target_names": " | ".join(names),
        "norad_ids": " | ".join(norads),
        "mother_dvs_m_s": " | ".join(f"{x:.1f}" for x in result["mother_dvs"]),
        "tug_dvs_m_s": " | ".join(f"{x:.1f}" for x in result["tug_dvs"]),
    }


# ============================================================
# Main script
# ============================================================

def main():
    print("\n=== REAVER Campaign Mass Optimiser ===\n")

    file_path = input("Enter candidate CSV file name [reaver_defunct_candidates_gcat.csv]: ").strip()
    if file_path == "":
        file_path = "defunct_satellites\\reaver_defunct_candidates_orbital_processed.csv"

    mother_dry = ask_float("Mothership dry mass excluding tugs [kg]", 1000)
    tug_dry = ask_float("Dry mass of each tug [kg]", 250)

    mother_isp = ask_float("Mothership Isp [s]", 350)
    tug_isp = ask_float("Tug Isp [s]", 1500)

    margin_percent = ask_float("Delta-v margin [%]", 20)
    mother_extra_dv = ask_float("Extra mothership delta-v per leg [m/s]", 0)
    tug_extra_dv = ask_float("Extra tug delta-v per return [m/s]", 100)

    number_of_samples = ask_int("Number of random 5-target campaigns to test", 2000)

    random_seed = ask_int("Random seed", 42)
    random.seed(random_seed)

    output_dir = Path("mothership\\campaign_mass_results")
    output_dir.mkdir(exist_ok=True)

    targets = load_targets(file_path)

    print(f"\nLoaded {len(targets)} targets.")

    if len(targets) < 5:
        raise ValueError("Need at least 5 target objects.")

    results_rows = []

    global_best = None
    global_worst = None

    # --------------------------------------------------------
    # Random campaign sampling
    # --------------------------------------------------------

    for sample_idx in range(number_of_samples):
        five_targets = random.sample(targets, 5)

        best_order, worst_order = evaluate_best_order_for_five_targets(
            five_targets=five_targets,
            mother_dry=mother_dry,
            tug_dry=tug_dry,
            mother_isp=mother_isp,
            tug_isp=tug_isp,
            margin_percent=margin_percent,
            mother_extra_dv=mother_extra_dv,
            tug_extra_dv=tug_extra_dv,
        )

        results_rows.append(campaign_to_row(best_order, f"sample_{sample_idx}_best_order"))

        if global_best is None or best_order["initial_stack_mass"] < global_best["initial_stack_mass"]:
            global_best = best_order

        if global_worst is None or worst_order["initial_stack_mass"] > global_worst["initial_stack_mass"]:
            global_worst = worst_order

    results_df = pd.DataFrame(results_rows)

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    average_mass = results_df["initial_stack_mass_kg"].mean()
    median_mass = results_df["initial_stack_mass_kg"].median()

    summary_rows = [
        campaign_to_row(global_best, "BEST_OPTIMAL_CAMPAIGN_FOUND"),
        campaign_to_row(global_worst, "WORST_CAMPAIGN_FOUND"),
        {
            "case": "AVERAGE_OF_SAMPLED_BEST_ORDER_CAMPAIGNS",
            "initial_stack_mass_kg": average_mass,
            "mother_propellant_kg": results_df["mother_propellant_kg"].mean(),
            "total_tug_propellant_kg": results_df["total_tug_propellant_kg"].mean(),
            "total_propellant_kg": results_df["total_propellant_kg"].mean(),
            "mother_total_dv_m_s": results_df["mother_total_dv_m_s"].mean(),
            "average_tug_dv_m_s": results_df["average_tug_dv_m_s"].mean(),
            "target_names": "",
            "norad_ids": "",
            "mother_dvs_m_s": "",
            "tug_dvs_m_s": "",
        },
        {
            "case": "MEDIAN_OF_SAMPLED_BEST_ORDER_CAMPAIGNS",
            "initial_stack_mass_kg": median_mass,
            "mother_propellant_kg": results_df["mother_propellant_kg"].median(),
            "total_tug_propellant_kg": results_df["total_tug_propellant_kg"].median(),
            "total_propellant_kg": results_df["total_propellant_kg"].median(),
            "mother_total_dv_m_s": results_df["mother_total_dv_m_s"].median(),
            "average_tug_dv_m_s": results_df["average_tug_dv_m_s"].median(),
            "target_names": "",
            "norad_ids": "",
            "mother_dvs_m_s": "",
            "tug_dvs_m_s": "",
        },
    ]

    summary_df = pd.DataFrame(summary_rows)

    # --------------------------------------------------------
    # Save files
    # --------------------------------------------------------

    results_path = output_dir / "all_sampled_campaigns.csv"
    summary_path = output_dir / "campaign_mass_summary.csv"

    results_df.to_csv(results_path, index=False)
    summary_df.to_csv(summary_path, index=False)

    # --------------------------------------------------------
    # Plot mass distribution
    # --------------------------------------------------------

    plt.figure(figsize=(8, 5))
    plt.hist(results_df["initial_stack_mass_kg"], bins=40)
    plt.xlabel("Initial stack mass [kg]")
    plt.ylabel("Number of sampled campaigns")
    plt.title("Distribution of campaign initial mass")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "campaign_initial_mass_distribution.png", dpi=300)
    plt.show()

    plt.figure(figsize=(8, 5))
    plt.hist(results_df["total_propellant_kg"], bins=40)
    plt.xlabel("Total propellant mass [kg]")
    plt.ylabel("Number of sampled campaigns")
    plt.title("Distribution of total propellant mass")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "campaign_propellant_distribution.png", dpi=300)
    plt.show()

    # --------------------------------------------------------
    # Print results
    # --------------------------------------------------------

    print("\n=== Summary ===")
    print(summary_df[[
        "case",
        "initial_stack_mass_kg",
        "mother_propellant_kg",
        "total_tug_propellant_kg",
        "total_propellant_kg",
        "mother_total_dv_m_s",
        "average_tug_dv_m_s",
        "target_names",
        "norad_ids",
    ]].round(2))

    print("\nFiles saved in:", output_dir)
    print("- all_sampled_campaigns.csv")
    print("- campaign_mass_summary.csv")
    print("- campaign_initial_mass_distribution.png")
    print("- campaign_propellant_distribution.png")


if __name__ == "__main__":
    main()