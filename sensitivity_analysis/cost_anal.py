import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# BASELINE INPUTS
# ============================================================

mothership_inputs = {
    "mass_ss": {
        "structure_thermal_mass_kg": 353.35,
        "adcs_mass_kg": 44.29,
        "eps_mass_kg": 90.69,
        "propulsion_chemical_dry_mass_kg": 131.96,
        "ttc_cdh_mass_kg": 21.236,
        "robotic_capture_payload_mass_kg": 154.86,
    },
    "trl": {
        "structure_thermal_trl": 8,
        "adcs_trl": 9,
        "eps_trl": 8,
        "propulsion_trl": 7,
        "ttc_cdh_trl": 8,
        "robotic_capture_payload_trl": 7,
    },
    "add_params": {
        "power_bol_w": 1818,
        "propellant_volume_cm3": 1_550_000,
        "total_burn_time_s": 60_000,
        "geo_orbit_flag": 1,
    },
}

tug_inputs = {
    "mass_ss": {
        "structure_mass_kg": 60.0,
        "thermal_control_mass_kg": 14.63,
        "adcs_mass_kg": 18.415,
        "eps_mass_kg": 24.643,
        "electric_propulsion_mass_kg": 35.019,
        "ttc_mass_kg": 1.806,
        "cdh_mass_kg": 9.46,
        "capture_mechanism_mass_kg": 20.52,
    },
    "trl": {
        "structure_trl": 8,
        "thermal_trl": 8,
        "adcs_trl": 9,
        "eps_trl": 8,
        "propulsion_trl": 7,
        "ttc_trl": 8,
        "cdh_trl": 8,
        "capture_mechanism_trl": 7,
    },
    "add_params": {
        "power_bol_w": 1952,
    },
}

software_inputs = {
    "mothership_flight_software_sloc": 7608,
    "tug_flight_software_sloc": 2838,
    "ground_software_sloc": 4080,
}

program_inputs = {
    "base_year": 2026,
    "fx_eur_per_usd": 0.86,
    "n_tugs": 5,
    "learning_curve_slope": 0.95,
    "tug_recurring_fraction": 0.40,
    "operations_duration_yr": 1.0,
    "programme_margin": 0.05,
    "contractor_fee": 0.00,
}

operations_inputs = {
    "mission_operations_engineers_fte": 8,
    "mission_operations_technicians_fte": 4,
    "facility_floor_area_m2": 1000,
}

launch_inputs = {
    "falcon9_price_kusd": 89_681.81818181818,
    "kick_stage_price_kusd": 5_000,
    "n_launches": 1,
}


# ============================================================
# VARIATION RANGES
# ============================================================

N_MONTE_CARLO = 20_000
RANDOM_SEED = 42

MASS_VARIATION_FRAC = 0.20
POWER_BOL_VARIATION_FRAC = 0.20
TRL_VARIATION = 1

# SLOC is more uncertain than mass/power at early design stage.
# Triangular: min = -30%, mode = baseline, max = +80%
SLOC_MIN_FRAC = 0.70
SLOC_MODE_FRAC = 1.00
SLOC_MAX_FRAC = 1.80

# Launch price triangular distribution.
# min = predicted cheaper launch in 10 years
# mode = current baseline / most likely value
# max = current listed high value
USE_LAUNCH_PRICE_VARIATION = False
LAUNCH_PRICE_MIN_10YR_KUSD = 64_564.5  # $75M USD https://amostech.com/TechnicalPapers/2023/Poster/Shahady.pdf
LAUNCH_PRICE_MAX_TODAY_KUSD = 89_681.81818181818
LAUNCH_PRICE_MODE_KUSD = (LAUNCH_PRICE_MAX_TODAY_KUSD - LAUNCH_PRICE_MIN_10YR_KUSD) * (2/3) + LAUNCH_PRICE_MIN_10YR_KUSD

# ============================================================
# CONSTANTS
# ============================================================

inflation_factors = {
    2010: 0.9701,
    2011: 0.9850,
    2012: 1.0000,
    2013: 1.0201,
    2014: 1.0413,
    2015: 1.0630,
    2016: 1.0857,
    2017: 1.1088,
    2018: 1.1326,
    2019: 1.1568,
    2020: 1.1815,
    2021: 1.2067,
    2022: 1.2324,
    2023: 1.2588,
    2024: 1.2856,
    2025: 1.3131,
    2026: 1.3411,
}

trl_factors = {
    1: 4.0,
    2: 3.0,
    3: 2.0,
    4: 1.5,
    5: 1.3,
    6: 1.0,
    7: 0.8,
    8: 0.7,
    9: 0.5,
}

mothership_make_flags = {
    "structure_thermal": 0,
    "adcs": 0,
    "eps": 0,
    "propulsion": 0,
    "ttc_cdh": 0,
    "robotic_capture_payload": 0,
}

tug_make_flags = {
    "structure": 0,
    "thermal_control": 0,
    "adcs": 0,
    "eps": 0,
    "electric_propulsion": 0,
    "ttc": 0,
    "cdh": 0,
    "capture_mechanism": 0,
}


# ============================================================
# SMAD CER OUTPUT UNCERTAINTY VALUES
# ============================================================

ms_cer_uncertainty = {
    "structure_thermal_nr": 0.22,
    "adcs_nr": 0.44,
    "eps_nr": 0.41,
    "propulsion_nr": 0.35,

    "structure_thermal_rec": 0.21,
    "adcs_rec": 0.36,
    "eps_rec": 0.31,
    "propulsion_rec": 0.22,
    "ttc_cdh_rec": 0.18,

    "iat_nr": 0.42,
    "program_nr": 0.50,
    "age_nr": 0.37,
    "iat_rec": 0.34,
    "program_rec": 0.40,
}

tug_cer_abs_uncertainty = {
    "structure": 1097,
    "thermal_control": 119,
    "adcs": 1113,
    "eps": 910,
    "electric_propulsion": 310,
    "ttc": 629,
    "cdh": 854,
}

software_cer_uncertainty = {
    "software_mothership": 0.30,
    "software_tug": 0.30,
    "software_ground": 0.30,
}

launch_cer_uncertainty = {
    "launch_falcon9": 0.0,
    "launch_kick_stage": 0.0,
}


# ============================================================
# RANDOM GENERATOR
# ============================================================

rng = np.random.default_rng(RANDOM_SEED)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def learning_curve_factor(n, slope):
    return n ** (1 + np.log(slope) / np.log(2))


def sample_uniform_fraction(value, frac):
    return rng.uniform(value * (1 - frac), value * (1 + frac))


def sample_trl(value):
    low = max(1, value - TRL_VARIATION)
    high = min(9, value + TRL_VARIATION)
    return int(rng.integers(low, high + 1))


def sample_sloc(value):
    return int(round(
        rng.triangular(
            SLOC_MIN_FRAC * value,
            SLOC_MODE_FRAC * value,
            SLOC_MAX_FRAC * value,
        )
    ))


def sample_launch_price(launch_base):
    if USE_LAUNCH_PRICE_VARIATION:
        return rng.triangular(
            LAUNCH_PRICE_MIN_10YR_KUSD,
            LAUNCH_PRICE_MODE_KUSD,
            LAUNCH_PRICE_MAX_TODAY_KUSD,
        )
    else:
        return launch_base["falcon9_price_kusd"]


def lognormal_multiplier(rel_sigma):
    if rel_sigma is None or rel_sigma == 0:
        return 1.0

    sigma_ln = np.sqrt(np.log(1 + rel_sigma**2))
    mu_ln = -0.5 * sigma_ln**2
    return rng.lognormal(mu_ln, sigma_ln)


def apply_relative_output_uncertainty(value, rel_sigma, output_uncertainty):
    if not output_uncertainty:
        return value
    return value * lognormal_multiplier(rel_sigma)


def apply_absolute_output_uncertainty(value, abs_sigma, output_uncertainty):
    if not output_uncertainty:
        return value

    if abs_sigma is None or abs_sigma == 0:
        return value

    rel_sigma = abs_sigma / max(value, 1e-9)
    return value * lognormal_multiplier(rel_sigma)


# ============================================================
# COST MODEL
# ============================================================

def compute_mothership_cost(inputs, output_uncertainty=False):
    m = inputs["mass_ss"]
    trl = inputs["trl"]
    p = inputs["add_params"]

    costs = {}

    costs["structure_thermal_nr"] = (
        646 * m["structure_thermal_mass_kg"]**0.684
        * trl_factors[trl["structure_thermal_trl"]]
        * mothership_make_flags["structure_thermal"]
    )

    costs["adcs_nr"] = (
        324 * m["adcs_mass_kg"]
        * trl_factors[trl["adcs_trl"]]
        * mothership_make_flags["adcs"]
    )

    costs["eps_nr"] = (
        64.3 * m["eps_mass_kg"]
        * trl_factors[trl["eps_trl"]]
        * mothership_make_flags["eps"]
    )

    costs["propulsion_nr"] = (
        20 * p["propellant_volume_cm3"]**0.485
        * trl_factors[trl["propulsion_trl"]]
        * mothership_make_flags["propulsion"]
    )

    costs["ttc_cdh_nr"] = (
        26916
        * trl_factors[trl["ttc_cdh_trl"]]
        * mothership_make_flags["ttc_cdh"]
    )

    costs["robotic_capture_payload_nr"] = (
        50 * m["robotic_capture_payload_mass_kg"]
        * trl_factors[trl["robotic_capture_payload_trl"]]
        * mothership_make_flags["robotic_capture_payload"]
    )

    costs["structure_thermal_rec"] = 22.6 * m["structure_thermal_mass_kg"]
    costs["adcs_rec"] = 795 * m["adcs_mass_kg"]**0.593
    costs["eps_rec"] = 32.4 * m["eps_mass_kg"]
    costs["propulsion_rec"] = (
        29 * m["propulsion_chemical_dry_mass_kg"]
        + 0.024 * p["total_burn_time_s"]
    )
    costs["ttc_cdh_rec"] = (
        883.7 * m["ttc_cdh_mass_kg"]**0.491 * 1.13**p["geo_orbit_flag"]
    )
    costs["robotic_capture_payload_rec"] = 50 * m["robotic_capture_payload_mass_kg"]

    for key in list(costs.keys()):
        costs[key] = apply_relative_output_uncertainty(
            costs[key],
            ms_cer_uncertainty.get(key, 0.0),
            output_uncertainty,
        )

    nr_subtotal = sum(v for k, v in costs.items() if k.endswith("_nr"))
    rec_subtotal = sum(v for k, v in costs.items() if k.endswith("_rec"))

    iat_nr = 0.195 * nr_subtotal
    iat_rec = 0.124 * rec_subtotal
    program_nr = 0.414 * (nr_subtotal + iat_nr)
    program_rec = 0.320 * (rec_subtotal + iat_rec)
    age_nr = 0.421 * nr_subtotal**0.907 * 2.244 if nr_subtotal > 0 else 0
    loos_rec = 5850

    iat_nr = apply_relative_output_uncertainty(
        iat_nr, ms_cer_uncertainty["iat_nr"], output_uncertainty
    )
    iat_rec = apply_relative_output_uncertainty(
        iat_rec, ms_cer_uncertainty["iat_rec"], output_uncertainty
    )
    program_nr = apply_relative_output_uncertainty(
        program_nr, ms_cer_uncertainty["program_nr"], output_uncertainty
    )
    program_rec = apply_relative_output_uncertainty(
        program_rec, ms_cer_uncertainty["program_rec"], output_uncertainty
    )
    age_nr = apply_relative_output_uncertainty(
        age_nr, ms_cer_uncertainty["age_nr"], output_uncertainty
    )

    dev = nr_subtotal + iat_nr + program_nr + age_nr
    prod = rec_subtotal + iat_rec + program_rec + loos_rec

    return {
        "mothership_dev": dev,
        "mothership_prod": prod,
        **costs,
    }


def compute_tug_cost(inputs, program_inputs, output_uncertainty=False):
    m = inputs["mass_ss"]
    trl = inputs["trl"]

    n_tugs = program_inputs["n_tugs"]
    rec_frac = program_inputs["tug_recurring_fraction"]
    lcf = learning_curve_factor(n_tugs, program_inputs["learning_curve_slope"])

    sscm = {}

    sscm["structure"] = (
        407 + 19.3 * m["structure_mass_kg"] * np.log(m["structure_mass_kg"])
    )

    sscm["thermal_control"] = (
        335 + 5.7 * m["thermal_control_mass_kg"]**2
    )

    sscm["adcs"] = (
        1850 + 11.7 * m["adcs_mass_kg"]**2
    )

    sscm["eps"] = (
        1261 + 539 * m["eps_mass_kg"]**0.72
    )

    sscm["electric_propulsion"] = (
        89 + 3.0 * m["electric_propulsion_mass_kg"]**1.261
    )

    sscm["ttc"] = (
        486 + 55.5 * m["ttc_mass_kg"]**1.35
    )

    sscm["cdh"] = (
        658 + 75 * m["cdh_mass_kg"]**1.35
    )

    sscm["capture_mechanism"] = (
        50 * m["capture_mechanism_mass_kg"]
    )

    for key in list(sscm.keys()):
        sscm[key] = apply_absolute_output_uncertainty(
            sscm[key],
            tug_cer_abs_uncertainty.get(key, 0.0),
            output_uncertainty,
        )

    trl_map = {
        "structure": trl["structure_trl"],
        "thermal_control": trl["thermal_trl"],
        "adcs": trl["adcs_trl"],
        "eps": trl["eps_trl"],
        "electric_propulsion": trl["propulsion_trl"],
        "ttc": trl["ttc_trl"],
        "cdh": trl["cdh_trl"],
        "capture_mechanism": trl["capture_mechanism_trl"],
    }

    dev = {}
    prod = {}

    for key, value in sscm.items():
        dev[key] = (
            value
            * (1 - rec_frac)
            * trl_factors[trl_map[key]]
            * tug_make_flags[key]
        )

        prod[key] = value * rec_frac * lcf

    bus_dev = sum(dev.values())
    bus_prod = sum(prod.values())
    bus_total = bus_dev + bus_prod

    iat = 0.139 * bus_total
    program = 0.229 * bus_total
    loos = 0.061 * bus_total
    gse = 0.066 * bus_total

    total = bus_dev + bus_prod + iat + program + loos + gse

    return {
        "tug_dev": bus_dev,
        "tug_prod": bus_prod,
        "tug_wrap": iat + program + loos + gse,
        "tug_total": total,
        **{f"tug_sscm_{k}": v for k, v in sscm.items()},
    }


def compute_software_cost(software_inputs, output_uncertainty=False):
    ms_sw = software_inputs["mothership_flight_software_sloc"] * 550 / 1000
    tug_sw = software_inputs["tug_flight_software_sloc"] * 550 / 1000
    ground_sw = software_inputs["ground_software_sloc"] * 130 / 1000

    ms_sw = apply_relative_output_uncertainty(
        ms_sw,
        software_cer_uncertainty["software_mothership"],
        output_uncertainty,
    )
    tug_sw = apply_relative_output_uncertainty(
        tug_sw,
        software_cer_uncertainty["software_tug"],
        output_uncertainty,
    )
    ground_sw = apply_relative_output_uncertainty(
        ground_sw,
        software_cer_uncertainty["software_ground"],
        output_uncertainty,
    )

    return {
        "software_mothership": ms_sw,
        "software_tug": tug_sw,
        "software_ground": ground_sw,
        "software_total": ms_sw + tug_sw + ground_sw,
    }


def compute_operations_cost(software_inputs, operations_inputs, program_inputs):
    fte_eng_cost = 200
    fte_tech_cost = 150

    flight_sloc_maintained = (
        software_inputs["mothership_flight_software_sloc"]
        + software_inputs["tug_flight_software_sloc"]
    )

    ground_sloc_maintained = software_inputs["ground_software_sloc"]

    space_sw_maintenance = flight_sloc_maintained / 16000 * fte_eng_cost
    mission_operations = (
        operations_inputs["mission_operations_engineers_fte"] * fte_eng_cost
        + operations_inputs["mission_operations_technicians_fte"] * fte_tech_cost
    )
    ground_sw_maintenance = ground_sloc_maintained / 28200 * fte_eng_cost
    ground_hw_maintenance = 0.0
    facilities = operations_inputs["facility_floor_area_m2"] * 1.25

    pmse = 0.10 * (
        space_sw_maintenance
        + mission_operations
        + ground_sw_maintenance
        + ground_hw_maintenance
        + facilities
    )

    annual = (
        pmse
        + space_sw_maintenance
        + mission_operations
        + ground_sw_maintenance
        + ground_hw_maintenance
        + facilities
    )

    return {
        "operations_annual": annual,
        "operations_total": annual * program_inputs["operations_duration_yr"],
    }


def compute_launch_cost(launch_inputs, program_inputs, output_uncertainty=False):
    infl = inflation_factors[program_inputs["base_year"]] / inflation_factors[2010]

    launch_2010 = (
        launch_inputs["falcon9_price_kusd"]
        * launch_inputs["n_launches"]
        / infl
    )

    kick_2010 = launch_inputs["kick_stage_price_kusd"] / infl

    launch_2010 = apply_relative_output_uncertainty(
        launch_2010,
        launch_cer_uncertainty["launch_falcon9"],
        output_uncertainty,
    )

    kick_2010 = apply_relative_output_uncertainty(
        kick_2010,
        launch_cer_uncertainty["launch_kick_stage"],
        output_uncertainty,
    )

    return {
        "launch_falcon9": launch_2010,
        "launch_kick_stage": kick_2010,
        "launch_total": launch_2010 + kick_2010,
    }


def compute_total_cost(
    mothership_inputs,
    tug_inputs,
    software_inputs,
    operations_inputs,
    launch_inputs,
    program_inputs,
    output_uncertainty=False,
):
    ms = compute_mothership_cost(mothership_inputs, output_uncertainty)
    tug = compute_tug_cost(tug_inputs, program_inputs, output_uncertainty)
    sw = compute_software_cost(software_inputs, output_uncertainty)
    launch = compute_launch_cost(launch_inputs, program_inputs, output_uncertainty)
    ops = compute_operations_cost(software_inputs, operations_inputs, program_inputs)

    acquisition_2010 = (
        ms["mothership_dev"]
        + ms["mothership_prod"]
        + tug["tug_dev"]
        + tug["tug_prod"]
        + tug["tug_wrap"]
        + sw["software_total"]
        + launch["launch_total"]
    )

    margin = program_inputs["programme_margin"] * acquisition_2010
    fee = program_inputs["contractor_fee"] * acquisition_2010

    total_2010 = acquisition_2010 + margin + fee + ops["operations_total"]

    infl = inflation_factors[program_inputs["base_year"]] / inflation_factors[2010]
    total_meur = total_2010 * infl * program_inputs["fx_eur_per_usd"] / 1000

    return {
        "total_fy2010_kusd": total_2010,
        "total_base_year_meur": total_meur,
        "acquisition_fy2010_kusd": acquisition_2010,
        "operations_fy2010_kusd": ops["operations_total"],
        "margin_fy2010_kusd": margin,
        **ms,
        **tug,
        **sw,
        **launch,
        **ops,
    }


# ============================================================
# INPUT UNCERTAINTY MONTE CARLO
# ============================================================

def sample_input_uncertainty(
    mothership_base,
    tug_base,
    software_base,
    launch_base,
):
    ms = copy.deepcopy(mothership_base)
    tug = copy.deepcopy(tug_base)
    sw = copy.deepcopy(software_base)
    launch = copy.deepcopy(launch_base)

    for key, value in ms["mass_ss"].items():
        ms["mass_ss"][key] = sample_uniform_fraction(value, MASS_VARIATION_FRAC)

    for key, value in tug["mass_ss"].items():
        tug["mass_ss"][key] = sample_uniform_fraction(value, MASS_VARIATION_FRAC)

    for key, value in ms["trl"].items():
        ms["trl"][key] = sample_trl(value)

    for key, value in tug["trl"].items():
        tug["trl"][key] = sample_trl(value)

    ms["add_params"]["power_bol_w"] = sample_uniform_fraction(
        ms["add_params"]["power_bol_w"],
        POWER_BOL_VARIATION_FRAC,
    )

    tug["add_params"]["power_bol_w"] = sample_uniform_fraction(
        tug["add_params"]["power_bol_w"],
        POWER_BOL_VARIATION_FRAC,
    )

    for key, value in sw.items():
        sw[key] = sample_sloc(value)

    launch["falcon9_price_kusd"] = sample_launch_price(launch_base)

    return ms, tug, sw, launch


def run_input_uncertainty_monte_carlo(n=N_MONTE_CARLO):
    rows = []

    for _ in range(n):
        ms_i, tug_i, sw_i, launch_i = sample_input_uncertainty(
            mothership_inputs,
            tug_inputs,
            software_inputs,
            launch_inputs,
        )

        result = compute_total_cost(
            ms_i,
            tug_i,
            sw_i,
            operations_inputs,
            launch_i,
            program_inputs,
            output_uncertainty=False,
        )

        rows.append(result)

    return pd.DataFrame(rows)


# ============================================================
# OUTPUT / SMAD SEE UNCERTAINTY MONTE CARLO
# ============================================================

def run_output_uncertainty_monte_carlo(n=N_MONTE_CARLO):
    rows = []

    for _ in range(n):
        result = compute_total_cost(
            mothership_inputs,
            tug_inputs,
            software_inputs,
            operations_inputs,
            launch_inputs,
            program_inputs,
            output_uncertainty=True,
        )

        rows.append(result)

    return pd.DataFrame(rows)


# ============================================================
# EXECUTE
# ============================================================

baseline = compute_total_cost(
    mothership_inputs,
    tug_inputs,
    software_inputs,
    operations_inputs,
    launch_inputs,
    program_inputs,
    output_uncertainty=False,
)

mc_input_uncertainty = run_input_uncertainty_monte_carlo(N_MONTE_CARLO)
mc_output_uncertainty = run_output_uncertainty_monte_carlo(N_MONTE_CARLO)

print("Baseline total cost [MEUR]:")
print(f"{baseline['total_base_year_meur']:.2f}")

print("\nInput uncertainty Monte Carlo [MEUR]:")
print(
    mc_input_uncertainty["total_base_year_meur"].describe(
        percentiles=[0.05, 0.10, 0.50, 0.90, 0.95]
    )
)

print("\nOutput / SMAD SEE uncertainty Monte Carlo [MEUR]:")
print(
    mc_output_uncertainty["total_base_year_meur"].describe(
        percentiles=[0.05, 0.10, 0.50, 0.90, 0.95]
    )
)


# ============================================================
# PLOTS
# ============================================================

from matplotlib.ticker import PercentFormatter

plt.figure(figsize=(8, 5))
plt.hist(
    mc_input_uncertainty["total_base_year_meur"],
    bins=50,
    weights=np.ones(len(mc_input_uncertainty)) / len(mc_input_uncertainty) * 100
)
plt.axvline(baseline["total_base_year_meur"], linestyle="--", label="Baseline", color="tab:orange")
plt.xlabel("Total mission cost [M€]")
plt.ylabel("Percentage of simulations [%]")
plt.title("Input Uncertainty Monte Carlo")
plt.gca().yaxis.set_major_formatter(PercentFormatter(xmax=100))
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(
    mc_output_uncertainty["total_base_year_meur"],
    bins=50,
    weights=np.ones(len(mc_output_uncertainty)) / len(mc_output_uncertainty) * 100
)
plt.axvline(baseline["total_base_year_meur"], linestyle="--", label="Baseline", color="tab:orange")
plt.xlabel("Total mission cost [M€]")
plt.ylabel("Percentage of simulations [%]")
plt.title("SMAD CER Output Uncertainty Monte Carlo")
plt.gca().yaxis.set_major_formatter(PercentFormatter(xmax=100))
plt.legend()
plt.tight_layout()
plt.show()