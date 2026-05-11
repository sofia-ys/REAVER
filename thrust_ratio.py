import math
from typing import Optional, Dict

# ============================================================
# CONSTANTS
# ============================================================
MU_EARTH = 3.986004418e14      # Earth gravitational parameter [m^3/s^2]
R_EARTH = 6378.137e3           # Earth equatorial radius [m]
G0 = 9.80665                   # Standard gravity [m/s^2]


# ============================================================
# SMALL UTILITY FUNCTIONS
# ============================================================
def clamp(value: float, lower: float = -1.0, upper: float = 1.0) -> float:
    """Avoid numerical issues in acos caused by tiny floating-point errors."""
    return max(lower, min(upper, value))


def angle_difference_rad(angle_1_deg: float, angle_2_deg: float) -> float:
    """
    Smallest absolute difference between two angles.

    Example:
        350 deg to 10 deg gives 20 deg, not 340 deg.

    Returns:
        angle difference [rad], in the range [0, pi]
    """
    diff_deg = (angle_2_deg - angle_1_deg + 180.0) % 360.0 - 180.0
    return math.radians(abs(diff_deg))


def vector_delta_v(v_from: float, v_to: float, angle_rad: float = 0.0) -> float:
    """
    Magnitude of the velocity change between two velocity vectors.

    v_from and v_to are the magnitudes of the initial and final velocity vectors.
    angle_rad is the angle between the two velocity directions.
    """
    value = v_from**2 + v_to**2 - 2.0 * v_from * v_to * math.cos(angle_rad)
    return math.sqrt(max(0.0, value))


def safe_ratio(numerator: float, denominator: float) -> Optional[float]:
    """Avoid division by zero."""
    if abs(denominator) < 1e-12:
        return None
    return numerator / denominator


# ============================================================
# ORBIT GEOMETRY
# ============================================================
def orbital_plane_angle_rad(
    i0_deg: float,
    if_deg: float,
    raan0_deg: float,
    raanf_deg: float,
) -> float:
    """
    Computes the true 3D angle between the initial and final orbital planes.

    This includes BOTH:
        - inclination change
        - RAAN change

    This is the correct angle to use for a combined plane change.
    """
    i0 = math.radians(i0_deg)
    i_f = math.radians(if_deg)
    delta_raan = angle_difference_rad(raan0_deg, raanf_deg)

    cos_theta = (
        math.cos(i0) * math.cos(i_f)
        + math.sin(i0) * math.sin(i_f) * math.cos(delta_raan)
    )

    return math.acos(clamp(cos_theta))


def hohmann_transfer_quantities(r1_m: float, r2_m: float) -> Dict[str, float]:
    """
    Circular-to-circular Hohmann transfer quantities.
    """
    a_t = 0.5 * (r1_m + r2_m)

    v_circ_1 = math.sqrt(MU_EARTH / r1_m)
    v_circ_2 = math.sqrt(MU_EARTH / r2_m)

    v_trans_1 = math.sqrt(MU_EARTH * (2.0 / r1_m - 1.0 / a_t))
    v_trans_2 = math.sqrt(MU_EARTH * (2.0 / r2_m - 1.0 / a_t))

    coast_time_s = math.pi * math.sqrt(a_t**3 / MU_EARTH)

    return {
        "a_transfer_m": a_t,
        "v_circ_1_m_s": v_circ_1,
        "v_circ_2_m_s": v_circ_2,
        "v_trans_1_m_s": v_trans_1,
        "v_trans_2_m_s": v_trans_2,
        "coast_time_s": coast_time_s,
    }


# ============================================================
# PROPULSION / MASS MODEL
# ============================================================
def rocket_equation_from_initial_mass(
    initial_mass_kg: float,
    delta_v_m_s: float,
    isp_s: Optional[float],
) -> Dict[str, Optional[float]]:
    """
    Uses the Tsiolkovsky rocket equation.

    If isp_s is None, propellant mass cannot be estimated and the function
    returns the initial mass as the representative mass.
    """
    if initial_mass_kg <= 0.0:
        raise ValueError("initial_mass_kg must be positive.")

    if delta_v_m_s < 0.0:
        raise ValueError("delta_v_m_s cannot be negative.")

    if isp_s is None:
        return {
            "initial_mass_kg": initial_mass_kg,
            "final_mass_kg": None,
            "propellant_mass_kg": None,
            "propellant_fraction": None,
            "average_mass_kg": initial_mass_kg,
        }

    if isp_s <= 0.0:
        raise ValueError("Isp must be positive if provided.")

    final_mass_kg = initial_mass_kg * math.exp(-delta_v_m_s / (isp_s * G0))
    propellant_mass_kg = initial_mass_kg - final_mass_kg
    average_mass_kg = 0.5 * (initial_mass_kg + final_mass_kg)

    return {
        "initial_mass_kg": initial_mass_kg,
        "final_mass_kg": final_mass_kg,
        "propellant_mass_kg": propellant_mass_kg,
        "propellant_fraction": propellant_mass_kg / initial_mass_kg,
        "average_mass_kg": average_mass_kg,
    }


# ============================================================
# HIGH-THRUST MODEL
# ============================================================
def high_thrust_optimised_hohmann(
    r1_m: float,
    r2_m: float,
    plane_change_rad: float,
) -> Dict[str, float]:
    """
    Impulsive high-thrust reference case.

    Model:
        - Hohmann transfer for altitude change.
        - Total plane change is split between burn 1 and burn 2.
        - The split is numerically optimised to minimise total Delta V.

    This is more accurate than forcing the full plane change to occur at
    only one burn.
    """
    v1 = math.sqrt(MU_EARTH / r1_m)

    # Same-radius case: no Hohmann coast, only pure plane change.
    if abs(r2_m - r1_m) < 1e-9:
        dv_total = 2.0 * v1 * math.sin(plane_change_rad / 2.0)
        return {
            "delta_v_1_m_s": 0.0,
            "delta_v_2_m_s": dv_total,
            "delta_v_total_m_s": dv_total,
            "plane_change_burn_1_deg": 0.0,
            "plane_change_burn_2_deg": math.degrees(plane_change_rad),
            "coast_time_s": 0.0,
        }

    h = hohmann_transfer_quantities(r1_m, r2_m)

    v_circ_1 = h["v_circ_1_m_s"]
    v_circ_2 = h["v_circ_2_m_s"]
    v_trans_1 = h["v_trans_1_m_s"]
    v_trans_2 = h["v_trans_2_m_s"]

    def total_delta_v_for_split(alpha: float) -> float:
        """
        alpha is the amount of plane change assigned to burn 1.
        The remaining plane change is assigned to burn 2.
        """
        beta = plane_change_rad - alpha

        dv_1 = vector_delta_v(v_circ_1, v_trans_1, alpha)
        dv_2 = vector_delta_v(v_trans_2, v_circ_2, beta)

        return dv_1 + dv_2

    # Golden-section search over alpha in [0, plane_change_rad].
    # This avoids scipy and keeps the script self-contained.
    lower = 0.0
    upper = plane_change_rad

    golden_ratio_part = (math.sqrt(5.0) - 1.0) / 2.0

    c = upper - golden_ratio_part * (upper - lower)
    d = lower + golden_ratio_part * (upper - lower)

    f_c = total_delta_v_for_split(c)
    f_d = total_delta_v_for_split(d)

    for _ in range(100):
        if upper - lower < 1e-12:
            break

        if f_c > f_d:
            lower = c
            c = d
            f_c = f_d
            d = lower + golden_ratio_part * (upper - lower)
            f_d = total_delta_v_for_split(d)
        else:
            upper = d
            d = c
            f_d = f_c
            c = upper - golden_ratio_part * (upper - lower)
            f_c = total_delta_v_for_split(c)

    alpha_opt = 0.5 * (lower + upper)
    beta_opt = plane_change_rad - alpha_opt

    dv_1 = vector_delta_v(v_circ_1, v_trans_1, alpha_opt)
    dv_2 = vector_delta_v(v_trans_2, v_circ_2, beta_opt)

    return {
        "delta_v_1_m_s": dv_1,
        "delta_v_2_m_s": dv_2,
        "delta_v_total_m_s": dv_1 + dv_2,
        "plane_change_burn_1_deg": math.degrees(alpha_opt),
        "plane_change_burn_2_deg": math.degrees(beta_opt),
        "coast_time_s": h["coast_time_s"],
    }


# ============================================================
# LOW-THRUST MODEL
# ============================================================
def low_thrust_edelbaum(
    r1_m: float,
    r2_m: float,
    plane_change_rad: float,
) -> Dict[str, float]:
    """
    First-order Edelbaum-style low-thrust estimate.

    Assumptions:
        - initial and final orbits are circular or near-circular
        - continuous low-thrust spiral
        - altitude change and plane change occur simultaneously
        - eccentricity is not targeted explicitly
        - perturbations and shadowing are not included

    Formula:
        Delta V = sqrt(V1^2 + V2^2 - 2 V1 V2 cos(pi/2 * theta))

    where theta is the true 3D plane-change angle.
    """
    v1 = math.sqrt(MU_EARTH / r1_m)
    v2 = math.sqrt(MU_EARTH / r2_m)

    delta_v = vector_delta_v(v1, v2, 0.5 * math.pi * plane_change_rad)

    return {
        "delta_v_total_m_s": delta_v,
    }


# ============================================================
# MAIN ESTIMATOR
# ============================================================
def estimate_transfer(
    h0_km: float,
    hf_km: float,
    i0_deg: float,
    if_deg: float,
    raan0_deg: float,
    raanf_deg: float,
    initial_mass_kg: float,
    low_thrust_N: float,
    low_isp_s: Optional[float] = None,
    high_isp_s: Optional[float] = None,
    high_thrust_N: Optional[float] = None,
) -> Dict[str, Optional[float]]:
    """
    Estimates low-thrust and high-thrust Delta V, transfer time, propellant
    mass, and useful scaling ratios.
    """
    if low_thrust_N <= 0.0:
        raise ValueError("low_thrust_N must be positive.")

    if high_thrust_N is not None and high_thrust_N <= 0.0:
        raise ValueError("high_thrust_N must be positive if provided.")

    r1_m = R_EARTH + h0_km * 1e3
    r2_m = R_EARTH + hf_km * 1e3

    if r1_m <= 0.0 or r2_m <= 0.0:
        raise ValueError("Orbit radius must be positive. Check h0_km and hf_km.")

    inclination_change_deg = abs(if_deg - i0_deg)
    raan_change_deg = math.degrees(angle_difference_rad(raan0_deg, raanf_deg))
    plane_change_rad = orbital_plane_angle_rad(i0_deg, if_deg, raan0_deg, raanf_deg)

    v1 = math.sqrt(MU_EARTH / r1_m)
    v2 = math.sqrt(MU_EARTH / r2_m)

    low = low_thrust_edelbaum(r1_m, r2_m, plane_change_rad)
    high = high_thrust_optimised_hohmann(r1_m, r2_m, plane_change_rad)

    low_mass = rocket_equation_from_initial_mass(
        initial_mass_kg=initial_mass_kg,
        delta_v_m_s=low["delta_v_total_m_s"],
        isp_s=low_isp_s,
    )

    high_mass = rocket_equation_from_initial_mass(
        initial_mass_kg=initial_mass_kg,
        delta_v_m_s=high["delta_v_total_m_s"],
        isp_s=high_isp_s,
    )

    # Low-thrust time estimate.
    # If Isp is provided, use average mass during the burn.
    # If Isp is not provided, average_mass_kg = initial_mass_kg.
    low_average_mass_kg = low_mass["average_mass_kg"]
    low_acceleration_average_m_s2 = low_thrust_N / low_average_mass_kg
    low_transfer_time_s = low["delta_v_total_m_s"] / low_acceleration_average_m_s2

    # High-thrust finite burn duration is optional.
    # The main high-thrust transfer time is still the Hohmann coast time.
    if high_thrust_N is None:
        high_average_acceleration_m_s2 = None
        high_burn_time_s = None
        high_total_time_with_burns_s = high["coast_time_s"]
    else:
        high_average_mass_kg = high_mass["average_mass_kg"]
        high_average_acceleration_m_s2 = high_thrust_N / high_average_mass_kg
        high_burn_time_s = high["delta_v_total_m_s"] / high_average_acceleration_m_s2
        high_total_time_with_burns_s = high["coast_time_s"] + high_burn_time_s

    return {
        # Input geometry
        "r1_m": r1_m,
        "r2_m": r2_m,
        "h0_km": h0_km,
        "hf_km": hf_km,
        "v1_circular_m_s": v1,
        "v2_circular_m_s": v2,
        "inclination_change_deg": inclination_change_deg,
        "raan_change_deg": raan_change_deg,
        "combined_plane_change_angle_deg": math.degrees(plane_change_rad),

        # Low-thrust
        "low_delta_v_total_m_s": low["delta_v_total_m_s"],
        "low_thrust_N": low_thrust_N,
        "low_isp_s": low_isp_s,
        "low_initial_mass_kg": low_mass["initial_mass_kg"],
        "low_final_mass_kg": low_mass["final_mass_kg"],
        "low_average_mass_kg": low_mass["average_mass_kg"],
        "low_propellant_mass_kg": low_mass["propellant_mass_kg"],
        "low_propellant_fraction": low_mass["propellant_fraction"],
        "low_average_acceleration_m_s2": low_acceleration_average_m_s2,
        "low_transfer_time_s": low_transfer_time_s,
        "low_transfer_time_hours": low_transfer_time_s / 3600.0,
        "low_transfer_time_days": low_transfer_time_s / 86400.0,

        # High-thrust
        "high_delta_v_1_m_s": high["delta_v_1_m_s"],
        "high_delta_v_2_m_s": high["delta_v_2_m_s"],
        "high_delta_v_total_m_s": high["delta_v_total_m_s"],
        "high_plane_change_burn_1_deg": high["plane_change_burn_1_deg"],
        "high_plane_change_burn_2_deg": high["plane_change_burn_2_deg"],
        "high_isp_s": high_isp_s,
        "high_initial_mass_kg": high_mass["initial_mass_kg"],
        "high_final_mass_kg": high_mass["final_mass_kg"],
        "high_average_mass_kg": high_mass["average_mass_kg"],
        "high_propellant_mass_kg": high_mass["propellant_mass_kg"],
        "high_propellant_fraction": high_mass["propellant_fraction"],
        "high_coast_time_s": high["coast_time_s"],
        "high_coast_time_hours": high["coast_time_s"] / 3600.0,
        "high_coast_time_days": high["coast_time_s"] / 86400.0,
        "high_thrust_N": high_thrust_N,
        "high_average_acceleration_m_s2": high_average_acceleration_m_s2,
        "high_burn_time_s": high_burn_time_s,
        "high_total_time_with_burns_s": high_total_time_with_burns_s,
        "high_total_time_with_burns_hours": high_total_time_with_burns_s / 3600.0,
        "high_total_time_with_burns_days": high_total_time_with_burns_s / 86400.0,

        # Ratios
        "delta_v_ratio_low_to_high": safe_ratio(
            low["delta_v_total_m_s"],
            high["delta_v_total_m_s"],
        ),
        "transfer_time_ratio_low_to_high_coast_only": safe_ratio(
            low_transfer_time_s,
            high["coast_time_s"],
        ),
        "propellant_mass_ratio_low_to_high": (
            safe_ratio(low_mass["propellant_mass_kg"], high_mass["propellant_mass_kg"])
            if low_mass["propellant_mass_kg"] is not None and high_mass["propellant_mass_kg"] is not None
            else None
        ),
    }


# ============================================================
# PRINTING
# ============================================================
def print_value(label: str, value: Optional[float], unit: str = "", decimals: int = 3) -> None:
    if value is None:
        print(f"{label}: not calculated")
    else:
        print(f"{label}: {value:.{decimals}f} {unit}".rstrip())


def print_results(result: Dict[str, Optional[float]]) -> None:
    print("\nTRANSFER GEOMETRY")
    print("-----------------")
    print_value("Initial circular velocity", result["v1_circular_m_s"], "m/s")
    print_value("Final circular velocity", result["v2_circular_m_s"], "m/s")
    print_value("Inclination change", result["inclination_change_deg"], "deg")
    print_value("RAAN change", result["raan_change_deg"], "deg")
    print_value("Combined 3D plane-change angle", result["combined_plane_change_angle_deg"], "deg")

    print("\nLOW-THRUST / ELECTRIC ESTIMATE")
    print("------------------------------")
    print_value("Low-thrust total Delta V", result["low_delta_v_total_m_s"], "m/s")
    print_value("Low-thrust Isp", result["low_isp_s"], "s")
    print_value("Low-thrust propellant mass", result["low_propellant_mass_kg"], "kg")
    print_value("Low-thrust propellant fraction", result["low_propellant_fraction"], "", decimals=4)
    print_value("Low-thrust final mass", result["low_final_mass_kg"], "kg")
    print_value("Mass used for time estimate", result["low_average_mass_kg"], "kg")
    print_value("Average acceleration", result["low_average_acceleration_m_s2"], "m/s^2", decimals=8)
    print_value("Low-thrust transfer time", result["low_transfer_time_days"], "days")

    print("\nHIGH-THRUST / CHEMICAL IMPULSIVE ESTIMATE")
    print("-----------------------------------------")
    print_value("High-thrust burn 1 Delta V", result["high_delta_v_1_m_s"], "m/s")
    print_value("High-thrust burn 2 Delta V", result["high_delta_v_2_m_s"], "m/s")
    print_value("High-thrust total Delta V", result["high_delta_v_total_m_s"], "m/s")
    print_value("Plane change assigned to burn 1", result["high_plane_change_burn_1_deg"], "deg")
    print_value("Plane change assigned to burn 2", result["high_plane_change_burn_2_deg"], "deg")
    print_value("High-thrust Isp", result["high_isp_s"], "s")
    print_value("High-thrust propellant mass", result["high_propellant_mass_kg"], "kg")
    print_value("High-thrust propellant fraction", result["high_propellant_fraction"], "", decimals=4)
    print_value("High-thrust final mass", result["high_final_mass_kg"], "kg")
    print_value("Hohmann coast time", result["high_coast_time_days"], "days")
    print_value("Optional finite burn time", result["high_burn_time_s"], "s")

    print("\nRATIOS")
    print("------")
    print_value("Delta V ratio, low / high", result["delta_v_ratio_low_to_high"], "")
    print_value(
        "Transfer time ratio, low / high coast-only",
        result["transfer_time_ratio_low_to_high_coast_only"],
        "",
    )
    print_value(
        "Propellant mass ratio, low / high",
        result["propellant_mass_ratio_low_to_high"],
        "",
    )


# ============================================================
# USER INPUTS
# ============================================================
if __name__ == "__main__":

    result = estimate_transfer(
        # Orbit geometry
        h0_km=35700.0,
        hf_km=36500.0,

        i0_deg=7.0,
        if_deg=15.0,

        raan0_deg=0.0,
        raanf_deg=70.0,

        # Spacecraft mass
        # Treat this as initial/wet mass if Isp is provided.
        # If Isp is None, it is treated as constant representative mass.
        initial_mass_kg=1500.0,

        # Low-thrust / electric propulsion
        low_thrust_N=0.23,
        low_isp_s=1600.0,      # Put None if you do not want propellant mass estimate

        # High-thrust / chemical reference
        high_isp_s=320.0,      # Put None if you do not want propellant mass estimate
        high_thrust_N=None,    # Optional; e.g. 400.0 N if you want finite burn time
    )

    print_results(result)