"""
=============================================================================
SPACECRAFT SOLAR ARRAY SIZING TOOL
=============================================================================
Based on standard spacecraft power system sizing methodology (SMAD / ECSS).

Sections:
  1. INPUTS          – all user-configurable parameters
  2. EQUATIONS       – each relation implemented as a named function
  3. MAIN SIZING     – orchestrates the calculation sequence
  4. OUTPUT          – structured print of results
=============================================================================
"""

import math


# =============================================================================
# 1. INPUTS
# =============================================================================

class SolarArrayInputs:
    """
    All inputs required for solar array and battery sizing.
    Units are SI unless otherwise noted.
    """

    # --- Solar / environment ---
    S_e: float = 1368.0          # Solar constant at 1 AU            [W/m²]
    theta_deg: float = 23.5      # Max Sun incidence angle (worst-case, SMAD) [deg]

    # --- Cell properties ---
    mu_cell: float = 0.33       # Cell efficiency (BOL, from datasheet) [-]
    I_d: float = 0.77            # Inherent degradation factor          [-]
                                 # (typical 0.77 for triple-junction SMAD)

    # --- Degradation factors ---
    # f_temp: float = 0.90         # Temperature & radiation loss factor (datasheet) [-]
    deg_per_year: float = 0.0275 # Annual radiation degradation rate (SMAD)          [-/yr]
    lifetime_yr: float = 5.0     # Mission lifetime                          [yr]

    # --- Mission power profile ---
    P_e: float = 150.0            # Power demand during eclipse              [W]
    P_d: float = 1300.0            # Power demand during daylight             [W]
    T_e: float = 72.0*60.0       # Eclipse duration (72 min)                [s]
    T_d: float = 24*60*60 - T_e     # Daylight duration                        [s]

    # --- Power conditioning efficiencies ---
    mu_har: float  = 0.98        # Harness efficiency         [-]
    mu_mppt: float = 1        # MPPT efficiency included in PCDU efficiency
    mu_pcdu: float = 0.95        # PCDU efficiency            [-]
    mu_bat: float  = 0.90        # Battery (discharge) [-]
    #TODO: REVISE

    # --- Battery ---
    DOD: float = 0.80            # Depth of discharge (conservative)    [-]

    # --- Solar array mass model (specific area densities) ---
    m_density = None # [kg/m&2]
    m_sp = 180 # [W/kg]
    #TODO: REVISE

class SolarArrayInputsMS(SolarArrayInputs):
    """
    All inputs required for solar array and battery sizing.
    Units are SI unless otherwise noted.
    """

    # --- Solar / environment ---
    S_e: float = 1368.0          # Solar constant at 1 AU            [W/m²]
    theta_deg: float = 23.5      # Max Sun incidence angle (worst-case, SMAD) [deg]

    # --- Cell properties ---
    mu_cell: float = 0.313       # Cell efficiency (BOL, from datasheet) [-]
    I_d: float = 0.77            # Inherent degradation factor          [-]
                                 # (typical 0.77 for triple-junction) #TODO: REVISE

    # --- Degradation factors ---
    deg_per_year: float = 0.0275 # Annual radiation degradation rate (SMAD)          [-/yr]
    lifetime_yr: float = 5.0     # Mission lifetime                          [yr]

    # --- Mission power profile ---
    P_e: float = 233            # Power demand during eclipse              [W]
    P_d: float = 1445.0            # Power demand during daylight             [W]
    T_e: float = 72*60       # Eclipse duration (72 min)                [s]
    T_d: float = 24*60*60 - T_e     # Daylight duration                        [s]
    #TODO: REVISE

    # --- Power conditioning efficiencies ---
    mu_har: float  = 0.98        # Harness efficiency         [-] #TODO REVISE
    mu_mppt: float = 1        # MPPT efficiency            [-]
    mu_pcdu: float = 0.97        # PCDU efficiency            [-]
    mu_bat: float  = 0.90        # Battery (charge/discharge) [-]

    # --- Battery ---
    DOD: float = 0.80            # Depth of discharge    [-]

    # --- Solar array mass model (specific area densities) ---
    m_sp = None # [W/kg]
    m_density = 7.9 # [kg/m2]


# =============================================================================
# 2. EQUATIONS  (each relation is an explicit, documented function)
# =============================================================================

def calc_lifetime_degradation(deg_per_year: float, lifetime_yr: float) -> float:
    """
    Lifetime degradation factor due to radiation.

    Relation:
        L_d = (1 - deg_per_year) ^ lifetime_yr

    Args:
        deg_per_year : Annual fractional degradation rate  [-/yr]
        lifetime_yr  : Mission lifetime                    [yr]

    Returns:
        L_d : Lifetime degradation factor  [-]
    """
    L_d = (1.0 - deg_per_year) ** lifetime_yr
    return L_d

def calc_eol_power_density(
    S_e: float,
    mu_cell: float,
    I_d: float,
    theta_deg: float,
    L_d: float
) -> float:
    """
    End-of-life (EOL) power density of the solar array.

    Relation:
        P_EOL = S_e × mu_cell × I_d × cos(theta) × L_d

    Args:
        S_e       : Solar constant                  [W/m²]
        mu_cell   : Cell efficiency                 [-]
        I_d       : Inherent degradation            [-]
        theta_deg : Sun incidence angle             [deg]
        L_d       : Lifetime degradation factor     [-]

    Returns:
        P_EOL : EOL power density  [W/m²]
    """
    cos_theta = math.cos(math.radians(theta_deg))
    P_EOL = S_e * mu_cell * I_d * cos_theta * L_d
    return P_EOL


def calc_path_efficiencies(
    mu_har: float,
    mu_mppt: float,
    mu_pcdu: float,
    mu_bat: float
) -> tuple[float, float]:
    """
    Compute end-to-end efficiencies for each load path.

    Daylight path  (array → MPPT → PCDU → load):
        mu_d = mu_har × mu_mppt × mu_pcdu

    Eclipse path   (array → PCDU → battery → PCDU → load):
        mu_e = mu_har × mu_pcdu × mu_bat

    Args:
        mu_har  : Harness efficiency   [-]
        mu_mppt : MPPT efficiency      [-]
        mu_pcdu : PCDU efficiency      [-]
        mu_bat  : Battery efficiency   [-]

    Returns:
        (mu_d, mu_e) : Daylight and eclipse efficiencies  [-]
    """
    mu_d = mu_har * mu_mppt * mu_pcdu
    mu_e = mu_har * mu_pcdu * mu_bat
    return mu_d, mu_e


def calc_required_array_power(
    P_e: float,
    T_e: float,
    mu_e: float,
    P_d: float,
    T_d: float,
    mu_d: float
) -> float:
    """
    Required solar array output power to sustain the mission power budget.

    The array must supply:
      - Daylight load directly
      - Eclipse load via the battery (charged during daylight)

    Relation:
        P_sa = (P_e × T_e) / (mu_e × T_d)  +  P_d / mu_d

    Args:
        P_e   : Eclipse power demand    [W]
        T_e   : Eclipse duration        [s]
        mu_e  : Eclipse path efficiency [-]
        P_d   : Daylight power demand   [W]
        T_d   : Daylight duration       [s]
        mu_d  : Daylight path efficiency[-]

    Returns:
        P_sa : Required array power  [W]
    """
    P_sa = (P_e * T_e) / (mu_e * T_d) + (P_d / mu_d)
    return P_sa


def calc_array_area(P_sa: float, P_EOL: float) -> float:
    """
    Minimum solar array area.

    Relation:
        A_sa = P_sa / P_EOL

    Args:
        P_sa  : Required array power   [W]
        P_EOL : EOL power density      [W/m²]

    Returns:
        A_sa : Solar array area  [m²]
    """
    A_sa = P_sa / P_EOL
    return A_sa


def calc_array_mass(
    # A_sa: float,
    # rho_cell: float,
    # m_yoke: float,
    # m_deploy: float
    A_sa: float,
    P_sa: float,
    m_density: float,
    m_sp: float,
) -> dict[str, float]:
    """
    Solar array mass breakdown.

    Relation:
        m_cells  = rho_cell × A_sa
        m_total  = m_cells + m_yoke + m_deploy

    Args:
        P_sa    : Solar array power   [W]
        m_sp    : specific mass [W/kg]

    Returns:
        dict with keys: m_cells, m_yoke, m_deploy, m_total  [kg]
    """
    if m_density:
        m_total = A_sa * m_density
    elif m_sp:
        m_total = P_sa / m_sp
    else:
        raise ValueError

    return {
        "m_total":  m_total,
    }


def calc_battery_energy(
    P_e: float,
    T_e: float,
    DOD: float,
    mu_bat: float
) -> float:
    """
    Required battery energy capacity.

    Relation:
        E_bat = (P_e × T_e) / (DOD × mu_bat)

    Args:
        P_e    : Eclipse power demand    [W]
        T_e    : Eclipse duration        [s]
        DOD    : Depth of discharge      [-]
        mu_bat : Battery efficiency      [-]

    Returns:
        E_bat : Battery energy  [J]  (divide by 3600 for Wh)
    """
    E_bat = (P_e * T_e) / (DOD * mu_bat)
    return E_bat


# =============================================================================
# 3. MAIN SIZING  – orchestrate all relations in correct dependency order
# =============================================================================

def run_sizing(inp: SolarArrayInputs) -> dict:
    """
    Execute the full solar array sizing sequence.

    Dependency order:
        L_d  →  P_EOL
        mu_d, mu_e  →  P_sa
        P_sa / P_EOL  →  A_sa
        A_sa  →  m_sa
        P_e, T_e, DOD, mu_bat  →  E_bat

    Args:
        inp : SolarArrayInputs instance

    Returns:
        results : dict containing all intermediate and final values
    """

    # Step 1 – Lifetime degradation
    L_d = calc_lifetime_degradation(inp.deg_per_year, inp.lifetime_yr)

    # Step 2 – EOL power density
    P_EOL = calc_eol_power_density(
        inp.S_e, inp.mu_cell, inp.I_d, inp.theta_deg, L_d
    )

    # Step 3 – Load-path efficiencies
    mu_d, mu_e = calc_path_efficiencies(
        inp.mu_har, inp.mu_mppt, inp.mu_pcdu, inp.mu_bat
    )

    # Step 4 – Required array power
    P_sa = calc_required_array_power(
        inp.P_e, inp.T_e, mu_e, inp.P_d, inp.T_d, mu_d
    )

    # Step 5 – Array area
    A_sa = calc_array_area(P_sa, P_EOL)

    # Step 6 – Array mass
    mass = calc_array_mass(A_sa, P_sa, inp.m_density, inp.m_sp)

    # Step 7 – Battery energy
    E_bat_J  = calc_battery_energy(inp.P_e, inp.T_e, inp.DOD, inp.mu_bat)
    E_bat_Wh = E_bat_J / 3600.0

    return {
        # Intermediate
        "L_d":   L_d,
        "P_EOL": P_EOL,
        "mu_d":  mu_d,
        "mu_e":  mu_e,
        # Primary outputs
        "P_sa":     P_sa,
        "A_sa":     A_sa,
        "mass":     mass,
        "E_bat_J":  E_bat_J,
        "E_bat_Wh": E_bat_Wh,
    }


# =============================================================================
# 4. OUTPUT
# =============================================================================

def print_results(inp: SolarArrayInputs, res: dict) -> None:
    """Pretty-print the sizing results with units and source equations."""

    sep = "=" * 60
    sub = "-" * 60

    print(f"\n{sep}")
    print("  SPACECRAFT SOLAR ARRAY SIZING RESULTS")
    print(sep)

    print("\n[INPUTS]")
    print(f"  Solar constant             S_e       = {inp.S_e:.1f}  W/m²")
    print(f"  Cell efficiency            mu_cell   = {inp.mu_cell*100:.1f}  %")
    print(f"  Inherent degradation       I_d       = {inp.I_d:.3f}")
    print(f"  Sun incidence angle        theta     = {inp.theta_deg:.1f}  deg")
    print(f"  Annual radiation deg.      deg/yr    = {inp.deg_per_year*100:.2f} %/yr")
    print(f"  Mission lifetime           L         = {inp.lifetime_yr:.1f}  yr")
    print(f"  Eclipse power demand       P_e       = {inp.P_e:.1f}  W")
    print(f"  Daylight power demand      P_d       = {inp.P_d:.1f}  W")
    print(f"  Eclipse duration           T_e       = {inp.T_e/60:.1f}  min")
    print(f"  Daylight duration          T_d       = {inp.T_d/60:.1f}  min")
    print(f"  Harness efficiency         mu_har    = {inp.mu_har:.3f}")
    print(f"  MPPT efficiency            mu_mppt   = {inp.mu_mppt:.3f}")
    print(f"  PCDU efficiency            mu_pcdu   = {inp.mu_pcdu:.3f}")
    print(f"  Battery efficiency         mu_bat    = {inp.mu_bat:.3f}")
    print(f"  Depth of discharge         DOD       = {inp.DOD:.2f}")
    # print(f"  Cell areal density         rho_cell  = {inp.rho_cell:.2f} kg/m²")
    # print(f"  Yoke mass                  m_yoke    = {inp.m_yoke:.2f} kg")
    # print(f"  Deployment mechanism mass  m_deploy  = {inp.m_deploy:.2f} kg")

    print(f"\n{sub}")
    print("[INTERMEDIATE RESULTS]")
    print(f"  Lifetime degradation       L_d  = (1 - {inp.deg_per_year})^{inp.lifetime_yr}")
    print(f"                                  = {res['L_d']:.4f}")
    print()
    print(f"  EOL power density          P_EOL = S_e × mu_cell × I_d × cos(theta) × f_temp × L_d")
    print(f"                                   = {res['P_EOL']:.2f}  W/m²")
    print()
    print(f"  Daylight path efficiency   mu_d  = mu_har × mu_mppt × mu_pcdu")
    print(f"                                   = {res['mu_d']:.4f}")
    print()
    print(f"  Eclipse path efficiency    mu_e  = mu_har × mu_pcdu × mu_bat")
    print(f"                                   = {res['mu_e']:.4f}")

    print(f"\n{sub}")
    print("[PRIMARY OUTPUTS]")
    print()
    print(f"  Required EOL array power       P_sa_EOL  = (P_e·T_e)/(mu_e·T_d) + P_d/mu_d")
    print(f"                                   = {res['P_sa']:.2f}  W")
    print()
    print(f"  Solar array area           A_sa  = P_sa / P_EOL")
    print(f"                                   = {res['A_sa']:.4f}  m²")
    print()
    print(f"  Required BOL array power       P_sa_BOL  = P_sa_EOL * L_d")
    print(f"                                   = {res['P_sa']/res['L_d']:.2f}  W")
    print()
    print(f"  Solar array mass breakdown:")
    print(f"    TOTAL                          = {res['mass']['m_total']:.3f}  kg")
    print()
    print(f"  Battery energy             E_bat = (P_e·T_e) / (DOD·mu_bat)")
    print(f"                                   = {res['E_bat_J']:.1f}  J")
    print(f"                                   = {res['E_bat_Wh']:.2f}  Wh")
    print(sep)

# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    inp_tug = SolarArrayInputs()        # ← Edit inputs here or subclass
    inp_ms = SolarArrayInputsMS()

    active_inp = inp_ms

    res = run_sizing(active_inp)
    print_results(active_inp, res)