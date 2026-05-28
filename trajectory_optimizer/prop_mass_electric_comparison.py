"""
Tug Propellant Mass Calculator
================================
Compares two methods for computing tug propellant mass:
  Method 1: Tsiolkovsky rocket equation (constant Isp assumption)
  Method 2: RK4 numerical integration (constant power, variable Isp)

The RK4 integrates the coupled ODE system:
  d(dv_acc)/dt = T(t) / m(t)          [acceleration]
  dm/dt        = -T(t) / vex(t)       [mass consumption]

With constant-power throttle model:
  vex(t) = vex0 * (m0/m)^alpha        [exhaust velocity varies with mass]
  T(t)   = 2 * eta * P / vex(t)       [thrust adjusted to maintain power]
  mdot(t)= T(t) / vex(t)             [mass flow rate]

alpha = 0   -> constant thrust  = rocket equation (baseline)
alpha = 0.1 -> mild throttle-up (realistic Hall thruster, default)
alpha = 0.5 -> significant throttle-up
alpha = 1.0 -> aggressive throttle-up (upper bound)

Author: REAVER Group 7
"""

import numpy as np
import sys

# =============================================================================
# CONSTANTS
# =============================================================================
MU  = 3.986004418e14   # Earth gravitational parameter [m^3/s^2]
G0  = 9.807            # Standard gravity [m/s^2]
DAY = 86400.0          # Seconds per day
D2R = np.pi / 180.0

# =============================================================================
# TUG & MISSION PARAMETERS  —  edit these
# =============================================================================
TUG_DRY   = 200.0      # Tug dry mass [kg]
TUG_ISP   = 1600.0     # Nominal Isp [s]
TUG_VEX   = TUG_ISP * G0   # Nominal exhaust velocity [m/s]
TUG_THR   = 0.05        # Nominal thrust [N]
ETA       = 0.60       # Thruster efficiency [-]

# Input power implied by nominal operating point
# P = T * vex / (2 * eta)
P_INPUT   = TUG_THR * TUG_VEX / (2 * ETA)  # [W]

# Throttle exponent for RK4 method
# 0.0 = constant thrust (= rocket eq), 0.1 = realistic Hall thruster
ALPHA     = 0.10

# RK4 time step
DT        = 0.5 *DAY        # 1 day in seconds

# Recycling Hub orbit
RH_SMA    = 42164.2e3  # [m]
RH_INC    = 7.0        # [deg]
RH_RAAN   = 10.0       # [deg]

# =============================================================================
# DEBRIS CATALOGUE  —  add or remove objects as needed
# =============================================================================
DEBRIS = [
    # name,                        mass_kg,   sma_m,         inc_deg,  raan_deg
    ('NSS 5 (Intelsat 803)',        2060.46,  42494.423e3,   10.4813,  48.1496),
    ('Sirius 2 (GE-1E)',            1762.14,  42410.084e3,   11.9727,  34.1101),
    ('EUTE W2',                     1793.86,  42451.133e3,   11.6171,  38.4600),
    ('Skynet 4E',                   1490.00,  42516.912e3,   12.0550,   9.2220),
    ('EUTE 16C (SESAT 1)',          2600.00,  42537.778e3,   10.8393,  46.2495),
    ('Skynet 4F',                   1489.00,  42482.488e3,   12.1030,  18.6123),
    ('EUTE 12 West A (AB 1)',       2700.00,  42728.199e3,    7.2278,  68.7938),
    ('EUTE 33A (EB 3)',             1552.00,  42558.620e3,    9.6938,  55.8840),
    ('SESAT 2 (Express AM-22)',     2542.00,  42420.904e3,    8.2239,  64.3413),
    ('Meteosat 9 (MSG 2)',          2054.00,  42163.863e3,    9.3136,  54.3880),
    ('EUTE 8A (Sinosat 3)',         2320.00,  42716.964e3,    9.5162,  48.3713),
    ('COMSATBW-1',                 2440.00,  42164.245e3,    0.0696,  89.6878),
    ('COMSATBW-2',                 2440.00,  42164.548e3,    0.0344, 111.5455),
    ('HYLAS 1',                    2570.00,  42164.346e3,    5.5899,  73.2864),
    ('Meteosat 10 (MSG 3)',         2035.00,  42164.226e3,    4.6216,  61.1229),
    ('Meteosat 11 (MSG 4)',         2043.00,  42164.271e3,    3.1471,  71.2262),
]

# =============================================================================
# ORBITAL MECHANICS
# =============================================================================

def edelbaum_dv(sma_a, inc_a, raan_a, sma_b, inc_b, raan_b):
    """
    Edelbaum low-thrust optimal transfer dV [m/s].
    Combines altitude change and plane change (inc + RAAN) simultaneously.
    """
    v1 = np.sqrt(MU / sma_a)
    v2 = np.sqrt(MU / sma_b)
    cos_dth = (np.cos(inc_a*D2R) * np.cos(inc_b*D2R) +
               np.sin(inc_a*D2R) * np.sin(inc_b*D2R) *
               np.cos((raan_b - raan_a) * D2R))
    dth = np.arccos(np.clip(cos_dth, -1.0, 1.0))   # combined plane angle [rad]
    dv  = np.sqrt(v1**2 + v2**2 - 2*v1*v2*np.cos(np.pi/2 * dth))
    return dv, np.degrees(dth)

def size_wet_mass(dv_target, m_payload, vex, thrust, alpha, p_input, eta,
                  tol=0.01, max_iter=80):
    """
    Iteratively size tug wet mass.
    m_wet = m_payload + m_prop, where m_prop depends on m_wet.
    Uses RK4 propellant estimate internally for consistency.
    """
    m_wet = m_payload * 1.35   # initial guess
    for _ in range(max_iter):
        _, mp, _ = rk4_integrate(dv_target, m_wet, vex, thrust,
                                  alpha, p_input, eta, verbose=False)
        m_new = m_payload + mp
        if abs(m_new - m_wet) < tol:
            m_wet = m_new
            break
        m_wet = 0.6 * m_wet + 0.4 * m_new
    return m_wet

# =============================================================================
# METHOD 1 — ROCKET EQUATION
# =============================================================================

def rocket_equation(dv, m_wet, vex):
    """
    Tsiolkovsky rocket equation.
    Assumes constant Isp (constant exhaust velocity).

    m_prop = m_wet * (1 - exp(-dv / vex))
    t_burn = (m_wet * vex / T) * (1 - exp(-dv / vex))
    """
    exp_term = np.exp(-dv / vex)
    m_prop   = m_wet * (1.0 - exp_term)
    t_burn   = (m_wet * vex / TUG_THR) * (1.0 - exp_term)
    return m_prop, t_burn

# =============================================================================
# METHOD 2 — RK4 NUMERICAL INTEGRATION
# =============================================================================

def derivatives(m, m0, vex0, p_input, eta, alpha):
    """
    Compute derivatives for the state [dv_acc, m].

    Parameters
    ----------
    m      : current mass [kg]
    m0     : initial wet mass [kg]
    vex0   : nominal exhaust velocity [m/s]
    p_input: input power [W]
    eta    : thruster efficiency [-]
    alpha  : throttle exponent [-]

    Returns
    -------
    ddv_dt : acceleration [m/s^2]
    dm_dt  : mass flow rate [kg/s] (negative)
    T_cur  : current thrust [N]
    vex_cur: current exhaust velocity [m/s]
    mdot   : current mass flow rate magnitude [kg/s]
    """
    vex_cur  = vex0 * (m0 / m) ** alpha        # throttled exhaust velocity
    T_cur    = 2.0 * eta * p_input / vex_cur    # constant-power thrust
    mdot     = T_cur / vex_cur                  # mass flow rate
    ddv_dt   = T_cur / m                        # acceleration
    dm_dt    = -mdot                            # mass rate of change
    return ddv_dt, dm_dt, T_cur, vex_cur, mdot


def rk4_integrate(dv_target, m_wet, vex0, thrust, alpha, p_input, eta,
                  dt=DT, verbose=True):
    """
    RK4 numerical integration of the constant-power EP equations of motion.

    Integrates until accumulated dV reaches dv_target.
    Uses linear interpolation for the final fractional step.

    Parameters
    ----------
    dv_target : target delta-V [m/s]
    m_wet     : initial wet mass [kg]
    vex0      : nominal exhaust velocity [m/s]
    thrust    : nominal thrust [N] (used only for display)
    alpha     : throttle exponent [-]
    p_input   : input power [W]
    eta       : thruster efficiency [-]
    dt        : time step [s]
    verbose   : print step-by-step table

    Returns
    -------
    t_final  : transfer time [s]
    mp_final : propellant consumed [kg]
    steps    : list of step data dicts
    """
    # ── Initial state ──────────────────────────────────────────────────────
    t      = 0.0
    dv_acc = 0.0
    m      = m_wet
    m0     = m_wet
    steps  = []

    if verbose:
        print()
        print(f"  RK4 Integration  |  alpha={alpha}  |  dt={dt/3600:.0f}h  "
              f"|  P={p_input:.1f}W  |  Target dV={dv_target:.4f} m/s")
        print()
        print(f"  {'Step':>5}  {'Time (d)':>9}  {'dv_acc (m/s)':>13}  "
              f"{'Mass (kg)':>10}  {'vex (m/s)':>10}  {'T (N)':>7}  "
              f"{'mdot (kg/s)':>12}  {'Prop used (kg)':>15}  "
              f"{'dV remain (m/s)':>16}")
        print("  " + "="*110)

        # print step 0
        a0, _, T0, vex0_disp, md0 = derivatives(m, m0, vex0,
                                                  p_input, eta, alpha)
        print(f"  {0:>5}  {0.0:>9.3f}  {dv_acc:>13.6f}  {m:>10.4f}  "
              f"{vex0_disp:>10.2f}  {T0:>7.4f}  {md0:>12.8f}  "
              f"{0.0:>15.6f}  {dv_target:>16.6f}")

    step = 0
    t_final = None
    mp_final = None

    while True:
        step += 1

        # ── RK4 slopes ────────────────────────────────────────────────────
        # k1: slopes at (t, dv_acc, m)
        k1_a, k1_m, T1, vex1, md1 = derivatives(m, m0, vex0,
                                                  p_input, eta, alpha)

        # k2: slopes at midpoint using k1 estimate
        m2   = m      + 0.5 * dt * k1_m
        k2_a, k2_m, T2, vex2, md2 = derivatives(m2, m0, vex0,
                                                  p_input, eta, alpha)

        # k3: slopes at midpoint using k2 estimate
        m3   = m      + 0.5 * dt * k2_m
        k3_a, k3_m, T3, vex3, md3 = derivatives(m3, m0, vex0,
                                                  p_input, eta, alpha)

        # k4: slopes at end of step using k3 estimate
        m4   = m      + dt * k3_m
        k4_a, k4_m, T4, vex4, md4 = derivatives(m4, m0, vex0,
                                                  p_input, eta, alpha)

        # ── Weighted update ───────────────────────────────────────────────
        dv_new = dv_acc + (dt / 6.0) * (k1_a + 2*k2_a + 2*k3_a + k4_a)
        m_new  = m      + (dt / 6.0) * (k1_m + 2*k2_m + 2*k3_m + k4_m)
        t_new  = t + dt

        prop_used = m0 - m_new

        # ── Display ───────────────────────────────────────────────────────
        if verbose:
            vex_cur = vex0 * (m0 / m) ** alpha
            T_cur   = 2.0 * eta * p_input / vex_cur
            md_cur  = T_cur / vex_cur
            print(f"  {step:>5}  {t_new/DAY:>9.3f}  {dv_new:>13.6f}  "
                  f"{m_new:>10.4f}  {vex_cur:>10.2f}  {T_cur:>7.4f}  "
                  f"{md_cur:>12.8f}  {prop_used:>15.6f}  "
                  f"{dv_target-dv_new:>16.6f}")

        steps.append({'t': t_new, 'dv_acc': dv_new, 'm': m_new,
                      'prop_used': prop_used})

        # ── Check if dV target crossed ────────────────────────────────────
        if dv_new >= dv_target:
            # Linear interpolation to find exact crossing
            frac     = (dv_target - dv_acc) / (dv_new - dv_acc)
            t_final  = t + frac * dt
            m_final  = m + frac * (m_new - m)
            mp_final = m0 - m_final

            if verbose:
                print()
                print(f"  --> dV target reached at t = {t_final/DAY:.4f} days"
                      f" (linear interp between step {step-1} and {step})")
                print(f"  --> Final mass        = {m_final:.4f} kg")
                print(f"  --> Propellant used   = {mp_final:.4f} kg")
            break

        # ── Safety: max steps ─────────────────────────────────────────────
        if step > 5000:
            print("  WARNING: max steps reached without convergence.")
            t_final  = t_new
            mp_final = prop_used
            break

        dv_acc = dv_new
        m      = m_new
        t      = t_new

    return t_final, mp_final, steps

# =============================================================================
# COMPARISON RUNNER
# =============================================================================

def compare_all(verbose_rk4=False):
    """
    Run both methods for all debris objects and print comparison table.
    Set verbose_rk4=True to print full RK4 step tables for every object.
    """
    sys.stdout.reconfigure(encoding='utf-8')

    print("=" * 80)
    print("  TUG PROPELLANT MASS — METHOD COMPARISON")
    print("  Rocket Equation  vs  RK4 Numerical Integration")
    print("=" * 80)
    print(f"  Tug dry mass  : {TUG_DRY} kg")
    print(f"  Nominal Isp   : {TUG_ISP} s")
    print(f"  Nominal thrust: {TUG_THR} N")
    print(f"  Input power   : {P_INPUT:.1f} W  ({P_INPUT/1000:.2f} kW)")
    print(f"  Efficiency    : {ETA}")
    print(f"  Throttle exp  : alpha = {ALPHA}  (0=const thrust, 0.1=Hall thruster)")
    print(f"  RK4 time step : {DT/3600:.0f} hours")
    print(f"  Destination   : RH orbit (sma={RH_SMA/1e3:.1f}km, "
          f"i={RH_INC}deg, RAAN={RH_RAAN}deg)")
    print()

    header = (f"  {'Name':<28}  {'dV':>7}  {'dTheta':>7}  {'mwet':>7}  "
              f"{'mp_RE':>7}  {'t_RE':>7}  "
              f"{'mp_RK4':>7}  {'t_RK4':>7}  "
              f"{'diff kg':>8}  {'diff %':>7}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    all_results = []

    for name, mass, sma, inc, raan in DEBRIS:

        # ── Edelbaum dV ───────────────────────────────────────────────────
        dv, dtheta = edelbaum_dv(sma, inc, raan, RH_SMA, RH_INC, RH_RAAN)

        # ── Size wet mass (using RK4 internally for consistency) ──────────
        m_payload = TUG_DRY + mass
        m_wet     = size_wet_mass(dv, m_payload, TUG_VEX, TUG_THR,
                                   ALPHA, P_INPUT, ETA)

        # ── Method 1: Rocket equation ─────────────────────────────────────
        mp_RE, t_RE = rocket_equation(dv, m_wet, TUG_VEX)

        # ── Method 2: RK4 ─────────────────────────────────────────────────
        if verbose_rk4:
            print(f"\n  {'='*70}")
            print(f"  {name}")
            print(f"  {'='*70}")

        t_RK4, mp_RK4, _ = rk4_integrate(dv, m_wet, TUG_VEX, TUG_THR,
                                           ALPHA, P_INPUT, ETA,
                                           dt=DT, verbose=verbose_rk4)

        # ── Comparison ────────────────────────────────────────────────────
        diff_kg  = mp_RK4 - mp_RE
        diff_pct = diff_kg / mp_RE * 100

        print(f"  {name:<28}  {dv:>7.1f}  {dtheta:>7.3f}  {m_wet:>7.1f}  "
              f"{mp_RE:>7.3f}  {t_RE/DAY:>7.2f}  "
              f"{mp_RK4:>7.3f}  {t_RK4/DAY:>7.2f}  "
              f"{diff_kg:>+8.4f}  {diff_pct:>+7.4f}%")

        all_results.append({
            'name': name, 'mass': mass, 'dv': dv, 'dtheta': dtheta,
            'm_wet': m_wet, 'mp_RE': mp_RE, 't_RE': t_RE,
            'mp_RK4': mp_RK4, 't_RK4': t_RK4,
            'diff_kg': diff_kg, 'diff_pct': diff_pct
        })

    # ── Summary statistics ────────────────────────────────────────────────
    diffs_pct = [r['diff_pct'] for r in all_results]
    diffs_kg  = [r['diff_kg']  for r in all_results]
    worst     = all_results[np.argmax(np.abs(diffs_pct))]
    best      = all_results[np.argmin(np.abs(diffs_pct))]

    print()
    print("  " + "-" * 80)
    print(f"  {'Column guide':}")
    print(f"    dV      : Edelbaum delta-V [m/s]")
    print(f"    dTheta  : combined plane separation angle [deg]")
    print(f"    mwet    : tug wet mass at capture [kg]")
    print(f"    mp_RE   : propellant — rocket equation [kg]")
    print(f"    t_RE    : transfer time — rocket equation [days]")
    print(f"    mp_RK4  : propellant — RK4 integration [kg]")
    print(f"    t_RK4   : transfer time — RK4 integration [days]")
    print(f"    diff    : mp_RK4 - mp_RE")
    print()
    print(f"  SUMMARY (alpha = {ALPHA}):")
    print(f"    Mean error  : {np.mean(diffs_pct):+.4f}%  "
          f"({np.mean(diffs_kg):+.4f} kg)")
    print(f"    Max error   : {max(diffs_pct, key=abs):+.4f}%  "
          f"-> {worst['name']}")
    print(f"    Min error   : {min(diffs_pct, key=abs):+.4f}%  "
          f"-> {best['name']}")
    print(f"    Std dev     : {np.std(diffs_pct):.5f}%")
    print()
    print(f"  INTERPRETATION:")
    print(f"    Negative diff% means RK4 uses LESS propellant than rocket eq.")
    print(f"    Rocket equation is conservative (safe for design use).")
    print(f"    At ~3% propellant fraction, error < 0.2% for alpha=0.1.")
    print()

    return all_results


def single_object(name, mass_kg, sma_m, inc_deg, raan_deg,
                  alpha=ALPHA, print_steps=True):
    """
    Run full comparison for a single debris object with step-by-step RK4 output.

    Example
    -------
    single_object('SESAT 2', 2542.0, 42420.904e3, 8.2239, 64.3413)
    """
    sys.stdout.reconfigure(encoding='utf-8')

    print("=" * 80)
    print(f"  SINGLE OBJECT ANALYSIS: {name}")
    print("=" * 80)

    dv, dtheta = edelbaum_dv(sma_m, inc_deg, raan_deg,
                              RH_SMA, RH_INC, RH_RAAN)
    m_payload  = TUG_DRY + mass_kg
    m_wet      = size_wet_mass(dv, m_payload, TUG_VEX, TUG_THR,
                                alpha, P_INPUT, ETA)

    print(f"  Debris mass       : {mass_kg:.2f} kg")
    print(f"  Tug payload mass  : {m_payload:.2f} kg  (dry {TUG_DRY}kg + debris)")
    print(f"  Tug wet mass      : {m_wet:.4f} kg")
    print(f"  Edelbaum dV       : {dv:.4f} m/s")
    print(f"  Plane angle       : {dtheta:.4f} deg")
    print(f"  Input power       : {P_INPUT:.2f} W")
    print(f"  Throttle exponent : alpha = {alpha}")

    # Method 1
    mp_RE, t_RE = rocket_equation(dv, m_wet, TUG_VEX)
    print()
    print(f"  METHOD 1 — Rocket equation:")
    print(f"    m_prop = {mp_RE:.4f} kg")
    print(f"    t      = {t_RE/DAY:.4f} days")

    # Method 2
    print()
    print(f"  METHOD 2 — RK4 numerical integration (alpha={alpha}):")
    t_RK4, mp_RK4, steps = rk4_integrate(dv, m_wet, TUG_VEX, TUG_THR,
                                           alpha, P_INPUT, ETA,
                                           dt=DT, verbose=print_steps)

    print()
    print("  " + "=" * 60)
    print(f"  COMPARISON SUMMARY:")
    print(f"  {'':30}  {'Rocket eq':>12}  {'RK4':>12}  {'Difference':>12}")
    print(f"  {'':30}  {'':->12}  {'':->12}  {'':->12}")
    print(f"  {'Propellant mass [kg]':<30}  {mp_RE:>12.4f}  {mp_RK4:>12.4f}  "
          f"{mp_RK4-mp_RE:>+12.4f}")
    print(f"  {'Transfer time [days]':<30}  {t_RE/DAY:>12.4f}  {t_RK4/DAY:>12.4f}  "
          f"{(t_RK4-t_RE)/DAY:>+12.4f}")
    print(f"  {'Error [%]':<30}  {'':>12}  {'':>12}  "
          f"{(mp_RK4-mp_RE)/mp_RE*100:>+12.4f}%")
    print()

    return mp_RE, mp_RK4, t_RE, t_RK4


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(
        description='Tug propellant: rocket equation vs RK4 integration')
    parser.add_argument('--mode', choices=['all', 'single', 'steps'],
                        default='all',
                        help=('all    : comparison table for all 16 debris\n'
                              'single : detailed output for one debris (SESAT 2)\n'
                              'steps  : all debris with full RK4 step tables'))
    parser.add_argument('--alpha', type=float, default=ALPHA,
                        help='Throttle exponent (0=const thrust, 0.1=Hall thruster)')
    parser.add_argument('--dt', type=float, default=1.0,
                        help='RK4 time step in days (default 1.0)')
    args = parser.parse_args()

    # Override globals if provided
    ALPHA = args.alpha
    DT    = args.dt * DAY

    if args.mode == 'all':
        compare_all(verbose_rk4=False)

    elif args.mode == 'single':
        # SESAT 2 — worst-case dV object
        single_object('SESAT 2 (Express AM-22)',
                      mass_kg=2542.00, sma_m=42420.904e3,
                      inc_deg=8.2239,  raan_deg=64.3413,
                      alpha=args.alpha, print_steps=True)

    elif args.mode == 'steps':
        compare_all(verbose_rk4=True)