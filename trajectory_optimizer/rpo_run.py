"""
REAVER RPO - Simulation Entry Point
===================================
Runs Part 1 (debris RPO) and Part 2 (RH docking cycles).

Prints three sections:
  CALC  -- step-by-step calculation trace for every phase (dV and time)
  RESULTS -- phase-by-phase summary table and per-leg propellant
  CONSTANTS -- five DV_RPO_* values ready to paste into config.py

    cd c:\\Projects\\DSE\\REAVER
    python trajectory_optimizer/rpo_run.py
"""
import numpy as np

from config import MS_DRY, MS_ISP, MS_VEX, G0
from rpo_config import (
    OMEGA_REQ_DPS, TUMBLE_SWEEP_DPS,
    N_MEAN, V_MAX, R_KOS2, R_KOS2_MEET, R_RETREAT,
    R_KOS1_DOCK, R_KOS2_DOCK, DOCK_PORT_HEIGHT,
    R_CUTOFF, INSPECT_STANDOFF, H_BAR_STANDOFF,
    DEBRIS_PANEL_SPAN_WC, KOS1_MARGIN,
    F_RCS, N_RCS_TOTAL, N_RCS_PER_EDGE, RCS_ARM, RCS_THR, MS_DIMS, TUG_DIMS,
    DEBRIS_BUS_DIMS, DEBRIS_PANEL_SPAN_WC, DEBRIS_PANEL_FRAC,
    M_DEBRIS_WC, M_TUG, ARM_LEN_EXTENDED,
    P_SOL, C_R, AOM_MS, A_TRIAX,
    T_DWELL_ANALYSIS, T_HOLD_RH, T_ARM_EXTEND, ABORT_MARGIN,
)
from rpo_dynamics import box_inertia, debris_inertia
from rpo_control import build_combined_body, build_ms_plus_tug, _angular_momentum
from rpo_debris_sim import simulate_debris_rpo, sweep_tumble_rates
from rpo_rh_sim import (simulate_all_rh_cycles,
                        R_MEET_CLOSE, R_MEET_KOS1, R_MEET_CAP,
                        R_CAPTURE, R_DOCK_CLOSE, R_DOCK_KOS1, R_DOCK_KOS2)

D2R = np.pi / 180.0

# COMSATBW-1 (catalogue index 17) is the reference target: confirmed panel span
# (17.2 m), part of the optimal mission sequence, and sizing driver for KOS1.
WORST_CASE_SEQ = [17, 19, 10, 1, 0]   # COMSATBW-1, METEOSAT10, GALAXY14, AMC-3, ECHOSTAR1


def debris_chaser_masses(seq=WORST_CASE_SEQ):
    """Real MS mass at each debris arrival from the optimizer forward pass."""
    import nominal_mission as nm
    r = nm.evaluate_sequence(list(seq))
    return np.asarray(r['mass_vec'][:5], float)


def _rule(char="-", n=72):
    print(char * n)


def _sec(title):
    print(f"\n{'=' * 72}")
    print(f"  {title}")
    print(f"{'=' * 72}")


# =============================================================================
# CALCULATION TRACE  — shows every intermediate step
# =============================================================================

def report_calculations(m0, omega_dps=OMEGA_REQ_DPS):
    """
    Print a step-by-step calculation trace for every phase.
    Two calculation types:
      TRANSLATION  -- CW two-impulse dV, time from chord / (0.5 * V_MAX)
      ATTITUDE     -- angular momentum H = J.omega, prop = H/(arm*vex)
    """
    _sec("CALCULATION TRACE  (how every dV and time is computed)")

    panel_span = DEBRIS_PANEL_SPAN_WC
    r_kos1 = (panel_span / 2.0) * KOS1_MARGIN
    omega_rad = omega_dps * D2R
    omega_vec = np.array([1., 1., 1.]) / np.sqrt(3.) * omega_rad

    print(f"""
  REFERENCE VALUES
  ----------------
  Chaser mass m0          = {m0:.0f} kg  (MS at first debris, from optimizer mass timeline)
  Exhaust velocity vex    = MS_ISP * g0 = {MS_ISP:.0f} * {G0:.4f} = {MS_VEX:.0f} m/s
  Mean motion n           = sqrt(mu / a^3) = {N_MEAN:.4e} rad/s  (GEO, a = {42878:.0f} km)
  Max approach velocity   = {V_MAX:.01f} m/s  (low-impact docking requirement)
  Panel span (COMSATBW-1) = {panel_span:.1f} m  (confirmed, BSS-702)
  KOS1 margin             = {KOS1_MARGIN:.2f}  (+20% safety on half-span)
  KOS1                    = ({panel_span:.1f}/2) x {KOS1_MARGIN:.2f} = {r_kos1:.2f} m
  KOS2                    = ARM ({ARM_LEN_EXTENDED:.1f} m) + tug/2 ({TUG_DIMS[2]/2:.2f} m) = {R_KOS2:.2f} m
  RCS arm (min)           = {RCS_ARM:.2f} m  (confirmed, worst-case moment arm)
  RCS thrusters           = {N_RCS_TOTAL} total  ({N_RCS_PER_EDGE} per active edge, 4 top edges + 4 bottom edges)
  Force per thruster      = {F_RCS:.1f} N  (confirmed)
  Force per edge          = {N_RCS_PER_EDGE} x {F_RCS:.1f} = {RCS_THR:.1f} N  (= RCS_THR, one side of torque couple)
  Torque (1 edge pair)    = 2 x {RCS_THR:.1f} x {RCS_ARM:.2f} = {2*RCS_THR*RCS_ARM:.2f} N.m
  Tumble rate             = {omega_dps:.1f} deg/s = {omega_rad:.4f} rad/s
""")

    # ── TRANSLATION helper ────────────────────────────────────────────────────
    def trace_translation(label, r0, rf, note=""):
        r0 = np.asarray(r0, float)
        rf = np.asarray(rf, float)
        chord = np.linalg.norm(rf - r0)
        t = max(chord / (0.5 * V_MAX), 1.0)

        # CW two-impulse
        from rpo_dynamics import cw_two_impulse
        dv1, dv2, dv = cw_two_impulse(r0, np.zeros(3), rf, np.zeros(3), N_MEAN, t)
        prop = float(m0 * (1. - np.exp(-dv / MS_VEX)))

        print(f"  {label}")
        if note:
            print(f"    Note    : {note}")
        print(f"    r0      = {r0}  m  (Hill frame: R-bar, V-bar, H-bar)")
        print(f"    rf      = {rf}  m")
        print(f"    Chord   = |rf - r0| = {chord:.2f} m")
        print(f"    Time    = chord / (0.5 * V_MAX)")
        print(f"            = {chord:.2f} / (0.5 * {V_MAX:.1f}) = {t:.1f} s = {t/60:.1f} min")
        print(f"    CW STM  : solve Phi_rv * v0+ = rf - Phi_rr * r0  (2-impulse targeting)")
        print(f"    dv1     = {np.round(dv1,4)} m/s   |dv1| = {np.linalg.norm(dv1):.4f} m/s")
        print(f"    dv2     = {np.round(dv2,4)} m/s   |dv2| = {np.linalg.norm(dv2):.4f} m/s")
        print(f"    dV      = |dv1| + |dv2| = {dv:.4f} m/s")
        print(f"    Prop    = m0 * (1 - exp(-dV/vex))")
        print(f"            = {m0:.0f} * (1 - exp(-{dv:.4f}/{MS_VEX:.0f})) = {prop:.4f} kg")
        print()

    # ── STATIONKEEPING helper ─────────────────────────────────────────────────
    def trace_stationkeeping(label, dwell_s):
        a_srp = C_R * AOM_MS * P_SOL
        a_dist = a_srp + A_TRIAX
        dv = a_dist * dwell_s
        prop = float(m0 * (1. - np.exp(-dv / MS_VEX)))

        print(f"  {label}")
        print(f"    Dwell time        = {dwell_s:.0f} s = {dwell_s/3600:.1f} h")
        print(f"    SRP accel         = C_R * (A/m) * P_sol")
        print(f"                      = {C_R} * {AOM_MS} * {P_SOL:.2e} = {a_srp:.2e} m/s^2")
        print(f"    Triaxiality accel = {A_TRIAX:.2e} m/s^2  (J22, GEO tangential)")
        print(f"    Total disturbance = {a_dist:.2e} m/s^2")
        print(f"    dV = a_dist * t   = {a_dist:.2e} * {dwell_s:.0f} = {dv:.4f} m/s")
        print(f"    Prop              = {m0:.0f} * (1 - exp(-{dv:.4f}/{MS_VEX:.0f})) = {prop:.4f} kg")
        print()

    # ── ATTITUDE helper ───────────────────────────────────────────────────────
    def trace_attitude(label, J, mass_body, note=""):
        H = _angular_momentum(J, omega_vec)
        prop = H / (RCS_ARM * MS_VEX)
        dv = float(MS_VEX * np.log(mass_body / (mass_body - prop))) if prop < mass_body else 0.
        t = H / (2. * RCS_THR * RCS_ARM)
        eig = np.linalg.eigvalsh(J)

        print(f"  {label}")
        if note:
            print(f"    Note       : {note}")
        print(f"    Inertia J  = diag({np.round(eig,0)}) kg.m^2  (principal moments)")
        print(f"    omega      = {omega_dps:.1f} deg/s = {omega_rad:.4f} rad/s")
        print(f"                 direction = [1,1,1]/sqrt(3) (arbitrary tumble axis)")
        print(f"    H = |J.omega|= {H:.2f} N.m.s")
        print(f"    Prop = H / (arm * vex)")
        print(f"         = {H:.2f} / ({RCS_ARM:.2f} * {MS_VEX:.0f}) = {prop:.4f} kg")
        print(f"    dV   = vex * ln(m / (m - prop))")
        print(f"         = {MS_VEX:.0f} * ln({mass_body:.0f} / ({mass_body:.0f} - {prop:.4f})) = {dv:.4f} m/s")
        print(f"    Time = H / (2 * F_edge * arm)  [1 opposed edge pair]")
        print(f"         = {H:.2f} / (2 * {RCS_THR:.1f} * {RCS_ARM:.2f}) = {t:.2f} s = {t/60:.2f} min")
        print()

    _rule()
    print("  PART 1 -- DEBRIS RPO  (phase by phase)")
    _rule()

    # P1 — Approach to the V-bar inspection hold
    r_cutoff = np.array([0., -R_CUTOFF, 0.])
    r_vhold  = np.array([0., -INSPECT_STANDOFF, 0.])
    r_hbar   = np.array([0., 0., H_BAR_STANDOFF])
    r_kos1v  = np.array([0., 0., r_kos1])
    trace_translation("P1   Approach  (-V-bar cut-off -> V-bar inspection hold)",
                      r_cutoff, r_vhold,
                      note="V-bar is the only fixed CW hold point; KOS1 keep-out active")

    # P2 — Inspection dwell (stable V-bar hold)
    trace_stationkeeping("P2   Inspection dwell  (V-bar hold: ID spin rate + tumble axis, go/no-go)",
                         T_DWELL_ANALYSIS)

    # P3 — Transition to the identified tumble axis
    trace_translation("P3a  Reposition  (V-bar hold -> tumble-axis / +H-bar standoff)",
                      r_vhold, r_hbar,
                      note="reposition around the keep-out sphere onto the measured tumble axis")
    trace_translation("P3b  Approach  (tumble-axis standoff -> KOS1)",
                      r_hbar, r_kos1v)

    # P4 — Arm extension (deploy the arm + tug at KOS1, before spin-up)
    J_ready, _, m_ready = build_ms_plus_tug()
    J_stowed = box_inertia(MS_DRY, MS_DIMS)
    omega_vec_n = omega_vec
    H_ready   = _angular_momentum(J_ready,   omega_vec_n)
    H_stowed  = _angular_momentum(J_stowed,  omega_vec_n)
    H_delta   = abs(H_ready - H_stowed)
    prop_p4   = H_delta / (RCS_ARM * MS_VEX)
    dv_p4     = float(MS_VEX * np.log(MS_DRY / (MS_DRY - prop_p4)))
    t_p4      = H_delta / (2. * RCS_THR * RCS_ARM)
    print(f"  P4   Arm extension  (deploy arm + tug at KOS1, before spin-up)")
    print(f"    Note    : debris NOT yet coupled; arm carries tug ({TUG_DIMS[0]:.2f}x{TUG_DIMS[1]:.2f}x{TUG_DIMS[2]:.2f} m) only")
    print(f"    J_stowed= diag({np.round(np.linalg.eigvalsh(J_stowed),0)}) kg.m^2  (MS alone)")
    print(f"    J_ready = diag({np.round(np.linalg.eigvalsh(J_ready),0)}) kg.m^2  (MS + tug on arm)")
    print(f"    H_stowed= {H_stowed:.2f} N.m.s")
    print(f"    H_ready = {H_ready:.2f} N.m.s")
    print(f"    dH      = |H_ready - H_stowed| = {H_delta:.2f} N.m.s  (extra momentum from inertia growth)")
    print(f"    Prop    = dH / (arm * vex) = {H_delta:.2f} / ({RCS_ARM:.2f} * {MS_VEX:.0f}) = {prop_p4:.4f} kg")
    print(f"    dV      = {dv_p4:.4f} m/s")
    print(f"    RCS time= dH / (2 * F_edge * arm) = {t_p4:.2f} s  (momentum impulse only)")
    print(f"    Time    = {T_ARM_EXTEND:.0f} s = {T_ARM_EXTEND/60:.2f} min  (mechanical boom deploy + grapple-lock; sets the phase duration)")
    print()

    # P5 — Spin sync (spin the deployed assembly to match the tumble)
    J_ms = box_inertia(MS_DRY, MS_DIMS)
    trace_attitude("P5   Spin sync  (spin the deployed assembly to match the tumble)",
                   J_ms, MS_DRY,
                   note="spin-up momentum + the P4 inertia growth sum to spinning the "
                        "deployed MS+tug to the tumble rate (order-independent)")

    # P6a — Final approach
    trace_translation("P6a  Final approach  (KOS1 -> KOS2, within approach cone)",
                      r_kos1v, np.array([0., 0., R_KOS2]),
                      note="5 deg terminal approach-cone corridor active")

    # P6b — Detumble
    J_comb, _, m_comb = build_combined_body()
    trace_attitude("P6b  Momentum dumping  (MS + arm + tug + debris, DETUMBLE)",
                   J_comb, m_comb,
                   note=f"combined body: MS {MS_DRY:.0f} kg + arm + tug {TUG_DIMS[0]*M_DEBRIS_WC**0:.0f} + "
                        f"debris {M_DEBRIS_WC:.0f} kg  -- AOCS dimensioning case (REQ-ACS-M4)")

    # P7 — Retreat
    trace_translation("P7   MS retreat  (KOS2 -> safe standoff after tug release)",
                      np.array([0., 0., R_KOS2]),
                      np.array([0., 0., R_RETREAT]),
                      note="step 4c of RPO Logistics sequence diagram (MD-MS-LRR)")

    print(f"  ABORT MARGIN  (+{ABORT_MARGIN*100:.0f}% on nominal dV sum, risk TR-MT-08)")
    print(f"    Nominal dV  = P1 + P2 + P3 + P4 + P5 + P6a + P7  (not P6b detumble)")
    print(f"    DV_RPO_DEBRIS = nominal_dV * (1 + {ABORT_MARGIN:.2f})")
    print(f"    DV_RPO_DETUMBLE = P6b dV  (extracted separately, no abort margin)")
    print()

    # =========================================================================
    _rule()
    print("  PART 2 -- RH DOCKING, PHASE A: TUG MEET  (cooperative close approach)")
    _rule()
    print(f"""
  PHASE A REFERENCE VALUES
  ------------------------
  Target            = cooperative tug+debris (known state, tug telemetry relay)
  Approach axis     = -V-bar (in-plane), frame centred on the tug+debris CoM
  Standoff start    = 30 m
  KOS1              = {r_kos1:.2f} m  (panel clearance, same sizing as the debris RPO)
  KOS2 (MEET)       = arm length only = {R_KOS2_MEET:.2f} m
                      bare MS arm reaches the tug; NO tug stowed on the arm, so
                      tug/2 is NOT added (contrast debris KOS2 = {R_KOS2:.2f} m)
  KOS1 hold         = {T_HOLD_RH:.0f} s = {T_HOLD_RH/3600:.1f} h  (synchronisation + comm relay, go/no-go)
""")

    # A1 — Close approach (30 m standoff -> KOS1), V-bar corridor
    trace_translation("A1   Close approach  (30 m standoff -> KOS1)",
                      R_MEET_CLOSE, R_MEET_KOS1,
                      note="cooperative target: direct -V-bar corridor, no H-bar detour")

    # A2 — KOS1 hold (synchronisation + comm relay)
    trace_stationkeeping("A2   KOS1 hold  (synchronisation + comm relay, go/no-go)",
                         T_HOLD_RH)

    # A3 — Capture run (KOS1 -> KOS2 = arm length only)
    trace_translation("A3   Capture run  (KOS1 -> KOS2 = arm length, no tug on arm)",
                      R_MEET_KOS1, R_MEET_CAP,
                      note="KOS2 = arm only; bare MS arm reaches to the cooperative tug")

    # A4 — Arm extension + dock (CoM shift, attitude impulse)
    omega_dock_dps = 0.2
    omega_dock = np.array([1., 1., 1.]) / np.sqrt(3.) * omega_dock_dps * D2R
    J_dock, _r_cm, m_dock = build_combined_body(ARM_LEN_EXTENDED)
    H_dock    = _angular_momentum(J_dock, omega_dock)
    prop_dock = H_dock / (RCS_ARM * MS_VEX)
    dv_dock   = float(MS_VEX * np.log(m_dock / (m_dock - prop_dock)))
    t_dock    = H_dock / (2. * RCS_THR * RCS_ARM)
    eig_d     = np.linalg.eigvalsh(J_dock)
    print("  A4   Arm extension + dock  (bare arm reaches the tug, clamps -> CoM shift)")
    print(f"    Note       : MS grabs the cooperative tug+debris at arm's length ({ARM_LEN_EXTENDED:.1f} m)")
    print(f"    Combined J = diag({np.round(eig_d,0)}) kg.m^2  (MS + arm + tug + debris)")
    print(f"    Combined m = {m_dock:.0f} kg  (MS {MS_DRY:.0f} + tug {M_TUG:.0f} + debris {M_DEBRIS_WC:.0f})")
    print(f"    settle w   = {omega_dock_dps:.1f} deg/s = {omega_dock_dps*D2R:.4f} rad/s  (small re-settle after clamp)")
    print(f"    H = |J.w|  = {H_dock:.2f} N.m.s")
    print(f"    Prop = H / (arm * vex)")
    print(f"         = {H_dock:.2f} / ({RCS_ARM:.2f} * {MS_VEX:.0f}) = {prop_dock:.4f} kg")
    print(f"    dV   = vex * ln(m / (m - prop))")
    print(f"         = {MS_VEX:.0f} * ln({m_dock:.0f} / ({m_dock:.0f} - {prop_dock:.4f})) = {dv_dock:.4f} m/s")
    print(f"    Time = H / (2 * F_edge * arm) = {H_dock:.2f} / (2 * {RCS_THR:.1f} * {RCS_ARM:.2f}) = {t_dock:.2f} s")
    print()
    print(f"  DV_RPO_TUG_MEET = A1 + A2 + A3 + A4  (cooperative target, no abort margin)")
    print()

    # =========================================================================
    _rule()
    print("  PART 2 -- RH DOCKING, PHASE B: DOCK  (top-port hard-dock with RH)")
    _rule()
    print(f"""
  PHASE B REFERENCE VALUES
  ------------------------
  Chaser            = MS + EXTENDED arm + tug + debris (off-balance combined body)
  Docking interface = dedicated port on TOP of the MS bus ({DOCK_PORT_HEIGHT:.2f} m tall);
                      the robotic arm stays extended, holding the payload to the side
  Approach axis     = -V-bar (in-plane), frame centred on the RH docking port
  Carry start       = 500 m hold  ->  30 m standoff
  KOS1 (DOCK)       = {R_KOS1_DOCK:.2f} m  (alignment hold; clears the extended
                      arm+tug+debris envelope = (arm {ARM_LEN_EXTENDED:.1f} + tug {TUG_DIMS[2]:.1f}
                      + bus {DEBRIS_BUS_DIMS[2]:.1f}) x {KOS1_MARGIN:.2f} margin)
  KOS2 (DOCK)       = {R_KOS2_DOCK:.2f} m  (top-port contact = MS height/2 {MS_DIMS[2]/2:.2f}
                      + port {DOCK_PORT_HEIGHT:.2f})
  KOS1 hold         = {T_HOLD_RH:.0f} s = {T_HOLD_RH/3600:.1f} h  (alignment + go/no-go)
""")

    # B1 — Carry (500 m hold -> 30 m standoff)
    trace_translation("B1   Carry  (500 m hold -> 30 m standoff)",
                      R_CAPTURE, R_DOCK_CLOSE,
                      note="combined body translates in toward the RH along -V-bar")

    # B2 — Re-point the off-balance combined body (dominant attitude driver)
    omega_pt_dps = 0.5
    omega_pt = np.array([1., 1., 1.]) / np.sqrt(3.) * omega_pt_dps * D2R
    J_pt, _r_cm_b, m_pt = build_combined_body(ARM_LEN_EXTENDED)
    H_pt    = _angular_momentum(J_pt, omega_pt)
    prop_pt = H_pt / (RCS_ARM * MS_VEX)
    dv_pt   = float(MS_VEX * np.log(m_pt / (m_pt - prop_pt)))
    t_pt    = H_pt / (2. * RCS_THR * RCS_ARM)
    eig_pt  = np.linalg.eigvalsh(J_pt)
    print("  B2   Re-point  (slew the off-balance body so the top port faces the RH)")
    print(f"    Note       : heavy debris held out on the {ARM_LEN_EXTENDED:.1f} m arm -> CoM far off the port axis")
    print(f"    Combined J = diag({np.round(eig_pt,0)}) kg.m^2  (MS + arm + tug + debris)")
    print(f"    Combined m = {m_pt:.0f} kg")
    print(f"    slew w     = {omega_pt_dps:.1f} deg/s = {omega_pt_dps*D2R:.4f} rad/s")
    print(f"    H = |J.w|  = {H_pt:.2f} N.m.s")
    print(f"    Prop = H / (arm * vex)")
    print(f"         = {H_pt:.2f} / ({RCS_ARM:.2f} * {MS_VEX:.0f}) = {prop_pt:.4f} kg")
    print(f"    dV   = vex * ln(m / (m - prop)) = {dv_pt:.4f} m/s")
    print(f"    Time = H / (2 * F_edge * arm) = {t_pt:.2f} s")
    print()

    # B3 — Close approach (30 m -> KOS1_DOCK)
    trace_translation("B3   Close approach  (30 m standoff -> KOS1)",
                      R_DOCK_CLOSE, R_DOCK_KOS1,
                      note="combined body; KOS1 keeps the extended payload clear of the RH")

    # B4 — KOS1 alignment hold
    trace_stationkeeping("B4   KOS1 hold  (final port alignment + go/no-go)", T_HOLD_RH)

    # B5 — Hard-dock (KOS1_DOCK -> KOS2_DOCK, top-port contact)
    trace_translation("B5   Hard-dock  (KOS1 -> KOS2, top-port contact)",
                      R_DOCK_KOS1, R_DOCK_KOS2,
                      note="MS top port (0.40 m) makes contact with the RH port at KOS2")

    print(f"  DV_RPO_RH_DOCK = B1 + B2 + B3 + B4 + B5  (no abort margin)")
    print()

    _rule()
    print("  PROPELLANT FORMULA  (Tsiolkovsky, used for every phase)")
    _rule()
    print(f"    m_prop = m0 * (1 - exp(-dV / vex))")
    print(f"    where  vex = ISP * g0 = {MS_ISP:.0f} * {G0:.4f} = {MS_VEX:.0f} m/s")
    print(f"    m0 = chaser wet mass at the START of each phase")
    print(f"         (decreases through the mission as prop is burned and tugs released)")
    print()


# =============================================================================
# RESULTS SUMMARY
# =============================================================================

def report_debris():
    _sec("PART 1 RESULTS - DEBRIS RPO  (COMSATBW-1 reference, 2440 kg)")

    masses = debris_chaser_masses()
    m_worst = float(masses[0])
    res = simulate_debris_rpo(omega_dps=OMEGA_REQ_DPS, chaser_mass=m_worst)

    print(f"\n  Tumble rate : {OMEGA_REQ_DPS:.1f} deg/s (1 rpm, REQ-MIS-06 ceiling)")
    print(f"  Chaser mass : {m_worst:.0f} kg  (first debris visit, from optimizer mass timeline)")
    print(f"  KOS1        : {res['r_kos1']:.2f} m  ({res['panel_span']:.1f}/2 x {KOS1_MARGIN:.2f} margin)")
    print(f"  KOS2        : {R_KOS2:.2f} m  (arm {ARM_LEN_EXTENDED:.1f} m + tug/2 {TUG_DIMS[2]/2:.2f} m)\n")
    _rule()
    print(f"  {'Phase':<36} {'dV [m/s]':>9} {'Prop [kg]':>10} {'Time':>9}")
    _rule()
    for p in res['phases']:
        print(f"  {p}")
    _rule()
    print(f"  {'DV_RPO_DEBRIS  (nominal x1.50 abort)':<36} {res['dv_debris']:9.3f} {res['prop_debris']:10.4f}")
    print(f"  {'DV_RPO_DETUMBLE  (no abort margin)':<36} {res['dv_detumble']:9.3f} {res['prop_detumble']:10.4f}")
    print(f"  {'Total wall-clock duration':<36} {'':>9} {'':>10} {res['duration']/3600:>8.2f} h")
    _rule()

    print("\n  Per-leg propellant  (real optimizer mass timeline):")
    print(f"  {'Leg':<6} {'Chaser [kg]':>12} {'dV [m/s]':>10} {'Prop [kg]':>10}")
    _rule()
    tot_prop = 0.0
    for i, m in enumerate(masses, 1):
        leg = simulate_debris_rpo(omega_dps=OMEGA_REQ_DPS, chaser_mass=float(m))
        prop = leg['prop_debris'] + leg['prop_detumble']
        tot_prop += prop
        print(f"  {i:<6} {m:12.0f} {leg['dv_debris']+leg['dv_detumble']:10.3f} {prop:10.3f}")
    _rule()
    print(f"  {'Total':>6} {'':>12} {'':>10} {tot_prop:10.3f} kg")

    print("\n  Tumble-rate sweep  (Conings & Mooij Table 5 benchmark):")
    _rule()
    print(f"  {'omega [d/s]':>11} {'DV_DEBRIS [m/s]':>16} {'DV_DETUMBLE [m/s]':>18} {'Prop [kg]':>10}")
    _rule()
    for r in sweep_tumble_rates(TUMBLE_SWEEP_DPS):
        print(f"  {r['omega_dps']:11.1f} {r['dv_debris']:16.3f} {r['dv_detumble']:18.3f} "
              f"{r['prop_debris']+r['prop_detumble']:10.3f}")
    return res


def report_rh():
    _sec("PART 2 RESULTS - RECYCLING-HUB RPO  (MS x5 docking cycles)")
    cycles = simulate_all_rh_cycles(ms_mass_floor=MS_DRY, n_cycles=5)
    r_kos1 = (DEBRIS_PANEL_SPAN_WC / 2.0) * KOS1_MARGIN

    # ── Phase-by-phase breakdown for the reference (heaviest) cycle ──────────
    # Same table format as Part 1: Phase / dV / Prop / Time, with the TUG-MEET
    # and DOCK subtotals (Phase A = A1..A4, Phase B = B1..B5).
    ref = cycles[0]
    print(f"\n  Reference cycle : cycle 1 (heaviest MS in the schedule)")
    print(f"  MS mass         : {ref['ms_mass']:.0f} kg")
    print(f"  Approach axis   : -V-bar (cooperative target, in-plane)")
    print(f"  TUG-MEET KOS1   : {r_kos1:.2f} m  (panel clearance, {DEBRIS_PANEL_SPAN_WC:.1f}/2 x {KOS1_MARGIN:.2f})")
    print(f"  TUG-MEET KOS2   : {R_KOS2_MEET:.2f} m  (bare arm length, no tug on the arm)")
    _rule()
    print(f"  DOCK KOS1       : {R_KOS1_DOCK:.2f} m  (extended-payload envelope, alignment hold)")
    print(f"  DOCK KOS2       : {R_KOS2_DOCK:.2f} m  (top-port contact = MS height/2 + {DOCK_PORT_HEIGHT:.2f} m port)\n")
    _rule()
    print(f"  {'Phase':<36} {'dV [m/s]':>9} {'Prop [kg]':>10} {'Time':>9}")
    _rule()
    for p in ref['phases']:
        print(f"  {p}")
    _rule()
    print(f"  {'DV_RPO_TUG_MEET (A1+A2+A3+A4)':<36} {ref['dv_meet']:9.3f} {ref['prop_meet']:10.4f}")
    print(f"  {'DV_RPO_RH_DOCK  (B1..B5)':<36} {ref['dv_dock']:9.3f} {ref['prop_dock']:10.4f}")
    print(f"  {'DV_RPO_RH       (MEET + DOCK)':<36} {ref['dv_rh']:9.3f} "
          f"{ref['prop_meet']+ref['prop_dock']:10.4f}")
    print(f"  {'Total wall-clock duration':<36} {'':>9} {'':>10} {ref['duration']/3600:>8.2f} h")
    _rule()

    # ── Per-cycle summary across the 5-cycle decreasing mass schedule ────────
    print("\n  Per-cycle summary  (decreasing MS mass schedule):")
    _rule()
    print(f"  {'Cycle':>5} {'MS [kg]':>8} {'MEET [m/s]':>11} {'DOCK [m/s]':>11} "
          f"{'RH [m/s]':>9} {'Prop [kg]':>10}")
    _rule()
    tot_meet = tot_dock = tot_prop = 0.
    for i, c in enumerate(cycles, 1):
        prop = c['prop_meet'] + c['prop_dock']
        tot_meet += c['dv_meet']; tot_dock += c['dv_dock']; tot_prop += prop
        print(f"  {i:>5} {c['ms_mass']:8.0f} {c['dv_meet']:11.3f} {c['dv_dock']:11.3f} "
              f"{c['dv_rh']:9.3f} {prop:10.4f}")
    _rule()
    n = len(cycles)
    meet_avg, dock_avg = tot_meet / n, tot_dock / n
    print(f"  {'avg':>5} {'':8} {meet_avg:11.3f} {dock_avg:11.3f} "
          f"{meet_avg+dock_avg:9.3f} {tot_prop:10.4f}")
    _rule()
    return {'meet_avg': meet_avg, 'dock_avg': dock_avg,
            'rh_avg': meet_avg + dock_avg, 'total_prop': tot_prop}


def emit_constants(debris, rh):
    _sec("PART 3 - RPO dV CONSTANTS  (paste into config.py)")
    dv_d = debris['dv_debris']
    dv_t = debris['dv_detumble']
    dv_m = rh['meet_avg']
    dv_k = rh['dock_avg']
    dv_rh = dv_m + dv_k

    print(f"  DV_RPO_DEBRIS   = {dv_d:6.2f}   # m/s")
    print(f"  DV_RPO_DETUMBLE = {dv_t:6.2f}   # m/s")
    print(f"  DV_RPO_TUG_MEET = {dv_m:6.2f}   # m/s")
    print(f"  DV_RPO_RH_DOCK  = {dv_k:6.2f}   # m/s")
    print(f"  DV_RPO_RH       = {dv_rh:6.2f}   # m/s  (= TUG_MEET + DOCK)")

    debris_total = dv_d + dv_t
    total_rpo_dv = 5 * debris_total + 5 * dv_rh
    m_prop = MS_DRY * (np.exp(total_rpo_dv / MS_VEX) - 1.)
    # Sanity: the combined-body momentum dump (DETUMBLE) is the AOCS-dimensioning
    # case and is expected to be the single largest line with the short RCS arm.
    ok = dv_t >= dv_d and dv_k >= dv_m
    print(f"\n  Consistency: DETUMBLE ({dv_t:.2f}) >= DEBRIS approach ({dv_d:.2f}) and "
          f"DOCK ({dv_k:.2f}) >= MEET ({dv_m:.2f})  -> {'PASS' if ok else 'FAIL'}")
    print(f"  Total RPO dV : 5x{debris_total:.1f} + 5x{dv_rh:.1f} = {total_rpo_dv:.1f} m/s")
    print(f"  Prop upper bound @ MS_DRY {MS_DRY:.0f} kg, Isp {MS_ISP:.0f} s : {m_prop:.1f} kg")
    _rule("=")


def main():
    masses = debris_chaser_masses()
    m_worst = float(masses[0])

    report_calculations(m_worst)
    debris = report_debris()
    rh     = report_rh()
    emit_constants(debris, rh)


if __name__ == "__main__":
    main()
