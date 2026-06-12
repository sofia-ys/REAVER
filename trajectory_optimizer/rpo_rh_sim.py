"""
REAVER RPO — Part 2: Recycling-Hub RPO Simulation
=================================================
After all five debris captures the mothership performs five sequential docking
cycles at the cooperative Recycling Hub (RH).  Each cycle has two phases:

  Phase A  TUG MEET  MS (alone) intercepts an arriving tug+debris within 500 m of
                     RH and clamps it with the robotic arm  -> DV_RPO_TUG_MEET
  Phase B  DOCK      MS carries the combined body (MS + extended arm + tug +
                     debris) in to the RH and hard-docks via a dedicated port on
                     TOP of the bus (the arm stays extended, holding the payload
                     out to the side)                       -> DV_RPO_RH_DOCK

RH is cooperative (known state, fixed docking axis) so there is no spin-
synchronisation phase — the cost driver in Phase B is the combined-body inertia
and the CoM offset from the nominal MS CoM, not tumble matching.
"""
import numpy as np

from config import MS_DRY, TUG_DRY
from rpo_config import (M_DEBRIS_WC, M_TUG, T_HOLD_RH, AOM_MS,
                        R_KOS2_MEET, R_KOS1_DOCK, R_KOS2_DOCK, ARM_LEN_EXTENDED,
                        DEBRIS_PANEL_SPAN_WC, KOS1_MARGIN)
from rpo_guidance import translational_phase, dwell_phase
from rpo_control import attitude_prop, _angular_momentum, build_combined_body
from rpo_dynamics import box_inertia
from rpo_budget import prop_from_dv, dv_from_prop, PhaseResult
from rpo_config import MS_DIMS

D2R = np.pi / 180.0

# Keep-out radii (shared with the figures; KOS1 from the cooperative target's
# panel span, KOS2 from arm + tug geometry).
R_KOS1 = (DEBRIS_PANEL_SPAN_WC / 2.0) * KOS1_MARGIN

# ── Phase geometry in the RH Hill frame ──────────────────────────────────────
# Phase B (DOCK) is referenced to the RH docking port at the origin.  The
# combined body is carried in from the 500 m hold to a 30 m standoff, then runs
# the same close-approach corridor as the other two RPOs (Fig 2b): 30 m ->
# KOS1_DOCK alignment hold -> KOS2_DOCK hard-dock (top-port contact).
R_CAPTURE     = np.array([0.0, -500.0,       0.0])   # tug+debris hold point — arm capture here
R_DOCK_CLOSE  = np.array([0.0, -30.0,        0.0])   # close-approach start, 30 m from RH port
R_DOCK_KOS1   = np.array([0.0, -R_KOS1_DOCK, 0.0])   # KOS1 hold (alignment + go/no-go)
R_DOCK_KOS2   = np.array([0.0, -R_KOS2_DOCK, 0.0])   # KOS2 — top-port hard-dock contact
R_PORT        = np.array([0.0,    0.0,       0.0])   # RH docking port (reference origin)

# Phase A (MEET) is a CLOSE cooperative approach expressed in the frame of the
# tug+debris CoM (Fig 2a). The MS comes in along -V-bar from a 30 m standoff,
# through the KOS1 hold (synchronisation + comm relay), to the KOS2 bare-arm
# capture range — the same close-approach corridor as the debris RPO (Fig 1),
# but for a cooperative target: no spin-sync detour and only a short (~30 min)
# go/no-go dwell. KOS2 here is the ARM LENGTH ONLY (no tug stowed on the MS arm).
R_MEET_CLOSE = np.array([0.0, -30.0,        0.0])   # close-approach start, 30 m from tug+debris
R_MEET_KOS1  = np.array([0.0, -R_KOS1,      0.0])   # KOS1 hold (panel clearance, synch + comm)
R_MEET_CAP   = np.array([0.0, -R_KOS2_MEET, 0.0])   # KOS2 — bare-arm capture range (arm length only)
R_MEET_START = R_MEET_CLOSE                         # back-compat alias for plotting code


def _com_offset_dv(m_chaser, arm_len, omega_settle_dps=0.5):
    """
    Attitude ΔV to slew/settle the off-balance combined body (MS + extended arm
    + tug + debris) about its shifted CoM.

    Model
    -----
    The maneuver is costed as the RCS angular impulse needed to drive the
    combined inertia to a peak rotation rate ``omega_settle`` about an arbitrary
    [1,1,1] axis:  H = |J·omega|,  prop = H/(arm·vex),  dV = vex·ln(m/(m-prop)).
    This is a single momentum-build proxy (spin-up); a full rest-to-rest slew
    would double it (spin-up + spin-down) — left single-sided as the budget value.

    Why these rates (the assumption being justified)
    ------------------------------------------------
    The peak slew rate is capped by the low-impact docking velocity limit
    (V_MAX = 0.05 m/s): a rotation whips the EXTENDED payload tip at v = omega·r,
    and that tip velocity must stay within the contact-velocity limit.  With the
    tug+debris held ~5-6 m out on the arm, omega ≈ V_MAX/r ≈ 0.5 deg/s — hence:
      - Phase B re-point  : 0.5 deg/s  (deliberate slew to aim the top port at RH)
      - Phase A arm clamp : 0.2 deg/s  (smaller — only arrests the tip-off
                                        disturbance from the low-velocity capture)
    Both are far below the 6 deg/s (1 rpm) tumble ceiling (REQ-MIS-06).
    """
    omega = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0) * omega_settle_dps * D2R
    J_comb, _r_cm, m_total = build_combined_body(arm_len)
    H = _angular_momentum(J_comb, omega)
    prop = attitude_prop(H)
    dv = dv_from_prop(prop, m_total)
    return dv, prop


def simulate_rh_cycle(ms_mass, tug_prop_remaining=0.0, arm_len=ARM_LEN_EXTENDED):
    """
    One RH docking cycle.

    Parameters
    ----------
    ms_mass            : MS mass at the start of this cycle [kg]
    tug_prop_remaining : leftover tug propellant shed with the tug this cycle [kg]
    arm_len            : arm extension while carrying the tug+debris to RH [m]
                         (held fully extended throughout the dock approach)

    Returns
    -------
    dict with per-phase results and the two cycle constants.
    """
    m_payload = M_TUG + M_DEBRIS_WC + tug_prop_remaining

    # ── Phase A — MEET (MS alone, cooperative close approach, Fig 2a) ────────
    # Manoeuvre sequence:
    #   A1 close approach  30 m standoff -> KOS1            (translation, V-bar)
    #   A2 KOS1 hold       synchronisation + comm relay      (stationkeeping)
    #   A3 capture run     KOS1 -> KOS2 (= arm length only)   (translation)
    #   A4 arm extend+dock bare arm reaches the tug, clamps   (attitude / CoM shift)
    # KOS-routing penalties use the same inner-sphere convention as the debris
    # RPO (clip-test against the sphere being approached).
    a_app1 = translational_phase("  A1 close approach (30 m -> KOS1)",
                                 R_MEET_CLOSE, R_MEET_KOS1, ms_mass, r_kos=R_KOS1)
    a_hold = dwell_phase("  A2 KOS1 hold (synch + comm)",
                         ms_mass, T_HOLD_RH)
    a_app2 = translational_phase("  A3 capture run (KOS1 -> KOS2)",
                                 R_MEET_KOS1, R_MEET_CAP, ms_mass, r_kos=R_KOS2_MEET)
    # A4 arm extension + dock: the bare arm reaches out to the tug and clamps,
    # producing a CoM shift that RCS nulls with a small settling rotation.
    dock_dv, dock_prop = _com_offset_dv(ms_mass, arm_len=ARM_LEN_EXTENDED,
                                        omega_settle_dps=0.2)
    a_dock = PhaseResult("  A4 arm extend + dock (CoM shift)", dock_dv, dock_prop,
                         duration=120.0)

    dv_meet = a_app1.dv + a_hold.dv + a_app2.dv + a_dock.dv
    prop_meet = a_app1.prop + a_hold.prop + a_app2.prop + a_dock.prop

    # ── Phase B — DOCK (combined body to RH, top-port hard-dock) ─────────────
    # Manoeuvre sequence (close-approach corridor, like the other two RPOs):
    #   B1 carry      500 m hold -> 30 m standoff               (translation)
    #   B2 re-point   slew the off-balance combined body so the top port lines
    #                 up with the RH port (dominant CoM-offset cost) (attitude)
    #   B3 approach   30 m -> KOS1_DOCK alignment hold           (translation)
    #   B4 align hold KOS1_DOCK go/no-go + final port alignment  (stationkeeping)
    #   B5 hard-dock  KOS1_DOCK -> KOS2_DOCK top-port contact    (translation)
    m_combined = ms_mass + m_payload
    b_carry = translational_phase("  B1 carry (hold -> 30 m standoff)",
                                  R_CAPTURE, R_DOCK_CLOSE, m_combined)
    # B2 attitude: re-point the large off-balance combined inertia (dominant driver).
    point_dv, point_prop = _com_offset_dv(m_combined, arm_len=arm_len, omega_settle_dps=0.5)
    b_point = PhaseResult("  B2 re-point (combined inertia)", point_dv, point_prop,
                          duration=120.0)
    b_app = translational_phase("  B3 close approach (30 m -> KOS1)",
                                R_DOCK_CLOSE, R_DOCK_KOS1, m_combined, r_kos=R_KOS1_DOCK)
    b_hold = dwell_phase("  B4 KOS1 hold (align + go/no-go)", m_combined, T_HOLD_RH)
    b_dock = translational_phase("  B5 hard-dock (KOS1 -> KOS2)",
                                 R_DOCK_KOS1, R_DOCK_KOS2, m_combined, r_kos=R_KOS2_DOCK)

    dv_dock = b_carry.dv + b_point.dv + b_app.dv + b_hold.dv + b_dock.dv
    prop_dock = b_carry.prop + b_point.prop + b_app.prop + b_hold.prop + b_dock.prop

    phases = [a_app1, a_hold, a_app2, a_dock,
              b_carry, b_point, b_app, b_hold, b_dock]
    duration = sum(p.duration for p in phases)

    return {
        'ms_mass': ms_mass,
        'phases': phases,
        'dv_meet': dv_meet,
        'dv_dock': dv_dock,
        'dv_rh': dv_meet + dv_dock,
        'prop_meet': prop_meet,
        'prop_dock': prop_dock,
        'duration': duration,
    }


def simulate_all_rh_cycles(ms_mass_floor=MS_DRY, n_cycles=5):
    """
    Run all five RH cycles.  The MS mass is heaviest on cycle 1 (just back from
    the return leg, still carrying the means to handle the heaviest payload) and
    lightest on cycle 5.  Without the optimizer's backward-pass mass timeline we
    bracket it with a simple decreasing schedule from floor+payload to floor;
    rpo_run.py can pass the true per-cycle masses when wired into the optimizer.
    """
    # Simple decreasing mass schedule (heaviest first), bracketed for budgeting.
    masses = np.linspace(ms_mass_floor + (TUG_DRY * (n_cycles - 1)),
                         ms_mass_floor, n_cycles)
    cycles = [simulate_rh_cycle(float(m)) for m in masses]
    return cycles
