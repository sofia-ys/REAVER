"""
REAVER RPO — Part 2: Recycling-Hub RPO Simulation
=================================================
After all five debris captures the mothership performs five sequential docking
cycles at the cooperative Recycling Hub (RH).  Each cycle has two phases:

  Phase A  MEET  MS (alone) intercepts an arriving tug+debris within 500 m of RH
                 and clamps it with the robotic arm        -> DV_RPO_RH_MEET
  Phase B  DOCK  MS carries the combined body (MS+tug+debris) to the RH docking
                 port and hands over                       -> DV_RPO_RH_DOCK

RH is cooperative (known state, fixed docking axis) so there is no spin-
synchronisation phase — the cost driver in Phase B is the combined-body inertia
and the CoM offset from the nominal MS CoM, not tumble matching.
"""
import numpy as np

from config import MS_DRY, TUG_DRY
from rpo_config import (M_DEBRIS_WC, M_TUG, T_HOLD_RH, AOM_MS)
from rpo_guidance import translational_phase, dwell_phase
from rpo_control import attitude_prop, _angular_momentum, build_combined_body
from rpo_dynamics import box_inertia
from rpo_budget import prop_from_dv, dv_from_prop, PhaseResult
from rpo_config import MS_DIMS

D2R = np.pi / 180.0

# Phase geometry in the RH Hill frame.
R_MEET_START = np.array([0.0, -500.0, 0.0])   # MS hold, 500 m on -V-bar
R_CAPTURE    = np.array([0.0, -15.0, 0.0])    # robotic-arm capture range
R_PORT       = np.array([0.0, 0.0, 0.0])      # RH docking port (reference origin)


def _com_offset_dv(m_chaser, arm_len, omega_settle_dps=0.5):
    """
    Attitude ΔV from the CoM offset when the combined body must be re-pointed
    for the docking approach.  A small settling rotation (``omega_settle``) of
    the large combined inertia is delivered by RCS — this is the dominant
    Phase-B cost driver.
    """
    omega = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0) * omega_settle_dps * D2R
    J_comb, _r_cm, m_total = build_combined_body(arm_len)
    H = _angular_momentum(J_comb, omega)
    prop = attitude_prop(H)
    dv = dv_from_prop(prop, m_total)
    return dv, prop


def simulate_rh_cycle(ms_mass, tug_prop_remaining=0.0, arm_len=2.0):
    """
    One RH docking cycle.

    Parameters
    ----------
    ms_mass            : MS mass at the start of this cycle [kg]
    tug_prop_remaining : leftover tug propellant shed with the tug this cycle [kg]
    arm_len            : arm extension during the docking carry [m]

    Returns
    -------
    dict with per-phase results and the two cycle constants.
    """
    m_payload = M_TUG + M_DEBRIS_WC + tug_prop_remaining

    # ── Phase A — MEET (MS alone, cooperative approach) ─────────────────────
    a_hold = dwell_phase("  A1 hold (await arrival)", ms_mass, T_HOLD_RH)
    a_app  = translational_phase("  A2 approach (500 m -> capture)",
                                 R_MEET_START, R_CAPTURE, ms_mass)
    # A3 robotic-arm clamp: small attitude correction for the new CoM.
    clamp_dv, clamp_prop = _com_offset_dv(ms_mass, arm_len=0.5, omega_settle_dps=0.2)
    a_clamp = PhaseResult("  A3 arm clamp (CoM shift)", clamp_dv, clamp_prop,
                          duration=60.0)

    dv_meet = a_hold.dv + a_app.dv + a_clamp.dv
    prop_meet = a_hold.prop + a_app.prop + a_clamp.prop

    # ── Phase B — DOCK (combined body to RH port) ───────────────────────────
    m_combined = ms_mass + m_payload
    b_carry = translational_phase("  B1 carry (capture -> port)",
                                  R_CAPTURE, R_PORT, m_combined)
    # B1 attitude: re-point the large combined inertia (dominant driver).
    point_dv, point_prop = _com_offset_dv(m_combined, arm_len=arm_len, omega_settle_dps=0.5)
    b_point = PhaseResult("  B1 re-point (combined inertia)", point_dv, point_prop,
                          duration=120.0)
    b_dock = translational_phase("  B2 hard-dock (fine approach)",
                                 R_PORT + np.array([0.0, -2.0, 0.0]), R_PORT,
                                 m_combined)

    dv_dock = b_carry.dv + b_point.dv + b_dock.dv
    prop_dock = b_carry.prop + b_point.prop + b_dock.prop

    phases = [a_hold, a_app, a_clamp, b_carry, b_point, b_dock]
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
