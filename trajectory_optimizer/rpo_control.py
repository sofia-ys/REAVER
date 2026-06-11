"""
REAVER RPO Control — attitude / rotational phases
=================================================
Costs the rotational manoeuvres of the RPO sequence in propellant and
equivalent ΔV:

  * Spin synchronisation  — MS (tug stowed, near-nominal inertia) spins up to
                            match the tumbling debris (Conings Eq. 19-25,
                            modelled here at the angular-impulse / budget level).
  * Momentum dumping      — combined MS + arm + tug + debris body de-tumbled to
                            zero rotation w.r.t. Earth (the AOCS dimensioning
                            case, REQ-ACS-M4 / risk TR-19) -> DV_RPO_DETUMBLE.

Budget model
------------
Reaction wheels carry the *smooth* tracking torque but cannot hold a net spin
indefinitely; the net angular-momentum change is therefore delivered (or
dumped) by the RCS thrusters.  For a momentum change ``H`` with thrusters at
effective moment arm ``L`` and exhaust velocity ``vex``:

    torque   tau   = 2 * F * L           (opposed thruster pair)
    burn time t    = H / tau
    propellant m_p = (2 F / vex) * t = H / (L * vex)

i.e. attitude propellant depends only on the angular momentum and the moment
arm, not on the thrust level.  The equivalent ΔV follows from Tsiolkovsky.
This is intentionally conservative (no reaction-wheel storage credit), which is
the correct posture for the worst-case AOCS sizing the checklist asks for.
"""
import numpy as np

from config import MS_VEX, MS_DRY
from rpo_config import (RCS_ARM, F_RCS, N_RCS_PER_EDGE, OMEGA_REQ_DPS, MS_DIMS, TUG_DIMS,
                        DEBRIS_BUS_DIMS, DEBRIS_PANEL_SPAN_WC, DEBRIS_PANEL_FRAC,
                        M_DEBRIS_WC, M_TUG, ARM_LEN_EXTENDED)
from rpo_dynamics import box_inertia, debris_inertia, combined_inertia
from rpo_budget import dv_from_prop, PhaseResult

D2R = np.pi / 180.0


def _angular_momentum(J, omega_vec):
    """Magnitude of the angular momentum H = J omega for a body."""
    return float(np.linalg.norm(J @ omega_vec))


def attitude_prop(H, arm=RCS_ARM, vex=MS_VEX):
    """Propellant [kg] to deliver/dump angular momentum H via RCS thrusters."""
    return float(H / (arm * vex))


def spin_sync_phase(name, omega_dps=OMEGA_REQ_DPS, m0=MS_DRY, arm=RCS_ARM):
    """
    Spin-synchronisation cost.  Per Rev 2 the tug is still stowed on the arm
    here, so the MS inertia is near nominal (arm length ~ 0) — this is the
    favourable result the RPO-Logistics clarification produces.
    """
    omega = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0) * omega_dps * D2R  # arbitrary axis
    J_ms = box_inertia(m0, MS_DIMS)
    H = _angular_momentum(J_ms, omega)
    prop = attitude_prop(H, arm)
    dv = dv_from_prop(prop, m0)
    # Conings: P-rate controller down to the 0.6 deg/s threshold, then INDI.
    # Budget duration: time to build the spin at the RCS torque (order of mins).
    return PhaseResult(name, dv, prop, duration=_spin_time(H, arm),
                       note="tug stowed, MS inertia near-nominal")


def _spin_time(H, arm, force_per_edge=None):
    """
    Burn time to build/dump angular momentum H using one opposed edge pair.
    tau = 2 * (N_per_edge * F_single) * arm  = 2 * RCS_THR * arm
    """
    if force_per_edge is None:
        force_per_edge = float(N_RCS_PER_EDGE * F_RCS)  # = 12.0 N
    tau = 2.0 * force_per_edge * arm
    return float(H / tau)


def build_ms_plus_tug(arm_len=ARM_LEN_EXTENDED):
    """
    MS + tug-on-arm inertia (NO debris).  This is the configuration during
    Phase 4 arm extension: the arm carries only the tug toward the LAE nozzle;
    the debris has not yet coupled to the MS.

    Returns (J, r_cm, m_total).
    """
    J_ms = box_inertia(MS_DRY, MS_DIMS)
    J_tug = box_inertia(M_TUG, TUG_DIMS)
    bodies = [
        (MS_DRY, J_ms, [0.0, 0.0, 0.0]),
        (M_TUG, J_tug, [0.0, 0.0, arm_len]),
    ]
    return combined_inertia(bodies)


def build_combined_body(arm_len=ARM_LEN_EXTENDED):
    """
    Worst-case combined inertia tensor for momentum dumping: MS at origin, arm
    extended along +z, tug at the arm tip, debris just beyond the tug.  This is
    the post-MEV-lock configuration (Phase 5) — the AOCS dimensioning case.

    Returns (J_combined, r_cm, m_total).
    """
    J_ms = box_inertia(MS_DRY, MS_DIMS)
    J_tug = box_inertia(M_TUG, TUG_DIMS)
    J_deb = debris_inertia(M_DEBRIS_WC, DEBRIS_BUS_DIMS,
                           DEBRIS_PANEL_SPAN_WC, DEBRIS_PANEL_FRAC)

    half_deb = DEBRIS_BUS_DIMS[2] / 2.0
    bodies = [
        (MS_DRY, J_ms, [0.0, 0.0, 0.0]),
        (M_TUG, J_tug, [0.0, 0.0, arm_len]),
        (M_DEBRIS_WC, J_deb, [0.0, 0.0, arm_len + half_deb + 0.5]),
    ]
    return combined_inertia(bodies)


def detumble_phase(name, omega_dps=OMEGA_REQ_DPS, arm=RCS_ARM, arm_len=ARM_LEN_EXTENDED):
    """
    Momentum-dumping (de-tumble) cost for the combined body at full arm
    extension — the dimensioning AOCS case -> DV_RPO_DETUMBLE.
    """
    omega = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0) * omega_dps * D2R
    J_comb, _r_cm, m_total = build_combined_body(arm_len)
    H = _angular_momentum(J_comb, omega)
    prop = attitude_prop(H, arm)
    dv = dv_from_prop(prop, m_total)
    return PhaseResult(name, dv, prop, duration=_spin_time(H, arm),
                       note=f"combined body {m_total:.0f} kg, arm extended {arm_len:.1f} m")


def arm_extension_phase(name, m0, omega_dps=OMEGA_REQ_DPS, arm=RCS_ARM,
                        arm_len=ARM_LEN_EXTENDED, duration=None):
    """
    Attitude-hold cost while the arm extends the tug from stowed to ready: the
    growing inertia must be continuously spun up to keep pace with the target
    rate.  Budgeted as the angular-momentum delta between the stowed (near-
    nominal MS) and ready (combined) configurations.

    The phase is mechanically (not momentum) limited: the RCS impulse takes only
    a few seconds, but the physical boom deploy + grapple-lock onto the tug sets
    the wall-clock time.  Pass ``duration`` to use that deploy time instead of
    the (much shorter) momentum-build time.
    """
    omega = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0) * omega_dps * D2R
    J_stowed = box_inertia(m0, MS_DIMS)
    J_ready, _r_cm, _m = build_ms_plus_tug(arm_len)   # tug only, NOT debris
    H_delta = abs(_angular_momentum(J_ready, omega) - _angular_momentum(J_stowed, omega))
    prop = attitude_prop(H_delta, arm)
    dv = dv_from_prop(prop, m0)
    t = _spin_time(H_delta, arm) if duration is None else float(duration)
    return PhaseResult(name, dv, prop, duration=t,
                       note="tug-on-arm inertia growth (debris not yet coupled)")
