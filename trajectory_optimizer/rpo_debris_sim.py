"""
REAVER RPO — Part 1: Debris RPO Simulation
==========================================
Close-range rendezvous with the worst-case tumbling debris target.  An
uncooperative target gives no state and no fixed docking axis, so the MS must
first observe it from a stable hold and only then approach along the tumble
axis it has measured.  Phase map:

  Phase 1  Approach ........ cut-off (150 m) -> V-bar inspection hold (20 m)
  Phase 2  Inspection ...... dwell on V-bar: ID spin rate + tumble axis, go/no-go
  Phase 3  Transition ...... reposition to the tumble axis (H-bar) -> KOS1
  Phase 4  Arm Extension ... deploy the arm + tug at KOS1 (before spin-up)
  Phase 5  Spin Sync ....... spin the deployed assembly to match the tumble
  Phase 6  MEV Lock + Detumble . final approach (KOS1->KOS2) + momentum dump
  Phase 7  Retreat ......... back off to a safe standoff after tug release

Why V-bar then H-bar: V-bar is the only naturally fixed CW hold point (an
along-track offset at rest stays put), so inspection is done there; the capture
corridor is along the tumble axis (modelled as H-bar), entered only after the
axis is known.

Outputs the two debris-RPO constants:
    DV_RPO_DEBRIS   = approach + inspection-SK + transition + spin-sync
                      + arm-ext + final-app + retreat (+ abort margin)
    DV_RPO_DETUMBLE = momentum dumping of the combined body (extracted separately)
"""
import numpy as np

from config import MS_DRY
from rpo_config import (R_CUTOFF, INSPECT_STANDOFF, H_BAR_STANDOFF, R_KOS2, R_RETREAT,
                        OMEGA_REQ_DPS, TUMBLE_SWEEP_DPS, T_DWELL_ANALYSIS, T_ARM_EXTEND,
                        ABORT_MARGIN, DEBRIS_PANEL_SPAN_WC, KOS1_MARGIN)
from rpo_guidance import translational_phase, dwell_phase
from rpo_control import spin_sync_phase, arm_extension_phase, detumble_phase
from rpo_budget import PhaseResult


def simulate_debris_rpo(omega_dps=OMEGA_REQ_DPS, chaser_mass=MS_DRY,
                        abort_margin=ABORT_MARGIN,
                        panel_span=DEBRIS_PANEL_SPAN_WC):
    """
    Run the five-phase debris RPO for one tumble rate.

    Parameters
    ----------
    panel_span : float
        Solar-panel tip-to-tip span of the TARGET debris [m].
        KOS1 = panel_span / 2 is computed here, mimicking the operational
        flow: GO uplinks the confirmed span to the MS at the cut-off
        handover (RF -> LiDAR), and the MS sets KOS1 accordingly.
        Defaults to the worst-case target (EUTE 12 WEST A, HS-601, 26.2 m).

    Returns
    -------
    dict with keys:
        phases          : list[PhaseResult]   (per-phase breakdown)
        dv_debris       : float   DV_RPO_DEBRIS [m/s]  (incl. abort margin)
        dv_detumble     : float   DV_RPO_DETUMBLE [m/s]
        prop_debris     : float   propellant for DV_RPO_DEBRIS phases [kg]
        prop_detumble   : float   propellant for de-tumble [kg]
        duration        : float   total phase wall-clock time [s]
        r_kos1          : float   KOS1 used [m]  (for reporting)
        panel_span      : float   panel span used [m]  (for reporting)
    """
    # KOS1 derived from GO-uplinked panel span at cut-off distance,
    # plus a 20% safety margin (KOS1_MARGIN) on the half-span.
    # COMSATBW-1: 17.2/2 * 1.20 = 10.3 m.
    r_kos1_val = (panel_span / 2.0) * KOS1_MARGIN

    m = chaser_mass

    # ── Phase geometry ──────────────────────────────────────────────────────
    # Hold + inspection are on V-bar (the only naturally fixed CW hold point);
    # the capture corridor is along the IDENTIFIED tumble axis (modelled as
    # +H-bar).  The MS only transitions to the tumble axis after the V-bar
    # inspection has measured the spin rate + axis and ground gives go/no-go.
    r_cutoff = np.array([0.0, -R_CUTOFF, 0.0])           # -V-bar sensor handover (150 m)
    r_vhold  = np.array([0.0, -INSPECT_STANDOFF, 0.0])   # V-bar inspection hold (KOS1-clear)
    r_hbar   = np.array([0.0, 0.0, H_BAR_STANDOFF])      # tumble-axis (H-bar) standoff
    r_kos1   = np.array([0.0, 0.0, r_kos1_val])          # KOS1 on the tumble axis
    r_kos2   = np.array([0.0, 0.0, R_KOS2])

    # ── Phase 1 — Approach to the V-bar inspection hold ─────────────────────
    p1 = translational_phase("P1 Approach (-Vbar cutoff -> Vbar hold)",
                             r_cutoff, r_vhold, m, r_kos=r_kos1_val)

    # ── Phase 2 — Inspection at the V-bar hold (ID spin rate + tumble axis) ──
    # Stable fixed V-bar hold: characterise the target and obtain go/no-go
    # BEFORE committing to the tumble-axis approach.
    p2 = dwell_phase("P2 Inspection (Vbar hold: spin+axis, go/no-go)",
                     m, T_DWELL_ANALYSIS)

    # ── Phase 3 — Transition to the identified tumble axis ──────────────────
    # Reposition from the V-bar hold around the keep-out sphere onto the tumble
    # axis, then close to KOS1.
    p3a = translational_phase("  P3 reposition (Vbar hold -> tumble-axis standoff)",
                              r_vhold, r_hbar, m, r_kos=r_kos1_val)
    p3b = translational_phase("  P3 approach (standoff -> KOS1)", r_hbar, r_kos1, m)
    p3 = PhaseResult("P3 Transition to tumble axis", p3a.dv + p3b.dv,
                     p3a.prop + p3b.prop, p3a.duration + p3b.duration)

    # ── Phase 4 — Arm extension (deploy the arm + tug at KOS1, before spin-up)
    # Mechanically limited: the boom deploy + grapple-lock onto the tug takes
    # ~1.5 min (T_ARM_EXTEND), far longer than the few-second RCS impulse.
    p4 = arm_extension_phase("P4 Arm extension", m0=m, omega_dps=omega_dps,
                             duration=T_ARM_EXTEND)

    # ── Phase 5 — Spin synchronisation (spin the deployed assembly to match) ─
    # Done AFTER extension: the spin-up momentum + the arm's inertia growth sum
    # to spinning the deployed MS+tug to the tumble rate (order-independent).
    p5 = spin_sync_phase("P5 Spin sync (match tumble)", omega_dps=omega_dps, m0=m)

    # ── Phase 6 — MEV lock (final approach) + momentum dumping ──────────────
    p6_app = translational_phase("  P6 final approach (KOS1 -> KOS2)",
                                 r_kos1, r_kos2, m, r_kos=R_KOS2)
    p6_det = detumble_phase("P6 Momentum dumping (DETUMBLE)", omega_dps=omega_dps)

    # ── Phase 7 — MS retreat after tug release (sequence-diagram step 4c) ────
    # MS backs from the capture position out to a safe standoff before the
    # long-range departure burn (which the optimizer already accounts for).
    p7 = translational_phase("P7 MS retreat (KOS2 -> safe)", r_kos2,
                             np.array([0.0, 0.0, R_RETREAT]), m)

    phases = [p1, p2, p3, p4, p5, p6_app, p6_det, p7]

    # ── Aggregate ───────────────────────────────────────────────────────────
    nominal_dv = p1.dv + p2.dv + p3.dv + p4.dv + p5.dv + p6_app.dv + p7.dv
    nominal_prop = p1.prop + p2.prop + p3.prop + p4.prop + p5.prop + p6_app.prop + p7.prop

    dv_debris = nominal_dv * (1.0 + abort_margin)
    prop_debris = nominal_prop * (1.0 + abort_margin)

    duration = sum(p.duration for p in phases)

    return {
        'omega_dps': omega_dps,
        'panel_span': panel_span,
        'r_kos1': r_kos1_val,
        'phases': phases,
        'dv_debris': dv_debris,
        'dv_detumble': p6_det.dv,
        'prop_debris': prop_debris,
        'prop_detumble': p6_det.prop,
        'duration': duration,
    }


def sweep_tumble_rates(rates=TUMBLE_SWEEP_DPS, chaser_mass=MS_DRY):
    """Re-run the debris RPO across the tumble-rate sweep (Conings Table 5)."""
    return [simulate_debris_rpo(omega_dps=r, chaser_mass=chaser_mass) for r in rates]
