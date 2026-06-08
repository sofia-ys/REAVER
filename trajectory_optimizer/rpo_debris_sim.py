"""
REAVER RPO — Part 1: Debris RPO Simulation
==========================================
Close-range rendezvous with the worst-case tumbling debris target
(EUTE 12 WEST A, 2700 kg, 1 rpm).  Implements the Rev 2 REAVER phase map:

  Phase 1  Approach .............. cut-off (150 m) -> H-bar (20 m) -> KOS1 (5 m)
  Phase 2  Debris Analysis ....... dwell at interaction distance (no net ΔV)
  Phase 3  Tug Deploy + Spin Sync. MS spins up to match target (tug stowed)
  Phase 4  Arm Extension ......... attitude hold while inertia grows
  Phase 5  MEV Lock + De-tumble .. final approach (KOS1->KOS2) + momentum dump

Outputs the two debris-RPO constants:
    DV_RPO_DEBRIS   = approach + analysis-SK + spin-sync + arm-ext + final-app
                      (+ abort margin)
    DV_RPO_DETUMBLE = momentum dumping of the combined body (extracted separately)
"""
import numpy as np

from config import MS_DRY
from rpo_config import (R_CUTOFF, H_BAR_STANDOFF, R_KOS1, R_KOS2,
                        OMEGA_REQ_DPS, TUMBLE_SWEEP_DPS, T_DWELL_ANALYSIS,
                        ABORT_MARGIN)
from rpo_guidance import translational_phase, dwell_phase
from rpo_control import spin_sync_phase, arm_extension_phase, detumble_phase
from rpo_budget import PhaseResult


def simulate_debris_rpo(omega_dps=OMEGA_REQ_DPS, chaser_mass=MS_DRY,
                        abort_margin=ABORT_MARGIN):
    """
    Run the five-phase debris RPO for one tumble rate.

    Returns
    -------
    dict with keys:
        phases          : list[PhaseResult]   (per-phase breakdown)
        dv_debris       : float   DV_RPO_DEBRIS [m/s]  (incl. abort margin)
        dv_detumble     : float   DV_RPO_DETUMBLE [m/s]
        prop_debris     : float   propellant for DV_RPO_DEBRIS phases [kg]
        prop_detumble   : float   propellant for de-tumble [kg]
        duration        : float   total phase wall-clock time [s]
    """
    m = chaser_mass

    # ── Phase 1 — Approach (Conings Ph.1+2) ─────────────────────────────────
    # -V-bar cut-off hold -> +H-bar standoff -> KOS1 (interaction distance).
    r_cutoff = np.array([0.0, -R_CUTOFF, 0.0])
    r_hbar   = np.array([0.0, 0.0, H_BAR_STANDOFF])
    r_kos1   = np.array([0.0, 0.0, R_KOS1])

    p1a = translational_phase("  P1 transfer (-Vbar -> +Hbar)", r_cutoff, r_hbar,
                              m, r_kos=R_KOS1)
    p1b = translational_phase("  P1 approach (Hbar -> KOS1)", r_hbar, r_kos1, m)
    p1 = PhaseResult("P1 Approach", p1a.dv + p1b.dv, p1a.prop + p1b.prop,
                     p1a.duration + p1b.duration)

    # ── Phase 2 — Debris analysis dwell (stationkeeping only) ───────────────
    p2 = dwell_phase("P2 Debris analysis (dwell)", m, T_DWELL_ANALYSIS)

    # ── Phase 3 — Tug deploy + spin synchronisation (tug stowed) ────────────
    p3 = spin_sync_phase("P3 Spin sync (tug stowed)", omega_dps=omega_dps, m0=m)

    # ── Phase 4 — Arm extension (attitude hold during inertia growth) ───────
    p4 = arm_extension_phase("P4 Arm extension", m0=m, omega_dps=omega_dps)

    # ── Phase 5 — MEV lock (final approach) + momentum dumping ──────────────
    p5_app = translational_phase("  P5 final approach (KOS1 -> KOS2)",
                                 r_kos1, np.array([0.0, 0.0, R_KOS2]), m)
    p5_det = detumble_phase("P5 Momentum dumping (DETUMBLE)", omega_dps=omega_dps)

    phases = [p1, p2, p3, p4, p5_app, p5_det]

    # ── Aggregate ───────────────────────────────────────────────────────────
    nominal_dv = p1.dv + p2.dv + p3.dv + p4.dv + p5_app.dv
    nominal_prop = p1.prop + p2.prop + p3.prop + p4.prop + p5_app.prop

    dv_debris = nominal_dv * (1.0 + abort_margin)
    prop_debris = nominal_prop * (1.0 + abort_margin)

    duration = sum(p.duration for p in phases)

    return {
        'omega_dps': omega_dps,
        'phases': phases,
        'dv_debris': dv_debris,
        'dv_detumble': p5_det.dv,
        'prop_debris': prop_debris,
        'prop_detumble': p5_det.prop,
        'duration': duration,
    }


def sweep_tumble_rates(rates=TUMBLE_SWEEP_DPS, chaser_mass=MS_DRY):
    """Re-run the debris RPO across the tumble-rate sweep (Conings Table 5)."""
    return [simulate_debris_rpo(omega_dps=r, chaser_mass=chaser_mass) for r in rates]
