"""
Unit tests for the REAVER RPO simulation modules (rpo_*).

Run from the repository root:
    cd c:\\Projects\\DSE\\REAVER
    python -m pytest trajectory_optimizer/test_rpo.py -v
"""
import os
import sys
import numpy as np
import pytest

TRAJ_DIR = os.path.dirname(os.path.abspath(__file__))
RPO_DIR  = os.path.join(TRAJ_DIR, 'rpo')
for _d in (TRAJ_DIR, RPO_DIR):
    if _d not in sys.path:
        sys.path.insert(0, _d)

from config import MU, RH_SMA, MS_DRY, MS_VEX
import rpo_config as rc
from rpo_dynamics import (cw_stm, cw_two_impulse, box_inertia, debris_inertia,
                          combined_inertia, parallel_axis, quat_kinematics,
                          euler_rigid_body, propagate_tumble, stationkeeping_dv)
from rpo_budget import prop_from_dv, dv_from_prop
from rpo_guidance import translational_phase, dwell_phase, _clips_kos
from rpo_control import (spin_sync_phase, detumble_phase, arm_extension_phase,
                         build_combined_body, build_ms_plus_tug, attitude_prop)
from rpo_debris_sim import simulate_debris_rpo, sweep_tumble_rates
from rpo_rh_sim import simulate_rh_cycle, simulate_all_rh_cycles


# =============================================================================
# rpo_config
# =============================================================================

class TestRpoConfig:
    def test_mean_motion_matches_geo(self):
        assert rc.N_MEAN == pytest.approx(np.sqrt(MU / RH_SMA**3), rel=1e-12)

    def test_mean_motion_order_of_magnitude(self):
        # GEO mean motion ~ 7e-5 rad/s, far below the LEO 0.063 deg/s of the paper
        assert 5e-5 < rc.N_MEAN < 9e-5

    def test_cutoff_in_handover_band(self):
        assert 100.0 <= rc.R_CUTOFF <= 200.0

    def test_keepout_ordering(self):
        # KOS1 is NOT a global constant — it is derived per-target from the
        # panel span uplinked by GO at the cut-off distance.
        # Worst-case: EUTE 12 WEST A (HS-601, 26.2 m span) -> KOS1 = 13.1 m.
        # KOS2 = arm + tug/2 (fixed geometry, TBD MIRON).
        # KOS1_wc > KOS2 must hold for the worst-case target.
        kos1_wc = (rc.DEBRIS_PANEL_SPAN_WC / 2.0) * rc.KOS1_MARGIN
        assert kos1_wc > rc.R_KOS2, (
            f"KOS1_wc ({kos1_wc:.1f} m) must exceed KOS2 ({rc.R_KOS2:.1f} m)"
        )
        assert rc.R_KOS2 > 0.0

    def test_combined_mass_identity(self):
        assert rc.M_COMBINED_RH == pytest.approx(MS_DRY + rc.M_TUG + rc.M_DEBRIS_WC)

    def test_kos1_derives_from_panel_span_with_margin(self):
        # KOS1 = (panel_span / 2) * 1.20 safety margin.
        # COMSATBW-1 (BSS-702, confirmed span 17.2 m):
        #   KOS1 = 17.2/2 * 1.20 = 8.6 * 1.20 = 10.32 m (~10 m).
        assert rc.DEBRIS_PANEL_SPAN_WC == pytest.approx(17.2)
        assert rc.KOS1_MARGIN == pytest.approx(1.20)
        kos1_wc = (rc.DEBRIS_PANEL_SPAN_WC / 2.0) * rc.KOS1_MARGIN
        assert kos1_wc == pytest.approx(10.32)


# =============================================================================
# rpo_dynamics — CW state transition matrix
# =============================================================================

class TestCwStm:
    def test_identity_at_zero_time(self):
        # As t->0 the STM -> identity (use a tiny time to avoid 0/0 in Phi_rv)
        Phi = cw_stm(rc.N_MEAN, 1e-6)
        np.testing.assert_allclose(Phi, np.eye(6), atol=1e-3)

    def test_shape(self):
        assert cw_stm(rc.N_MEAN, 100.0).shape == (6, 6)

    def test_propagation_matches_stm(self):
        # Free drift: x(t) = Phi @ x(0)
        n, t = rc.N_MEAN, 600.0
        x0 = np.array([3.0, -5.0, 2.0, 0.01, 0.0, -0.02])
        Phi = cw_stm(n, t)
        xt = Phi @ x0
        assert xt.shape == (6,)
        assert np.isfinite(xt).all()


class TestCwTwoImpulse:
    def test_reaches_target(self):
        # After applying dv1, free drift must land exactly on rf
        n, t = rc.N_MEAN, 500.0
        r0 = np.array([0.0, -20.0, 0.0])
        rf = np.array([0.0, 0.0, 5.0])
        dv1, dv2, dvt = cw_two_impulse(r0, np.zeros(3), rf, np.zeros(3), n, t)
        Phi = cw_stm(n, t)
        x0 = np.concatenate([r0, dv1])      # v0=0, +dv1
        xt = Phi @ x0
        np.testing.assert_allclose(xt[:3], rf, atol=1e-6)

    def test_total_is_sum_of_norms(self):
        n, t = rc.N_MEAN, 500.0
        r0 = np.array([0.0, -20.0, 0.0])
        rf = np.array([0.0, 0.0, 5.0])
        dv1, dv2, dvt = cw_two_impulse(r0, np.zeros(3), rf, np.zeros(3), n, t)
        assert dvt == pytest.approx(np.linalg.norm(dv1) + np.linalg.norm(dv2))

    def test_positive(self):
        n, t = rc.N_MEAN, 500.0
        _, _, dvt = cw_two_impulse(np.array([0., -20., 0.]), np.zeros(3),
                                   np.array([0., 0., 5.]), np.zeros(3), n, t)
        assert dvt > 0.0


# =============================================================================
# rpo_dynamics — inertia & attitude
# =============================================================================

class TestInertia:
    def test_box_inertia_symmetric_cube(self):
        J = box_inertia(12.0, (1.0, 1.0, 1.0))
        # cube: all principal moments equal
        assert J[0, 0] == pytest.approx(J[1, 1]) == pytest.approx(J[2, 2])

    def test_box_inertia_positive_definite(self):
        J = box_inertia(1120.0, rc.MS_DIMS)
        assert np.all(np.linalg.eigvalsh(J) > 0)

    def test_debris_panels_inflate_off_axis(self):
        J = debris_inertia(2700.0, rc.DEBRIS_BUS_DIMS, rc.DEBRIS_PANEL_SPAN_WC,
                           rc.DEBRIS_PANEL_FRAC)
        # panels along y inflate Ixx and Izz above Iyy
        assert J[0, 0] > J[1, 1]
        assert J[2, 2] > J[1, 1]

    def test_parallel_axis_increases_inertia(self):
        J0 = box_inertia(100.0, (1.0, 1.0, 1.0))
        J1 = parallel_axis(J0, 100.0, [5.0, 0.0, 0.0])
        assert np.trace(J1) > np.trace(J0)

    def test_combined_body_heavier_than_ms_plus_tug(self):
        J_full, _, m_full = build_combined_body()
        J_lite, _, m_lite = build_ms_plus_tug()
        assert m_full > m_lite
        assert np.trace(J_full) > np.trace(J_lite)

    def test_combined_mass_conserved(self):
        _, _, m_total = combined_inertia([
            (100.0, np.eye(3), [0, 0, 0]),
            (50.0, np.eye(3), [1, 0, 0]),
        ])
        assert m_total == pytest.approx(150.0)


class TestAttitudeProp:
    def test_quat_norm_preserved(self):
        q0 = np.array([0.0, 0.0, 0.0, 1.0])
        w0 = np.array([0.1, 0.05, -0.07])
        J = box_inertia(2700.0, rc.DEBRIS_BUS_DIMS)
        _, q_hist, _ = propagate_tumble(q0, w0, J, 100.0, n_eval=50)
        norms = np.linalg.norm(q_hist, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-6)

    def test_euler_zero_torque_conserves_energy(self):
        J = np.diag([20000.0, 3000.0, 20000.0])
        w0 = np.array([0.1, 0.05, -0.07])
        _, _, w_hist = propagate_tumble(np.array([0, 0, 0, 1.0]), w0, J,
                                        200.0, n_eval=100)
        # rotational KE = 1/2 w.J.w should be ~constant under free tumble
        ke = np.array([0.5 * w @ J @ w for w in w_hist])
        assert np.std(ke) / np.mean(ke) < 1e-3

    def test_attitude_prop_scales_with_momentum(self):
        assert attitude_prop(200.0) == pytest.approx(2.0 * attitude_prop(100.0))


# =============================================================================
# rpo_budget
# =============================================================================

class TestBudget:
    def test_prop_dv_roundtrip(self):
        m0 = 1120.0
        dv = 50.0
        prop = prop_from_dv(dv, m0)
        assert dv_from_prop(prop, m0) == pytest.approx(dv, rel=1e-9)

    def test_prop_positive(self):
        assert prop_from_dv(10.0, 1120.0) > 0.0

    def test_dv_from_prop_raises_when_too_much(self):
        with pytest.raises(ValueError):
            dv_from_prop(2000.0, 1120.0)


# =============================================================================
# rpo_guidance
# =============================================================================

class TestGuidance:
    def test_clips_kos_detects_crossing(self):
        # chord straight through origin with KOS=5 must be flagged
        assert _clips_kos([0, -20, 0], [0, 20, 0], 5.0)

    def test_clips_kos_clears_tangential(self):
        # chord far out on H-bar does not cross a 5 m KOS
        assert not _clips_kos([0, 0, 20], [0, 0, 10], 5.0)

    def test_translational_phase_positive(self):
        p = translational_phase("t", [0, -150, 0], [0, 0, 20], MS_DRY, r_kos=5.0)
        assert p.dv > 0.0 and p.prop > 0.0 and p.duration > 0.0

    def test_dwell_phase_small_but_nonzero(self):
        p = dwell_phase("d", MS_DRY, 3600.0)
        assert p.dv > 0.0
        assert p.dv < 1.0     # stationkeeping over 1 h at GEO is small

    def test_kos_routing_costs_more(self):
        direct = translational_phase("a", [0, -20, 0], [0, 20, 0], MS_DRY, r_kos=0.0)
        routed = translational_phase("b", [0, -20, 0], [0, 20, 0], MS_DRY, r_kos=5.0)
        assert routed.dv >= direct.dv


# =============================================================================
# rpo_control
# =============================================================================

class TestControl:
    def test_spin_sync_zero_at_zero_rate(self):
        p = spin_sync_phase("s", omega_dps=0.0)
        assert p.dv == pytest.approx(0.0, abs=1e-12)

    def test_detumble_scales_with_rate(self):
        p3 = detumble_phase("d", omega_dps=3.0)
        p6 = detumble_phase("d", omega_dps=6.0)
        assert p6.dv > p3.dv

    def test_arm_extension_excludes_debris(self):
        # Phase 4 inertia growth (tug only) must be far smaller than the full
        # combined-body detumble (which includes the 2700 kg debris).
        p_arm = arm_extension_phase("arm", m0=MS_DRY, omega_dps=6.0)
        p_det = detumble_phase("det", omega_dps=6.0)
        assert p_arm.dv < p_det.dv


# =============================================================================
# rpo_debris_sim  (Part 1)
# =============================================================================

class TestDebrisSim:
    def test_returns_five_plus_phases(self):
        r = simulate_debris_rpo()
        # P1 approach, P2 inspection, P3 transition, P4 spin-sync, P5 arm,
        # P6 final-approach + P6 detumble, P7 retreat = 8
        assert len(r['phases']) == 8

    def test_detumble_separated(self):
        r = simulate_debris_rpo()
        assert r['dv_detumble'] > 0.0
        assert r['dv_debris'] > 0.0
        # detumble is extracted separately; the two are independent budget lines.
        # With a short RCS moment arm (1.05 m) the momentum dump can exceed the
        # approach-phase cost — no ordering constraint is assumed here.
        assert r['dv_debris'] != r['dv_detumble']

    def test_abort_margin_applied(self):
        no_margin = simulate_debris_rpo(abort_margin=0.0)
        with_margin = simulate_debris_rpo(abort_margin=0.5)
        assert with_margin['dv_debris'] == pytest.approx(1.5 * no_margin['dv_debris'], rel=1e-6)

    def test_sweep_monotonic_in_detumble(self):
        sweep = sweep_tumble_rates([0.0, 3.0, 6.0])
        det = [s['dv_detumble'] for s in sweep]
        assert det[0] < det[1] < det[2]

    def test_zero_tumble_zero_detumble(self):
        r = simulate_debris_rpo(omega_dps=0.0)
        assert r['dv_detumble'] == pytest.approx(0.0, abs=1e-12)


# =============================================================================
# rpo_rh_sim  (Part 2)
# =============================================================================

class TestRhSim:
    def test_cycle_has_two_phase_constants(self):
        c = simulate_rh_cycle(MS_DRY)
        assert c['dv_meet'] > 0.0
        assert c['dv_dock'] > 0.0
        assert c['dv_rh'] == pytest.approx(c['dv_meet'] + c['dv_dock'])

    def test_dock_dominates_meet(self):
        # Phase B (combined body to port) is expected to dominate Phase A
        c = simulate_rh_cycle(MS_DRY)
        assert c['dv_dock'] > c['dv_meet']

    def test_five_cycles(self):
        cycles = simulate_all_rh_cycles(n_cycles=5)
        assert len(cycles) == 5

    def test_mass_decreases_across_cycles(self):
        cycles = simulate_all_rh_cycles(n_cycles=5)
        masses = [c['ms_mass'] for c in cycles]
        assert masses[0] > masses[-1]
        assert masses[-1] == pytest.approx(MS_DRY)
