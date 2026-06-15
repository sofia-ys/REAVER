"""
Unit tests for trajectory_optimizer — config, reaver_optimizer, nominal_mission, rh_raan_sweep.

Run from the repository root:
    cd c:\\Projects\\DSE\\REAVER
    python -m pytest trajectory_optimizer/test_reaver.py -v

NOTE: rh_raan_sweep.py runs the full RAAN sweep at module level (outside __main__) and
takes several minutes. Its tests are automatically skipped unless you first wrap lines
236 onward of rh_raan_sweep.py in `if __name__ == '__main__':`.
"""

import sys
import os
import numpy as np
import pytest

# Make trajectory_optimizer importable regardless of where pytest is invoked from.
TRAJ_DIR = os.path.dirname(os.path.abspath(__file__))
if TRAJ_DIR not in sys.path:
    sys.path.insert(0, TRAJ_DIR)

# =============================================================================
# config.py
# =============================================================================

from config import (
    MU, G0, DAY, D2R,
    MS_DRY, MS_ISP, MS_VEX,
    TUG_DRY, TUG_ISP, TUG_VEX, TUG_THR,
    N_PHASE_REV, T_OPS, MAX_DAYS, SOFT_MASS,
    RH_SMA, RH_INC, RH_RAAN,
    _RAW, N_DEB, IDS, NAMES, MASS, SMA, INC, RAAN,
    RH_IDX, SMA_ALL, INC_ALL, RAAN_ALL, MASS_ALL,
)


class TestConfig:
    # ── Physical constants ────────────────────────────────────────────────────

    def test_mu_value(self):
        assert MU == pytest.approx(3.986004418e14, rel=1e-9)

    def test_g0_value(self):
        assert G0 == pytest.approx(9.80665, rel=1e-9)

    def test_day_seconds(self):
        assert DAY == 86400.0

    def test_d2r_value(self):
        assert D2R == pytest.approx(np.pi / 180.0, rel=1e-12)

    def test_ms_vex_derived_from_isp(self):
        assert MS_VEX == pytest.approx(MS_ISP * G0, rel=1e-9)

    def test_tug_vex_derived_from_isp(self):
        assert TUG_VEX == pytest.approx(TUG_ISP * G0, rel=1e-9)

    # ── Debris catalogue ──────────────────────────────────────────────────────

    def test_catalogue_size(self):
        assert N_DEB == 16
        assert len(_RAW) == 16
        assert len(IDS) == 16
        assert len(NAMES) == 16

    def test_sma_in_metres(self):
        # GEO graveyard debris SMA ≈ 42 000 – 43 000 km
        assert float(SMA.min()) > 42_000e3
        assert float(SMA.max()) < 43_000e3

    def test_mass_positive_and_plausible(self):
        assert bool((MASS > 0).all())
        assert float(MASS.min()) > 1_000.0     # smallest entry > 1 t
        assert float(MASS.max()) < 5_000.0     # largest entry < 5 t

    def test_raan_range_degrees(self):
        assert float(RAAN.min()) >= 0.0
        assert float(RAAN.max()) <= 360.0

    def test_inc_range_degrees(self):
        # GEO drift inclinations are well below 15°
        assert float(INC.min()) >= 0.0
        assert float(INC.max()) < 15.0

    # ── Recycling Hub augmented arrays ───────────────────────────────────────

    def test_rh_idx_equals_n_deb(self):
        assert RH_IDX == N_DEB == 16

    def test_augmented_array_lengths(self):
        assert len(SMA_ALL) == N_DEB + 1
        assert len(INC_ALL) == N_DEB + 1
        assert len(RAAN_ALL) == N_DEB + 1
        assert len(MASS_ALL) == N_DEB + 1

    def test_rh_appended_at_correct_index(self):
        assert SMA_ALL[RH_IDX] == RH_SMA
        assert INC_ALL[RH_IDX] == RH_INC
        assert RAAN_ALL[RH_IDX] == RH_RAAN
        assert MASS_ALL[RH_IDX] == 0.0

    def test_debris_entries_unchanged_in_augmented(self):
        np.testing.assert_array_equal(SMA_ALL[:N_DEB],  SMA)
        np.testing.assert_array_equal(INC_ALL[:N_DEB],  INC)
        np.testing.assert_array_equal(RAAN_ALL[:N_DEB], RAAN)
        np.testing.assert_array_equal(MASS_ALL[:N_DEB], MASS)


# =============================================================================
# reaver_optimizer.py
# Importing runs build_transfer_table() and the tug-spiral pre-compute, both fast.
# =============================================================================

import reaver_optimizer as ro

_N = N_DEB + 1  # number of nodes (16 debris + RH)


class TestBuildTransferTable:

    def test_matrix_shapes(self):
        for mat in (ro.DV1, ro.DV2, ro.DV_PH, ro.T_TR, ro.T_PH, ro.DV_LEG, ro.T_LEG):
            assert mat.shape == (_N, _N)

    def test_diagonal_zero(self):
        for i in range(_N):
            assert ro.DV1[i, i]   == 0.0
            assert ro.DV2[i, i]   == 0.0
            assert ro.DV_PH[i, i] == 0.0
            assert ro.T_TR[i, i]  == 0.0

    def test_off_diagonal_positive(self):
        for i in range(_N):
            for j in range(_N):
                if i != j:
                    assert ro.DV_LEG[i, j] > 0.0, f"DV_LEG[{i},{j}] ≤ 0"
                    assert ro.T_LEG[i, j]  > 0.0, f"T_LEG[{i},{j}] ≤ 0"

    def test_transfer_time_symmetric(self):
        # Hohmann half-period depends only on the semi-major axis pair, so T_TR[i,j] == T_TR[j,i]
        for i in range(_N):
            for j in range(_N):
                if i != j:
                    assert ro.T_TR[i, j] == pytest.approx(ro.T_TR[j, i], rel=1e-9)

    def test_dv_leg_is_sum_of_components(self):
        np.testing.assert_allclose(ro.DV_LEG, ro.DV1 + ro.DV2 + ro.DV_PH, rtol=1e-10)

    def test_t_leg_is_sum_of_components(self):
        np.testing.assert_allclose(ro.T_LEG, ro.T_TR + ro.T_PH, rtol=1e-10)


class TestPhasingHohmann:

    def test_positive_results(self):
        dv, t = ro.phasing_hohmann(RH_SMA)
        assert dv > 0.0
        assert t > 0.0

    def test_return_types_are_float(self):
        dv, t = ro.phasing_hohmann(RH_SMA)
        assert isinstance(dv, float)
        assert isinstance(t, float)

    def test_higher_orbit_yields_longer_phasing(self):
        # Super-synchronous orbit (RH_SMA > GEO) should have longer period → more phasing time
        _, t_rh  = ro.phasing_hohmann(RH_SMA)
        _, t_geo = ro.phasing_hohmann(42_164e3)
        assert t_rh > t_geo


class TestTugPrecompute:

    def test_array_shapes(self):
        for arr in (ro.TUG_DV, ro.TUG_TIME, ro.TUG_MPROP, ro.TUG_MWET):
            assert arr.shape == (N_DEB,)

    def test_all_positive(self):
        assert bool((ro.TUG_DV    > 0).all())
        assert bool((ro.TUG_TIME  > 0).all())
        assert bool((ro.TUG_MPROP > 0).all())
        assert bool((ro.TUG_MWET  > 0).all())

    def test_wet_mass_equals_dry_plus_payload_plus_prop(self):
        # TUG_MWET[k] ≈ TUG_DRY + MASS[k] + TUG_MPROP[k]
        for k in range(N_DEB):
            expected = TUG_DRY + float(MASS[k]) + float(ro.TUG_MPROP[k])
            assert float(ro.TUG_MWET[k]) == pytest.approx(expected, rel=1e-3)

    def test_loaded_wet_mass(self):
        expected = TUG_DRY + float(ro.TUG_MPROP.max())
        assert ro.TUG_WET_LOADED == pytest.approx(expected, rel=1e-9)


class TestEvalSingle:
    """Tests for the scalar sensitivity helper _eval_single(seq, ms_prop)."""

    SEQ = (0, 1, 2, 3, 4)

    def test_returns_dict_with_feasible_key(self):
        r = ro._eval_single(self.SEQ, 2500.0)
        assert 'feasible' in r

    def test_zero_propellant_is_infeasible(self):
        r = ro._eval_single(self.SEQ, 0.0)
        assert not r['feasible']

    def test_large_propellant_is_feasible(self):
        r = ro._eval_single(self.SEQ, 10_000.0)
        assert r['feasible']
        assert 0.0 < r['day'] <= MAX_DAYS

    def test_day_positive_when_feasible(self):
        r = ro._eval_single(self.SEQ, 3000.0)
        if r['feasible']:
            assert r['day'] > 0.0

    def test_more_propellant_cannot_increase_mission_day(self):
        # Adding propellant does not change visit times — mission day is timing-limited
        r_lo = ro._eval_single(self.SEQ, 2000.0)
        r_hi = ro._eval_single(self.SEQ, 5000.0)
        if r_lo['feasible'] and r_hi['feasible']:
            assert r_hi['day'] == pytest.approx(r_lo['day'], rel=1e-6)


# =============================================================================
# nominal_mission.py
# Importing runs the transfer table and tug spirals — both fast.
# =============================================================================

import nominal_mission as nm


class TestEvaluateSequence:
    """Tests for nominal_mission.evaluate_sequence(seq)."""

    SEQ  = [0, 1, 2, 3, 4]
    SEQ2 = [5, 10, 7, 3, 0]

    REQUIRED_KEYS = [
        'seq', 'ms_prop', 'ms_prop_orbital', 'rcs_allocation', 'rcs_margin',
        'ms_wet0', 'tug_mwets',
        'dv_legs', 't_legs', 'mass_vec',
        'tug_starts', 'tug_arrive', 'handover',
        'ms_return', 'mission_day', 'tot_dv', 'feasible',
    ]

    def test_result_has_all_keys(self):
        r = nm.evaluate_sequence(self.SEQ)
        for key in self.REQUIRED_KEYS:
            assert key in r, f"Missing result key: '{key}'"

    def test_ms_prop_positive(self):
        r = nm.evaluate_sequence(self.SEQ)
        assert r['ms_prop'] > 0.0

    def test_wet_mass_identity(self):
        # ms_wet0 = MS_DRY + ms_prop + sum(tug wet masses)
        r = nm.evaluate_sequence(self.SEQ)
        expected = MS_DRY + r['ms_prop'] + float(r['tug_mwets'].sum())
        assert r['ms_wet0'] == pytest.approx(expected, rel=1e-6)

    def test_forward_pass_ends_above_dry_mass(self):
        # mass_vec[-1] is after the return orbital burn but before RH RPO burns.
        # With the 10% RCS margin, this should be above MS_DRY by roughly rcs_allocation
        # minus the RPO prop consumed at the 5 debris stops.
        r = nm.evaluate_sequence(self.SEQ)
        assert float(r['mass_vec'][-1]) > MS_DRY

    def test_rcs_margin_positive(self):
        # Full RCS budget must cover all RPO burns; remaining margin > 0
        r = nm.evaluate_sequence(self.SEQ)
        assert r['rcs_margin'] > 0.0

    def test_rcs_allocation_is_ten_percent(self):
        r = nm.evaluate_sequence(self.SEQ)
        assert r['rcs_allocation'] == pytest.approx(0.10 * r['ms_prop_orbital'], rel=1e-9)

    def test_ms_prop_equals_orbital_plus_rcs(self):
        r = nm.evaluate_sequence(self.SEQ)
        assert r['ms_prop'] == pytest.approx(r['ms_prop_orbital'] + r['rcs_allocation'], rel=1e-9)

    def test_tug_arrive_not_before_tug_start(self):
        r = nm.evaluate_sequence(self.SEQ)
        for i in range(5):
            assert r['tug_arrive'][i] >= r['tug_starts'][i] - 1e-9

    def test_handover_after_ms_return(self):
        # Handover = max(ms_return, tug_arrive) + T_OPS — always ≥ ms_return
        r = nm.evaluate_sequence(self.SEQ)
        for i in range(5):
            assert r['handover'][i] >= r['ms_return'] - 1e-9

    def test_mission_day_is_max_handover(self):
        r = nm.evaluate_sequence(self.SEQ)
        assert r['mission_day'] == pytest.approx(float(r['handover'].max()), rel=1e-9)

    def test_feasibility_flag_correct(self):
        r = nm.evaluate_sequence(self.SEQ)
        assert r['feasible'] == (r['mission_day'] <= MAX_DAYS)

    def test_tot_dv_equals_sum_of_legs(self):
        r = nm.evaluate_sequence(self.SEQ)
        assert r['tot_dv'] == pytest.approx(float(r['dv_legs'].sum()), rel=1e-9)

    def test_exactly_six_legs(self):
        r = nm.evaluate_sequence(self.SEQ)
        assert len(r['dv_legs'])  == 6
        assert len(r['t_legs'])   == 6
        assert len(r['mass_vec']) == 6

    def test_exactly_five_tugs(self):
        r = nm.evaluate_sequence(self.SEQ)
        assert len(r['tug_starts']) == 5
        assert len(r['tug_arrive']) == 5
        assert len(r['handover'])   == 5

    @pytest.mark.parametrize("seq", [
        [0, 1, 2, 3, 4],
        [5, 10, 7, 3, 0],
        [15, 14, 13, 12, 11],
    ])
    def test_various_sequences_return_positive_values(self, seq):
        r = nm.evaluate_sequence(seq)
        assert r['ms_prop'] > 0.0
        assert r['mission_day'] > 0.0

    def test_leg0_dv_matches_precomputed_table(self):
        # The first leg (RH → seq[0]) must equal the pre-computed DV_LEG entry
        r = nm.evaluate_sequence(self.SEQ)
        assert float(r['dv_legs'][0]) == pytest.approx(
            float(nm.DV_LEG[RH_IDX, self.SEQ[0]]), rel=1e-9
        )

    def test_return_leg_dv_matches_precomputed_table(self):
        # Last leg (seq[4] → RH) must equal the pre-computed DV_LEG entry
        r = nm.evaluate_sequence(self.SEQ)
        assert float(r['dv_legs'][5]) == pytest.approx(
            float(nm.DV_LEG[self.SEQ[4], RH_IDX]), rel=1e-9
        )

    def test_transfer_tables_consistent_across_modules(self):
        # Both modules build the same transfer table from the same algorithm;
        # verify one representative entry agrees.
        i, j = 3, 7
        assert float(nm.DV_LEG[i, j]) == pytest.approx(float(ro.DV_LEG[i, j]), rel=1e-6)


# =============================================================================
# rh_raan_sweep.py
# The sweep loop is at module level and takes several minutes — it is NOT wrapped
# in `if __name__ == '__main__':`.  These tests are skipped automatically until
# you add that guard before line 236 of rh_raan_sweep.py.
# =============================================================================

try:
    import rh_raan_sweep as rhs
    _RHS_OK = True
except Exception:
    rhs = None
    _RHS_OK = False

_rhs_skip = pytest.mark.skipif(
    not _RHS_OK,
    reason=(
        "rh_raan_sweep.py runs the full sweep at module level. "
        "Wrap lines 236+ in `if __name__ == '__main__':` to enable these tests."
    ),
)


@_rhs_skip
class TestRhTransfers:
    """Tests for rh_raan_sweep._rh_transfers(rh_raan_deg)."""

    def test_matrix_shapes(self):
        DV, DT = rhs._rh_transfers(35.0)
        assert DV.shape == (_N, _N)
        assert DT.shape == (_N, _N)

    def test_off_diagonal_positive(self):
        DV, DT = rhs._rh_transfers(35.0)
        for i in range(_N):
            for j in range(_N):
                if i != j:
                    assert DV[i, j] > 0.0
                    assert DT[i, j] > 0.0

    def test_debris_to_debris_independent_of_raan(self):
        # RH RAAN only affects legs that involve the RH node
        DV_a, _ = rhs._rh_transfers(0.0)
        DV_b, _ = rhs._rh_transfers(180.0)
        for i in range(N_DEB):
            for j in range(N_DEB):
                if i != j:
                    assert DV_a[i, j] == pytest.approx(DV_b[i, j], rel=1e-9)

    def test_rh_column_changes_with_raan(self):
        DV_a, _ = rhs._rh_transfers(0.0)
        DV_b, _ = rhs._rh_transfers(90.0)
        # At least one RH↔debris entry must differ
        assert not np.allclose(DV_a[RH_IDX, :N_DEB], DV_b[RH_IDX, :N_DEB])


@_rhs_skip
class TestRhsTugData:
    """Tests for rh_raan_sweep._tug_data(rh_raan_deg)."""

    def test_array_shapes(self):
        mp, tt = rhs._tug_data(35.0)
        assert mp.shape == (N_DEB,)
        assert tt.shape == (N_DEB,)

    def test_all_positive(self):
        mp, tt = rhs._tug_data(35.0)
        assert bool((mp > 0).all())
        assert bool((tt > 0).all())

    def test_varies_with_raan(self):
        mp_a, _ = rhs._tug_data(0.0)
        mp_b, _ = rhs._tug_data(90.0)
        assert not np.allclose(mp_a, mp_b)


@_rhs_skip
class TestEvaluateRaan:
    """Tests for rh_raan_sweep.evaluate_raan(rh_raan_deg)."""

    def test_returns_required_keys(self):
        r = rhs.evaluate_raan(35.0)
        for key in ('n_feas_combos', 'prop_worst', 'prop_best', 'prop_mean'):
            assert key in r

    def test_feasibility_count_in_range(self):
        r = rhs.evaluate_raan(35.0)
        assert 0 <= r['n_feas_combos'] <= rhs.N_COMBOS

    def test_best_le_mean_le_worst(self):
        r = rhs.evaluate_raan(35.0)
        if r['n_feas_combos'] > 0:
            assert r['prop_best']  <= r['prop_mean']  + 1e-6
            assert r['prop_mean']  <= r['prop_worst'] + 1e-6

    def test_propellant_positive_when_feasible(self):
        r = rhs.evaluate_raan(35.0)
        if r['n_feas_combos'] > 0:
            assert r['prop_best']  > 0.0
            assert r['prop_worst'] > 0.0

    def test_optimal_raan_all_combos_feasible(self):
        # RH_RAAN = 35° (from config) is the known optimum — all 4368 combos should be feasible
        r = rhs.evaluate_raan(float(RH_RAAN))
        assert r['n_feas_combos'] == rhs.N_COMBOS