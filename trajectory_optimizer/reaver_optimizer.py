"""
REAVER Mission Trajectory Optimizer  —  Group 7
=================================================
Vectorised NumPy brute-force over all C(16,5)x5! = 524,160 sequences.

Mission architecture:
  - Mothership (chemical, Isp=253s) visits 5 debris sequentially
  - Each capture hands off to one tug (E-prop, Isp=1850s, 0.2N)
  - Tugs spiral independently back to RH (parallel)
  - Mothership returns to RH after all 5 captures
  - Mothership meets each tug at RH for final handover
  - All 5 handovers within 365 days

Transfer model per mothership leg:
  Burn 1: Hohmann departure (pure, no plane change)
  Burn 2: Combined circularisation + plane change at apoapsis
          (plane angle via spherical law of cosines, RAAN+inc combined)
  Burn 3: Double-Hohmann phasing (15 revs, closes 90° avg phase gap)

Tug model: Edelbaum low-thrust optimal transfer (iterative mass sizing)
"""

import numpy as np
from itertools import combinations, permutations
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from collections import Counter
import sys, os, time, warnings
warnings.filterwarnings('ignore')

from config import *
from reaver_core import (build_transfer_table, compute_tug_spirals, phasing_hohmann,
                         DV1, DV2, DV_PH, T_TR, T_PH, DV_LEG, T_LEG, T_LEG_FINITE,
                         finite_burn_time, finite_burn_info,
                         TUG_DV, TUG_TIME, TUG_MPROP, TUG_MWET, TUG_WET_LOADED)

# =============================================================================
# VECTORISED SEQUENCE EVALUATOR
# =============================================================================

def evaluate_all_sequences():
    """
    Evaluate all C(16,5)*5! ordered sequences using pre-computed tables.
    Returns arrays of results for feasible sequences.
    """
    # Generate all ordered 5-tuples
    print(f"\n  Building sequence index...", end=' ', flush=True)
    sequences = []
    for combo in combinations(range(N_DEB), 5):
        for perm in permutations(combo):
            sequences.append(perm)
    sequences = np.array(sequences, dtype=np.int32)   # (524160, 5)
    n_seq = len(sequences)
    print(f"{n_seq:,} sequences")

    # ── Phase 1: Mothership chain RH->D1->D2->D3->D4->D5 ──────────────────
    # Track mass and time through 5 legs + return leg
    print("  Evaluating mothership chains...", end=' ', flush=True)
    t1 = time.time()

    # from_node for each step: RH=16, then d[0..4]
    # shape: (n_seq, 6) for legs 0..5 (leg5 = return to RH)
    # leg from_nodes: [RH, d0, d1, d2, d3, d4]
    # leg to_nodes:   [d0, d1, d2, d3, d4, RH]

    from_nodes = np.zeros((n_seq, 6), dtype=np.int32)
    to_nodes   = np.zeros((n_seq, 6), dtype=np.int32)
    from_nodes[:, 0] = RH_IDX
    from_nodes[:, 1:5] = sequences[:, :4]
    from_nodes[:, 5] = sequences[:, 4]
    to_nodes[:, :5]   = sequences
    to_nodes[:, 5]    = RH_IDX

    # ΔV and time for each of the 6 legs
    dv_legs = DV_LEG[from_nodes, to_nodes]   # (n_seq, 6)
    t_legs  = T_LEG_FINITE[from_nodes, to_nodes]   # (n_seq, 6)

    # Per-sequence tug wet masses: all 5 tugs sized to the worst-case debris in each sequence
    tug_prop_uniform = TUG_MPROP[sequences].max(axis=1, keepdims=True)   # (n_seq, 1)
    tug_mwet_seq = (TUG_DRY + np.repeat(tug_prop_uniform, 5, axis=1)).astype(np.float64)  # (n_seq, 5)

    # Backward pass: orbital prop sized so final mass = MS_DRY (transfers only)
    m_req = np.full(n_seq, MS_DRY, dtype=np.float64)
    m_req = m_req * np.exp(dv_legs[:, 5] / MS_VEX)       # un-burn return leg
    for leg in range(4, -1, -1):
        m_req += tug_mwet_seq[:, leg]                      # re-attach tug
        m_req = m_req * np.exp(dv_legs[:, leg] / MS_VEX)  # un-burn this leg
    ms_prop_seq  = m_req - MS_DRY - tug_mwet_seq.sum(axis=1)  # (n_seq,) orbital prop
    rcs_alloc_seq = MS_RCS_MARGIN * ms_prop_seq                 # 10 % RCS margin
    ms_wet0      = m_req + rcs_alloc_seq                        # updated initial mass

    # Forward pass: orbital burns + RPO burns interleaved
    cum_dv = np.cumsum(dv_legs, axis=1)       # (n_seq, 6) — ΔV totals for reporting
    mass   = ms_wet0.copy()
    mass_after = np.zeros((n_seq, 6))
    for leg in range(6):
        mass = mass * np.exp(-dv_legs[:, leg] / MS_VEX)
        mass_after[:, leg] = mass
        if leg < 5:
            # RPO at debris (inspect+capture + detumble) during T_OPS
            mass = mass * np.exp(-(DV_RPO_DEBRIS + DV_RPO_DETUMBLE) / MS_VEX)
            mass -= tug_mwet_seq[:, leg]      # tug detaches after ops
    # 5 × RH proximity docking (tug-meet + dock) after MS return
    mass = mass * np.exp(-5.0 * DV_RPO_RH / MS_VEX)
    mass_final = mass   # true final mass after all orbital + RPO burns

    # Cumulative time with T_OPS added after each of the 5 captures
    # t_ops only added after legs 0-4 (not the return leg 5)
    t_ops_arr = np.array([T_OPS]*5 + [0.0])
    cum_time  = np.cumsum(t_legs + t_ops_arr[None, :], axis=1)   # (n_seq, 6)

    # Handoff day for each tug = time after leg i + T_OPS
    # tug_start[i] = cum_time[:, i]  (after leg i including T_OPS)
    tug_start = cum_time[:, :5]   # (n_seq, 5)  days when each tug starts spiral

    # MS return day = cum_time[:, 5]  (no T_OPS after return)
    ms_return = cum_time[:, 5]    # (n_seq,)

    # Total mothership ΔV
    tot_dv   = cum_dv[:, 5]        # (n_seq,)
    tot_prop = ms_prop_seq   # orbital propellant per sequence

    print(f"done ({time.time()-t1:.2f}s)")

    # ── Phase 2: Tug spirals (parallel) ────────────────────────────────────
    print("  Evaluating tug spirals...", end=' ', flush=True)
    t2 = time.time()

    # Tug arrival = tug_start + TUG_TIME[debris_index]
    tug_spiral = TUG_TIME[sequences]             # (n_seq, 5)
    tug_arrive = tug_start + tug_spiral          # (n_seq, 5)

    # Sequential RH handovers: process tugs in arrival order per sequence;
    # each handover waits for the previous one to finish AND for the tug to arrive.
    sort_idx = np.argsort(tug_arrive, axis=1)          # (n_seq, 5)
    handover  = np.zeros((n_seq, 5), dtype=np.float64)
    t_rh      = ms_return.copy()                        # (n_seq,) current RH free time
    row_idx   = np.arange(n_seq)
    for rank in range(5):
        col    = sort_idx[:, rank]                      # (n_seq,) which tug is rank-th to arrive
        ta     = tug_arrive[row_idx, col]               # (n_seq,) its arrival day
        t_start = np.maximum(t_rh, ta)                  # wait for tug if needed
        ho      = t_start + T_OPS
        handover[row_idx, col] = ho                     # place back in original tug order
        t_rh    = ho                                    # next handover can't start before this

    # Mission completion = last handover (= t_rh after loop)
    mission_day = handover.max(axis=1)           # (n_seq,)

    print(f"done ({time.time()-t2:.2f}s)")

    # ── Feasibility & scoring ───────────────────────────────────────────────
    time_ok  = mission_day <= MAX_DAYS
    feasible = time_ok
    f_idx    = np.where(feasible)[0]

    print(f"  Feasible: {len(f_idx):,}  |  Infeasible: {n_seq-len(f_idx):,}")

    if len(f_idx) == 0:
        return None

    # Extract feasible subset
    f_seq    = sequences[f_idx]           # (nf, 5)
    f_dv     = tot_dv[f_idx]
    f_prop   = tot_prop[f_idx]
    f_day    = mission_day[f_idx]
    f_mass   = MASS[f_seq].sum(axis=1)
    f_mret   = ms_return[f_idx]
    f_tstart = tug_start[f_idx]
    f_tarr   = tug_arrive[f_idx]
    f_ho     = handover[f_idx]
    f_heavy  = (MASS[f_seq] > SOFT_MASS).any(axis=1)
    f_mfinal = mass_final[f_idx]          # true final mass after all orbital + RPO burns
    f_rcs    = rcs_alloc_seq[f_idx]       # RCS budget per sequence

    # Pareto score on all feasible permutations (used to pick best ordering per combo)
    dv_n_all   = (f_dv  - f_dv.min())  / (f_dv.max()  - f_dv.min()  + 1e-9)
    time_n_all = (f_day - f_day.min()) / (f_day.max() - f_day.min() + 1e-9)
    score_all  = 0.50*dv_n_all + 0.50*time_n_all

    # Reduce to best permutation per combination (C(16,5)=4368 unique target sets)
    # Each combo has up to 120 orderings; keep only the best-scoring one so that
    # worst-case analysis reflects hard target sets, not bad visit orders.
    sorted_combos = np.sort(f_seq, axis=1)
    combo_keys = [tuple(r) for r in sorted_combos]
    best_per_combo = {}
    for i, key in enumerate(combo_keys):
        if key not in best_per_combo or score_all[i] < score_all[best_per_combo[key]]:
            best_per_combo[key] = i
    c_idx = np.array(list(best_per_combo.values()))
    print(f"  Unique combinations: {len(c_idx):,}  "
          f"(best ordering per combo, from {len(f_idx):,} feasible permutations)")

    # Re-normalise score on the reduced set so ranking reflects the combo space
    c_dv   = f_dv[c_idx];  c_day = f_day[c_idx]
    dv_n   = (c_dv  - c_dv.min())  / (c_dv.max()  - c_dv.min()  + 1e-9)
    time_n = (c_day - c_day.min()) / (c_day.max() - c_day.min() + 1e-9)
    score  = 0.50*dv_n + 0.50*time_n
    order  = np.argsort(score)
    sel    = c_idx[order]

    tug_prop_sel = TUG_MPROP[f_seq[sel]].max(axis=1) * 5   # uniform sizing: 5× worst tug in sequence
    ms_prop_sel  = ms_prop_seq[f_idx][sel]
    rcs_sel      = f_rcs[sel]
    results = {
        'sequences':   f_seq[sel],
        'score':       score[order],
        'total_dv':    f_dv[sel],
        'prop_used':   f_prop[sel],              # orbital prop (sizing basis)
        'rcs_alloc':   rcs_sel,                  # RCS budget (10 % of orbital prop)
        'prop_margin': f_mfinal[sel] - MS_DRY,  # remaining RCS at end of mission
        'mission_day': f_day[sel],
        'ms_return':   f_mret[sel],
        'mass_removed':f_mass[sel],
        'tug_start':   f_tstart[sel],
        'tug_arrive':  f_tarr[sel],
        'handover':    f_ho[sel],
        'has_heavy':   f_heavy[sel],
        'n_feasible':  len(c_idx),
        'n_total':     n_seq,
        'ms_prop':     ms_prop_sel,
        'tug_prop':    tug_prop_sel,
        'total_prop':  ms_prop_sel + tug_prop_sel,
        'dv_legs':     dv_legs[f_idx][sel],
        't_legs':      t_legs[f_idx][sel],
        'cum_time':    cum_time[f_idx][sel],
    }
    return results

# =============================================================================
# RESULTS PRINTER
# =============================================================================

def print_top(res, n=3, start=0, label=None):
    if label is None:
        label = f"TOP {n}" if start == 0 else f"WORST-CASE (rank {start+1}/{res['n_feasible']})"
    print("\n"+"="*70)
    print(f"  {label}  (50% ΔV + 50% time)")
    print("="*70)

    for rank in range(start, min(start + n, res['n_feasible'])):
        seq   = res['sequences'][rank]
        score = res['score'][rank]
        heavy = '  ⚠ contains debris >2000kg' if res['has_heavy'][rank] else ''
        print(f"\n{'─'*70}")
        print(f"  RANK {rank+1}  |  Score {score:.4f}{heavy}")
        print(f"{'─'*70}")
        print("  Visit order:")
        for i, idx in enumerate(seq):
            flag = ' ⚠' if MASS[idx] > SOFT_MASS else ''
            print(f"    {i+1}. {NAMES[idx]:<33} {MASS[idx]:>7.0f} kg{flag}")

        dv_l = res['dv_legs'][rank]    # 6 legs
        t_l  = res['t_legs'][rank]
        leg_names = [NAMES[j] for j in seq] + ['→ Recycling Hub']
        print(f"\n  {'Leg':<3} {'To':<28} {'ΔV tot':>8} {'ΔV Hohm':>8} "
              f"{'ΔV plane':>9} {'ΔV ph':>7} {'Transf':>7} {'Phase':>7}")
        print(f"  {'':─<3} {'':─<28} {'':─>8} {'':─>8} "
              f"{'':─>9} {'':─>7} {'':─>7} {'':─>7}")
        nodes_from = [RH_IDX] + list(seq)
        nodes_to   = list(seq) + [RH_IDX]
        for li in range(6):
            fi = nodes_from[li]; ti = nodes_to[li]
            nm = leg_names[li][:28]
            dv1_v  = DV1[fi,ti]; dv2_v = DV2[fi,ti]
            dvph_v = DV_PH[fi,ti]
            print(f"  {li+1 if li<5 else 'R':<3} {nm:<28} "
                  f"{dv_l[li]:>7.1f}m {dv1_v:>7.1f}m "
                  f"{dv2_v:>8.1f}m {dvph_v:>6.1f}m "
                  f"{T_TR[fi,ti]:>6.1f}d {T_PH[fi,ti]:>6.1f}d")

        print(f"\n  {'Tug':<3} {'Debris':<28} "
              f"{'ΔV Edel':>8} {'Prop kg':>8} {'Spiral d':>9} "
              f"{'Arrive d':>9} {'Handover':>9}")
        print(f"  {'':─<3} {'':─<28} {'':─>8} {'':─>8} "
              f"{'':─>9} {'':─>9} {'':─>9}")
        tug_prop_uniform = TUG_MPROP[list(seq)].max()
        for i, idx in enumerate(seq):
            worst_flag = '  ← sizing driver' if TUG_MPROP[idx] == tug_prop_uniform else ''
            print(f"  {i+1:<3} {NAMES[idx]:<28} "
                  f"{TUG_DV[idx]:>7.1f}m "
                  f"{tug_prop_uniform:>7.1f}kg "
                  f"{TUG_TIME[idx]:>8.1f}d "
                  f"{res['tug_arrive'][rank,i]:>8.1f}d "
                  f"{res['handover'][rank,i]:>8.1f}d{worst_flag}")

        rcs_pct = 100.0 * res['prop_margin'][rank] / res['rcs_alloc'][rank]
        print(f"\n  ┌──────────────────────────────────────────────────────────────┐")
        print(f"  │  Mothership ΔV total       : {res['total_dv'][rank]:>8.1f} m/s               │")
        print(f"  │  MS orbital propellant     : {res['prop_used'][rank]:>8.1f} kg                │")
        print(f"  │  RCS margin (+10% budget)  : {res['rcs_alloc'][rank]:>8.1f} kg                │")
        print(f"  │  MS total initial prop     : {res['prop_used'][rank]+res['rcs_alloc'][rank]:>8.1f} kg                │")
        print(f"  │  RCS margin remaining      : {res['prop_margin'][rank]:>8.1f} kg  ({rcs_pct:5.1f}% of budget)  │")
        print(f"  │  MS returns to RH day      : {res['ms_return'][rank]:>8.1f}                  │")
        print(f"  │  Total debris mass         : {res['mass_removed'][rank]:>8.0f} kg                │")
        print(f"  │  Mission completion        : {res['mission_day'][rank]:>8.1f} days  ✓         │")
        print(f"  └──────────────────────────────────────────────────────────────┘")

# =============================================================================
# TARGET FREQUENCY
# =============================================================================

def print_target_frequency(res, top_n=50):
    print("\n" + "="*70)
    print(f"  TARGET SELECTION FREQUENCY — Top {top_n} Sequences (with positions)")
    print("="*70)

    n = min(top_n, res['n_feasible'])
    pos_counts  = np.zeros((N_DEB, 5), dtype=int)
    total_counts = np.zeros(N_DEB, dtype=int)

    for r in range(n):
        for pos, idx in enumerate(res['sequences'][r]):
            pos_counts[idx, pos] += 1
            total_counts[idx] += 1

    appeared = [(i, total_counts[i]) for i in range(N_DEB) if total_counts[i] > 0]
    appeared.sort(key=lambda x: -x[1])

    print(f"\n  {'Debris':<33} {'Total':>5}  "
          f"{'P1':>4} {'P2':>4} {'P3':>4} {'P4':>4} {'P5':>4}")
    print(f"  {'':─<33} {'':─>5}  "
          f"{'':─>4} {'':─>4} {'':─>4} {'':─>4} {'':─>4}")
    for idx, total in appeared:
        pc = pos_counts[idx]
        print(f"  {NAMES[idx]:<33} {total:>5}  "
              f"{pc[0]:>4} {pc[1]:>4} {pc[2]:>4} {pc[3]:>4} {pc[4]:>4}")

# =============================================================================
# SENSITIVITY
# =============================================================================

def print_mass_timeline(res, rank, label):
    seq        = res['sequences'][rank]
    ms_prop    = res['ms_prop'][rank]         # orbital prop
    rcs_alloc  = res['rcs_alloc'][rank]       # RCS budget (10 %)
    tug_mwets  = np.full(5, TUG_DRY + TUG_MPROP[list(seq)].max())
    nf_        = [RH_IDX] + list(seq)
    nt_        = list(seq) + [RH_IDX]
    m0         = MS_DRY + ms_prop + rcs_alloc + tug_mwets.sum()

    W = 64
    print("\n" + "="*98)
    print(f"  MOTHERSHIP MASS TIMELINE — {label}")
    print("="*98)
    print("  " + "  ->  ".join(NAMES[k].split('(')[0].strip()[:14] for k in seq))
    print(f"\n  Initial: {MS_DRY:.0f} kg dry  +  {ms_prop:.1f} kg orbital prop"
          f"  +  {rcs_alloc:.1f} kg RCS (+10%)  +  {tug_mwets.sum():.1f} kg tugs  =  {m0:.1f} kg")
    print(f"\n  {'Day':>7}  {'Event':<{W}}  {'D Mass kg':>10}  {'Mass kg':>9}")
    print(f"  {'':->7}  {'':─<{W}}  {'':->10}  {'':->9}")

    m = m0
    t = 0.0

    def prow(day, event, delta=None):
        nonlocal m
        if delta is not None:
            m += delta
            ds = f"{delta:+.1f}"
        else:
            ds = "—"
        print(f"  {day:>7.1f}  {event:<{W}}  {ds:>10}  {m:>9.1f}")

    prow(t, "Depart Recycling Hub")

    for i in range(6):
        fi, ti = nf_[i], nt_[i]
        dest = NAMES[ti].split('(')[0].strip()[:28] if ti < N_DEB else 'Recycling Hub'
        asc  = SMA_ALL[ti] >= SMA_ALL[fi]
        d1, d2, dph = DV1[fi,ti], DV2[fi,ti], DV_PH[fi,ti]

        # Burn 1 — departure
        m_b1 = m
        dm = -m * (1 - np.exp(-d1 / MS_VEX))
        nf1, no1 = finite_burn_info(d1, m_b1)
        b1 = (f"  Burn 1  Hohmann dep       ({d1:6.1f} m/s | {nf1} fires, {no1} orbits)" if asc else
              f"  Burn 1  dep + plane chg   ({d1:6.1f} m/s | {nf1} fires, {no1} orbits)")
        prow(t, b1, dm)
        t += finite_burn_time(d1, m_b1, SMA_ALL[fi])

        t += T_TR[fi, ti]

        # Burn 2 — at apoapsis / arrival
        m_b2 = m
        dm = -m * (1 - np.exp(-d2 / MS_VEX))
        nf2, no2 = finite_burn_info(d2, m_b2)
        b2 = (f"  Burn 2  circ + plane chg  ({d2:6.1f} m/s | {nf2} fires, {no2} orbits)" if asc else
              f"  Burn 2  circularise        ({d2:6.1f} m/s | {nf2} fires, {no2} orbits)")
        prow(t, b2, dm)
        t += finite_burn_time(d2, m_b2, SMA_ALL[ti])

        # Burn 3 — phasing
        m_b3 = m
        dm_ph = -m * (1 - np.exp(-dph / MS_VEX))
        nf3, no3 = finite_burn_info(dph, m_b3)
        prow(t, f"  Burn 3  phasing           ({dph:6.1f} m/s | {nf3} fires, {no3} orbits)", dm_ph)
        t += finite_burn_time(dph, m_b3, SMA_ALL[ti])
        t += T_PH[fi, ti]
        prow(t, f"  Phase coast               ({N_PHASE_REV} rev × {T_PH[fi,ti]/N_PHASE_REV:.2f} d = {T_PH[fi,ti]:.1f} d)", None)

        if i < 5:
            prow(t, f"  Rendezvous: {dest}", None)
            # Inspect + capture on arrival day; detumble ~1 day later after ground relay
            dm_rpo1 = -m * (1 - np.exp(-DV_RPO_DEBRIS   / MS_VEX))
            prow(t,       f"  RPO inspect + capture     ({DV_RPO_DEBRIS:.2f} m/s, +50% abort)", dm_rpo1)
            dm_rpo2 = -m * (1 - np.exp(-DV_RPO_DETUMBLE / MS_VEX))
            prow(t + 1.0, f"  RPO detumble              ({DV_RPO_DETUMBLE:.2f} m/s, ACS body dump)", dm_rpo2)
            t += T_OPS
            prow(t, f"  Tug {i+1} released ({tug_mwets[i]:.1f} kg)", -tug_mwets[i])
        else:
            prow(t, "  Arrive Recycling Hub", None)
            # Sequential handovers: process tugs in arrival order; each waits for
            # the previous handover to complete AND for the tug to arrive at RH.
            arrive_order_rh = np.argsort(res['tug_arrive'][rank])
            t_rh = t
            for j_rank in range(5):
                j = int(arrive_order_rh[j_rank])
                ta = float(res['tug_arrive'][rank, j])
                t_next = max(t_rh, ta)
                if t_next > t_rh + 0.1:
                    prow(t_next, f"  Wait: tug {j+1} arriving at RH (d {ta:.1f})", None)
                t_rh = t_next
                dm_meet = -m * (1 - np.exp(-DV_RPO_TUG_MEET / MS_VEX))
                prow(t_rh,       f"  RPO tug-meet tug {j+1}        ({DV_RPO_TUG_MEET:.2f} m/s)", dm_meet)
                dm_dock = -m * (1 - np.exp(-DV_RPO_RH_DOCK  / MS_VEX))
                prow(t_rh + 0.5, f"  RPO dock to RH tug {j+1}      ({DV_RPO_RH_DOCK:.2f} m/s)", dm_dock)
                t_rh += T_OPS
                prow(t_rh, f"  Handover {j+1} complete", None)

    print(f"  {'':->7}  {'':─<{W}}  {'':->10}  {'':->9}")
    rcs_remaining = m - MS_DRY
    rcs_pct = 100.0 * rcs_remaining / rcs_alloc
    print(f"\n  Orbital prop    : {ms_prop:.1f} kg"
          f"  |  RCS budget (+10%) : {rcs_alloc:.1f} kg"
          f"  |  Total init prop : {ms_prop+rcs_alloc:.1f} kg")
    print(f"  RCS remaining   : {rcs_remaining:.1f} kg  ({rcs_pct:.1f}% of budget)"
          f"  |  Mission day : {res['mission_day'][rank]:.1f} d")


# =============================================================================
# PROPELLANT COMPARISON TABLE
# =============================================================================

def print_prop_comparison(res, wp, wtp):
    nf = res['n_feasible']
    cases = [
        ('Best case',           0),
        ('Worst score',         nf - 1),
        ('Worst MS prop',       wp),
        ('Worst total prop',    wtp),
    ]
    W = 18
    print("\n" + "="*74)
    print("  PROPELLANT COMPARISON — key design cases")
    print("="*74)
    header = f"  {'Metric':<22}" + "".join(f"  {lbl:>{W}}" for lbl, _ in cases)
    print(header)
    print("  " + "─"*22 + ("  " + "─"*W) * len(cases))

    def row(label, vals):
        print(f"  {label:<22}" + "".join(f"  {v:>{W}}" for v in vals))

    row("Rank",              [f"{r+1}/{nf}" for _, r in cases])
    row("MS orbital [kg]",  [f"{res['ms_prop'][r]:.1f}"   for _, r in cases])
    row("RCS budget [kg]",  [f"{res['rcs_alloc'][r]:.1f}" for _, r in cases])
    row("MS total init [kg]",[f"{res['ms_prop'][r]+res['rcs_alloc'][r]:.1f}" for _, r in cases])
    row("RCS remaining [kg]",[f"{res['prop_margin'][r]:.1f}" for _, r in cases])
    row("Tug prop [kg]",    [f"{res['tug_prop'][r]:.1f}"  for _, r in cases])
    row("Total prop [kg]",  [f"{res['total_prop'][r]:.1f}" for _, r in cases])
    row("Total DV [m/s]",   [f"{res['total_dv'][r]:.1f}"  for _, r in cases])
    row("Mission day",      [f"{res['mission_day'][r]:.1f}" for _, r in cases])
    row("Debris mass [kg]", [f"{res['mass_removed'][r]:.0f}" for _, r in cases])
    print()
    for lbl, r in cases:
        seq_names = " -> ".join(NAMES[k].split('(')[0].strip()[:10] for k in res['sequences'][r])
        print(f"  {lbl:<22}  {seq_names}")


# =============================================================================
# TUG PROPELLANT ANALYSIS  (16 individual debris → RH transfers)
# =============================================================================

def print_tug_analysis():
    """Print propellant and timing for each of the 16 tug debris→RH spirals."""
    order = np.argsort(TUG_MPROP)[::-1]   # descending by propellant

    print("\n" + "="*80)
    print("  TUG PROPELLANT ANALYSIS  —  all 16 debris → Recycling Hub spirals")
    print(f"  Tug dry {TUG_DRY} kg  |  Isp {TUG_ISP} s  |  Thrust {TUG_THR} N  "
          f"|  RH RAAN {RH_RAAN}°  inc {RH_INC}°")
    print(f"  Worst-case loaded tug wet mass (mothership carry): "
          f"{TUG_WET_LOADED:.1f} kg  (TUG_DRY + max prop)")
    print("="*80)
    print(f"\n  {'#':<3} {'Debris':<33} {'Debris kg':>9}  "
          f"{'ΔV m/s':>7}  {'Prop kg':>7}  {'Dry+Prop kg':>11}  "
          f"{'Wet+Deb kg':>10}  {'Time d':>7}")
    print(f"  {'':─<3} {'':─<33} {'':─>9}  "
          f"{'':─>7}  {'':─>7}  {'':─>11}  {'':─>10}  {'':─>7}")

    for rank, k in enumerate(order):
        flag = '  ← worst' if rank == 0 else ''
        print(f"  {rank+1:<3} {NAMES[k]:<33} {MASS[k]:>9.0f}  "
              f"{TUG_DV[k]:>7.1f}  {TUG_MPROP[k]:>7.1f}  "
              f"{TUG_DRY + TUG_MPROP[k]:>11.1f}  "
              f"{TUG_MWET[k]:>10.1f}  {TUG_TIME[k]:>7.1f}{flag}")

    wk = order[0]
    print(f"\n  Worst-case tug  : {NAMES[wk]}")
    print(f"    Debris mass   : {MASS[wk]:.0f} kg")
    print(f"    Spiral ΔV     : {TUG_DV[wk]:.1f} m/s")
    print(f"    Propellant    : {TUG_MPROP[wk]:.1f} kg")
    print(f"    Spiral time   : {TUG_TIME[wk]:.1f} days")
    print(f"    Wet mass (w/ debris) : {TUG_MWET[wk]:.1f} kg")


# =============================================================================
# TUG PROPELLANT MARGINS
# =============================================================================

def print_tug_margins(res, rank):
    """For each tug in the sequence, show propellant required vs uniform budget."""
    seq              = res['sequences'][rank]
    tug_prop_uniform = float(TUG_MPROP[list(seq)].max())
    nf               = res['n_feasible']

    print("\n" + "="*76)
    print(f"  TUG PROPELLANT MARGINS — rank {rank+1}/{nf}  "
          f"(uniform budget: {tug_prop_uniform:.1f} kg/tug)")
    print("="*76)
    print(f"  {'Tug':<3} {'Debris':<33} {'Req kg':>8} {'Budget kg':>10} "
          f"{'Margin kg':>10} {'Margin %':>9}")
    print(f"  {'':─<3} {'':─<33} {'':─>8} {'':─>10} {'':─>10} {'':─>9}")

    for i, idx in enumerate(seq):
        req    = float(TUG_MPROP[idx])
        margin = tug_prop_uniform - req
        pct    = 100.0 * margin / tug_prop_uniform
        flag   = '  ← sizing driver' if req == tug_prop_uniform else ''
        print(f"  {i+1:<3} {NAMES[idx]:<33} {req:>8.1f} {tug_prop_uniform:>10.1f} "
              f"{margin:>10.1f} {pct:>8.1f}%{flag}")

    total_loaded = 5 * tug_prop_uniform
    total_needed = float(TUG_MPROP[list(seq)].sum())
    print(f"\n  Total loaded  : {total_loaded:>8.1f} kg  (5 × {tug_prop_uniform:.1f} kg)")
    print(f"  Total needed  : {total_needed:>8.1f} kg")
    print(f"  Total margin  : {total_loaded - total_needed:>8.1f} kg")


# =============================================================================
# PLOTTING
# =============================================================================

def make_plots(res, save_path=None, wp=None):
    if save_path is None:
        save_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 'figures', 'reaver_optimizer', 'reaver_optimizer_results.png')
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    BG='#0d1117'; CB='#161b22'; TC='#c9d1d9'; MU_='#8b949e'; GR='#21262d'
    A1='#58a6ff'; A2='#3fb950'; A3='#f78166'; A4='#d2a8ff'; A5='#ffa657'
    LC=[A1,A2,A4,A3,A5,'#79c0ff']; TC_=[A1,A2,A4,A3,A5]

    nf = res['n_feasible']
    _wp = int(np.argmax(res['prop_used'])) if wp is None else wp   # worst MS prop rank

    fig = plt.figure(figsize=(20,11)); fig.patch.set_facecolor(BG)
    gs  = GridSpec(2,3,figure=fig,hspace=0.52,wspace=0.36,
                   height_ratios=[1,1.3])

    def sax(ax,title=''):
        ax.set_facecolor(CB)
        [s.set_edgecolor(GR) for s in ax.spines.values()]
        ax.tick_params(colors=MU_,labelsize=8)
        ax.xaxis.label.set_color(MU_); ax.yaxis.label.set_color(MU_)
        ax.grid(True,color=GR,lw=0.5,alpha=0.7)
        if title: ax.set_title(title,color=TC,fontsize=9,fontweight='bold',pad=8)

    fig.text(0.5,0.980,'REAVER Mission Optimizer — Results Dashboard (Worst-Case MS Propellant)',
             ha='center',color=TC,fontsize=15,fontweight='bold')
    fig.text(0.5,0.958,
             f'MS: dry {MS_DRY}kg, Isp={MS_ISP}s  |  '
             f'Tug: dry {TUG_DRY}kg, Isp={TUG_ISP}s, T={TUG_THR}N  |  '
             f'T_ops={T_OPS}d/debris  |  Hohmann phasing {N_PHASE_REV} revs  |  ≤{MAX_DAYS:.0f}d  |  '
             f'Rank {_wp+1}/{nf} (worst MS prop)',
             ha='center',color=MU_,fontsize=8)

    seqW = res['sequences'][_wp]

    # ── 1. Mission timeline ribbon (full width) ───────────────────────────
    ax5 = fig.add_subplot(gs[0,:])
    sax(ax5,f'Worst-case MS prop (rank {_wp+1}/{nf}) — full mission timeline')
    events=[(0,'Depart RH',MU_)]
    cum=0.0
    for i,idx in enumerate(seqW):
        fi=[RH_IDX]+list(seqW)
        ti=list(seqW)+[RH_IDX]
        cum+=T_TR[fi[i],ti[i]]+T_PH[fi[i],ti[i]]
        events.append((cum,f'Arr.{i+1}:{NAMES[idx].split("(")[0][:9]}',LC[i]))
        cum+=T_OPS
        events.append((cum,f'Handoff Tug{i+1}',TC_[i]))
    cum+=T_TR[seqW[-1],RH_IDX]+T_PH[seqW[-1],RH_IDX]
    events.append((cum,'MS@RH',A2))
    for i in range(5):
        events.append((res['handover'][_wp,i],
                       f'HO:{NAMES[seqW[i]].split("(")[0][:9]}',A4))
    events.append((res['mission_day'][_wp],'✓ DONE',A3))
    events.sort(key=lambda x:x[0])
    for j,(day,label,col) in enumerate(events):
        ax5.axvline(day,color=col,alpha=0.5,lw=0.9)
        yo=0.65 if j%2==0 else 0.22
        ax5.text(day,yo,f'{label}\nd{day:.0f}',ha='center',va='center',
                 fontsize=6.5,color=col,
                 bbox=dict(fc=CB,ec=col,boxstyle='round,pad=0.25',alpha=0.92))
    ax5.axvline(MAX_DAYS,color=A3,lw=1.5,ls='--',alpha=0.85)
    ax5.text(MAX_DAYS+1.5,0.42,'365d',color=A3,fontsize=7.5,va='center')
    ax5.set_xlim(-8,MAX_DAYS+28); ax5.set_ylim(0,1)
    ax5.set_yticks([]); ax5.set_xlabel('Mission day')

    # ── 2. Tug spiral Gantt ───────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[1,0])
    sax(ax4,f'Worst-case MS prop (rank {_wp+1}/{nf}) — tug spiral Gantt')
    for i,idx in enumerate(seqW):
        s=res['tug_start'][_wp,i]; e=res['tug_arrive'][_wp,i]
        ax4.barh(i,e-s,left=s,color=TC_[i],alpha=0.8,ec='none',height=0.55)
        ax4.text(e+1.5,i,f'd{e:.0f}',va='center',fontsize=7,color=TC_[i])
    ax4.axvline(res['ms_return'][_wp],color='white',lw=1.2,ls=':',alpha=0.7,
                label=f"MS@RH d{res['ms_return'][_wp]:.0f}")
    ax4.axvline(MAX_DAYS,color=A3,lw=1.3,ls='--',alpha=0.85,label='365d')
    ax4.set_yticks(range(5))
    ax4.set_yticklabels([NAMES[i].split('(')[0].strip()[:18] for i in seqW],fontsize=7)
    ax4.set_xlabel('Mission day')
    ax4.legend(fontsize=7,facecolor=CB,labelcolor=TC,edgecolor=GR)

    # ── 3. Mothership mass vs time ────────────────────────────────────────
    ax_m = fig.add_subplot(gs[1,1])
    sax(ax_m,f'Worst-case MS prop (rank {_wp+1}/{nf}) — mothership mass vs time')
    tug_mwets0 = np.full(5, TUG_DRY + TUG_MPROP[list(seqW)].max())
    ms_wet0 = MS_DRY + res['ms_prop'][_wp] + tug_mwets0.sum()
    nf_m = [RH_IDX] + list(seqW)
    nt_m = list(seqW) + [RH_IDX]
    m = ms_wet0
    t = 0.0
    t_pts=[t]; m_pts=[m]
    ms_t=[t]; ms_m=[m]; ms_lbl=['RH']
    for i in range(6):
        fi=nf_m[i]; tj=nt_m[i]
        for dv_b in [DV1[fi,tj], DV2[fi,tj], DV_PH[fi,tj]]:
            m_pre=m; m*=np.exp(-dv_b/MS_VEX)
            t_pts+=[t,t]; m_pts+=[m_pre,m]
        t+=T_TR[fi,tj]+T_PH[fi,tj]
        t_pts.append(t); m_pts.append(m)
        lbl=NAMES[tj].split('(')[0][:12] if tj<N_DEB else 'RH'
        ms_t.append(t); ms_m.append(m); ms_lbl.append(lbl)
        if i<5:
            t+=T_OPS; t_pts.append(t); m_pts.append(m)
            m_pre=m; m-=tug_mwets0[i]        # tug detaches with debris after ops
            t_pts+=[t,t]; m_pts+=[m_pre,m]
    ax_m.plot(t_pts,m_pts,color=A1,lw=1.5,zorder=3)
    ax_m.scatter(ms_t,ms_m,color=A2,s=40,zorder=5,lw=0)
    for k,(tt,mm,lb) in enumerate(zip(ms_t,ms_m,ms_lbl)):
        yo=8 if k%2==0 else -12
        ax_m.annotate(lb,(tt,mm),xytext=(0,yo),textcoords='offset points',
                      fontsize=6.5,color=A2,ha='center')
    ax_m.axhline(MS_DRY,color=A3,lw=1.2,ls='--',alpha=0.85,
                 label=f'Dry mass {MS_DRY:.0f} kg')
    ax_m.set_xlabel('Mission day'); ax_m.set_ylabel('Mothership mass [kg]')
    ax_m.legend(fontsize=7,facecolor=CB,labelcolor=TC,edgecolor=GR)

    # ── 4. All-combinations scatter: ΔV vs mission day ───────────────────
    ax6 = fig.add_subplot(gs[1,2])
    sax(ax6, 'All combinations — total ΔV vs mission day')
    sc = ax6.scatter(res['mission_day'], res['total_dv'],
                     c=res['prop_used'], cmap='plasma',
                     s=5, alpha=0.45, lw=0, zorder=2)
    for rank, col, marker, lbl in [
        (0,      A2, '*', f'Best (rank 1)'),
        (nf - 1, A3, '*', f'Worst score (rank {nf})'),
        (_wp,    A5, 'D', f'Worst MS prop (rank {_wp+1})'),
    ]:
        ax6.scatter(res['mission_day'][rank], res['total_dv'][rank],
                    color=col, s=130, marker=marker, zorder=5,
                    edgecolors='white', lw=0.5, label=lbl)
    cb = plt.colorbar(sc, ax=ax6, pad=0.02)
    cb.set_label('Prop used [kg]', color=MU_, fontsize=7)
    cb.ax.yaxis.set_tick_params(color=MU_, labelsize=7)
    plt.setp(cb.ax.yaxis.get_ticklabels(), color=MU_)
    ax6.set_xlabel('Mission day')
    ax6.set_ylabel('Total ΔV [m/s]')
    ax6.legend(fontsize=7, facecolor=CB, labelcolor=TC, edgecolor=GR,
               loc='lower right')

    plt.savefig(save_path,dpi=150,bbox_inches='tight',
                facecolor=BG,edgecolor='none')
    print(f"\n  Plot saved → {save_path}")

# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    print("="*70)
    print("  REAVER TRAJECTORY OPTIMIZER  —  Group 7")
    print("="*70)
    print(f"  MS dry {MS_DRY}kg | MS Isp {MS_ISP}s | prop sized per sequence")
    print(f"  Tug dry {TUG_DRY}kg | Tug Isp {TUG_ISP}s | Thrust {TUG_THR}N")
    print(f"  T_ops {T_OPS}d | Constraint ≤{MAX_DAYS}d | Hohmann phasing {N_PHASE_REV} revs")

    print_tug_analysis()

    t0=time.time()
    res = evaluate_all_sequences()
    print(f"  Total runtime: {time.time()-t0:.2f}s")

    if res:
        nf  = res['n_feasible']
        wp  = int(np.argmax(res['prop_used']))
        wtp = int(np.argmax(res['total_prop']))

        print_top(res, n=1, start=wp,
                  label=f"WORST CASE — MS propellant (rank {wp+1}/{nf})")

        print_prop_comparison(res, wp, wtp)

        print_mass_timeline(res, wp, f"WORST CASE — MS propellant (rank {wp+1}/{nf})")

        # print_target_frequency(res, top_n=res['n_feasible'])
        print_tug_margins(res, wp)
        make_plots(res, wp=wp)

    print("\n  Done.")
