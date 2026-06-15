"""
REAVER — Nominal Mission Evaluator
====================================
Interactive: enter 5 spacecraft IDs in visit order.
Produces the full mission profile (ΔV breakdown, mass timeline, standalone figures).
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import permutations
import time, warnings
warnings.filterwarnings('ignore')

from config import *
from reaver_core import (DV1, DV2, DV_PH, T_TR, T_PH, DV_LEG, T_LEG, T_LEG_FINITE,
                         finite_burn_time, finite_burn_info,
                         TUG_DV, TUG_TIME, TUG_MPROP, TUG_MWET, TUG_WET_LOADED)

# =============================================================================
# USER INPUT
# =============================================================================

def display_catalogue():
    print("\n" + "="*70)
    print("  AVAILABLE SPACECRAFT")
    print("="*70)
    print(f"  {'ID':>4}  {'Name':<35}  {'Mass kg':>8}  {'SMA km':>9}  "
          f"{'Inc °':>6}  {'RAAN °':>7}")
    print(f"  {'':─>4}  {'':─<35}  {'':─>8}  {'':─>9}  {'':─>6}  {'':─>7}")
    for r in _RAW:
        print(f"  {r[0]:>4}  {r[1]:<35}  {r[2]:>8.0f}  {r[3]:>9.3f}  "
              f"{r[4]:>6.4f}  {r[5]:>7.4f}")

def get_debris_set():
    """Ask user for 5 spacecraft IDs (order-free). Returns list of 5 catalogue indices."""
    display_catalogue()
    id_map = {r[0]: i for i, r in enumerate(_RAW)}
    print("\n" + "="*70)
    print("  DEBRIS SELECTION  (order does not matter — optimiser finds best visit order)")
    print("="*70)
    print("  Enter any 5 spacecraft IDs space-separated.")
    print("  Example:  556 628 489 705 834\n")
    while True:
        raw = input("  5 IDs > ").strip().split()
        if len(raw) != 5:
            print(f"  ✗  Need exactly 5 IDs, got {len(raw)}. Try again.")
            continue
        try:
            ids_in = [int(x) for x in raw]
        except ValueError:
            print("  ✗  IDs must be integers. Try again.")
            continue
        bad = [x for x in ids_in if x not in id_map]
        if bad:
            print(f"  ✗  Unknown ID(s): {bad}. Try again.")
            continue
        if len(set(ids_in)) != 5:
            print("  ✗  All 5 IDs must be distinct. Try again.")
            continue
        indices = [id_map[x] for x in ids_in]
        print("\n  Selected debris:")
        for idx in indices:
            print(f"    [{IDS[idx]}]  {NAMES[idx]}  ({MASS[idx]:.0f} kg)")
        confirm = input("\n  Confirm selection? [Y/n] > ").strip().lower()
        if confirm in ('', 'y', 'yes'):
            return indices
        print()


def optimise_ordering(debris_indices):
    """
    Evaluate all 120 orderings of the 5 given debris.
    Returns results sorted by Pareto score (feasible first, then 50% ΔV + 50% time).
    """
    print("\n  Evaluating all 120 orderings...", end=' ', flush=True)
    all_results = []
    for perm in permutations(debris_indices):
        r = evaluate_sequence(list(perm))
        all_results.append(r)
    print("done")

    feasible   = [r for r in all_results if r['feasible']]
    infeasible = [r for r in all_results if not r['feasible']]

    if feasible:
        dvs    = np.array([r['tot_dv']      for r in feasible])
        days   = np.array([r['mission_day'] for r in feasible])
        dv_n   = (dvs  - dvs.min())  / (dvs.max()  - dvs.min()  + 1e-9)
        time_n = (days - days.min()) / (days.max() - days.min() + 1e-9)
        scores = 0.5*dv_n + 0.5*time_n
        order  = np.argsort(scores)
        sorted_results = [feasible[i] for i in order] + infeasible
    else:
        sorted_results = infeasible   # all infeasible — sort by mission day
        sorted_results.sort(key=lambda r: r['mission_day'])

    # Print ranked table
    n_feas = len(feasible)
    print(f"\n  {n_feas}/120 orderings feasible within {MAX_DAYS:.0f} days.\n")
    print(f"  {'Rank':<5} {'Visit order':<60} {'ΔV m/s':>8} {'Day':>7} "
          f"{'Prop kg':>8} {'Feas':>5}")
    print(f"  {'':─<5} {'':─<60} {'':─>8} {'':─>7} {'':─>8} {'':─>5}")
    for rank, r in enumerate(sorted_results[:10]):
        visit = " -> ".join(NAMES[k].split('(')[0].strip()[:9] for k in r['seq'])
        fstr  = '✓' if r['feasible'] else '✗'
        print(f"  {rank+1:<5} {visit:<60} {r['tot_dv']:>8.1f} "
              f"{r['mission_day']:>7.1f} {r['ms_prop']:>8.1f} {fstr:>5}")

    # Let user pick
    print(f"\n  Best ordering is rank 1.  Enter rank to select [1–{min(10,len(sorted_results))}] "
          f"or press Enter to accept rank 1:")
    while True:
        choice = input("  > ").strip()
        if choice == '':
            return sorted_results[0]['seq']
        try:
            n = int(choice)
            if 1 <= n <= min(10, len(sorted_results)):
                return sorted_results[n-1]['seq']
        except ValueError:
            pass
        print(f"  ✗  Enter a number between 1 and {min(10, len(sorted_results))}.")

# =============================================================================
# SINGLE-SEQUENCE EVALUATOR
# =============================================================================

def evaluate_sequence(seq):
    nf = [RH_IDX] + list(seq)
    nt = list(seq) + [RH_IDX]
    tug_prop_uniform = float(TUG_MPROP[list(seq)].max())   # worst-case tug sets the standard
    tug_mwets = np.full(5, TUG_DRY + tug_prop_uniform)

    # Backward pass: orbital prop sized so final mass = MS_DRY (transfers only)
    m = MS_DRY
    m = m * np.exp(DV_LEG[seq[-1], RH_IDX] / MS_VEX)
    for i in range(3, -1, -1):
        m += tug_mwets[i+1]
        m  = m * np.exp(DV_LEG[seq[i], seq[i+1]] / MS_VEX)
    m += tug_mwets[0]
    m  = m * np.exp(DV_LEG[RH_IDX, seq[0]] / MS_VEX)
    ms_prop_orbital = m - MS_DRY - tug_mwets.sum()

    # 10 % RCS margin added at launch — covers all proximity operations (RPO)
    rcs_allocation = MS_RCS_MARGIN * ms_prop_orbital
    ms_prop        = ms_prop_orbital + rcs_allocation   # total initial propellant
    ms_wet0        = MS_DRY + tug_mwets.sum() + ms_prop

    # Forward pass: mass + timing (orbital burns + RPO burns interleaved)
    dv_legs = np.array([DV_LEG[nf[i], nt[i]] for i in range(6)])
    t_legs  = np.array([T_LEG_FINITE[nf[i], nt[i]] for i in range(6)])

    mass = ms_wet0
    mass_vec = []
    t = 0.0
    tug_starts = []
    for i in range(6):
        mass = mass * np.exp(-dv_legs[i] / MS_VEX)
        mass_vec.append(mass)
        t += t_legs[i]
        if i < 5:
            # RPO at debris: V-bar inspect + tumble-axis capture + detumble
            mass = mass * np.exp(-(DV_RPO_DEBRIS + DV_RPO_DETUMBLE) / MS_VEX)
            t += T_OPS
            tug_starts.append(t)
            mass -= tug_mwets[i]
    # 5 × RH proximity docking (tug-meet + dock) as each tug+debris arrives
    for _ in range(5):
        mass = mass * np.exp(-DV_RPO_RH / MS_VEX)
    ms_return = t

    rcs_margin  = mass - MS_DRY   # remaining RCS propellant at end of mission
    tug_arrive  = np.array([tug_starts[i] + TUG_TIME[seq[i]] for i in range(5)])

    # Sequential RH handovers: each must wait for (a) the previous handover to finish
    # and (b) the tug to arrive — no two handovers overlap.
    arrive_order = np.argsort(tug_arrive)
    handover = np.zeros(5)
    t_rh = ms_return
    for rank in range(5):
        i = int(arrive_order[rank])
        t_start = max(t_rh, float(tug_arrive[i]))
        handover[i] = t_start + T_OPS
        t_rh = handover[i]
    mission_day = handover.max()
    tot_dv      = dv_legs.sum()

    return dict(
        seq             = seq,
        ms_prop_orbital = ms_prop_orbital,
        rcs_allocation  = rcs_allocation,
        ms_prop         = ms_prop,         # total initial prop (orbital + RCS)
        ms_wet0         = ms_wet0,
        tug_mwets       = tug_mwets,
        dv_legs         = dv_legs,
        t_legs          = t_legs,
        mass_vec        = np.array(mass_vec),
        tug_starts      = np.array(tug_starts),
        tug_arrive      = tug_arrive,
        handover        = handover,
        ms_return       = ms_return,
        mission_day     = mission_day,
        tot_dv          = tot_dv,
        feasible        = mission_day <= MAX_DAYS,
        rcs_margin      = rcs_margin,
    )

# =============================================================================
# PRINT RESULTS
# =============================================================================

def print_results(r):
    seq = r['seq']
    print("\n" + "="*72)
    print("  NOMINAL MISSION PROFILE")
    print("="*72)

    # Visit order
    print("\n  Visit order:")
    for i, idx in enumerate(seq):
        flag = '  ⚠ >2000 kg' if MASS[idx] > SOFT_MASS else ''
        print(f"    {i+1}. [{IDS[idx]}]  {NAMES[idx]:<35}  {MASS[idx]:>7.0f} kg{flag}")

    # Leg ΔV table
    leg_names = [NAMES[j] for j in seq] + ['Recycling Hub']
    nf = [RH_IDX] + list(seq); nt = list(seq) + [RH_IDX]
    print(f"\n  {'Leg':<3} {'To':<28} {'ΔV tot':>8} {'ΔV Hohm':>8} "
          f"{'ΔV plane':>9} {'ΔV ph':>7} {'Transf d':>9} {'Phase d':>8}")
    print(f"  {'':─<3} {'':─<28} {'':─>8} {'':─>8} "
          f"{'':─>9} {'':─>7} {'':─>9} {'':─>8}")
    for li in range(6):
        fi=nf[li]; ti=nt[li]; nm=leg_names[li][:28]
        dv_tot = r['dv_legs'][li]
        print(f"  {str(li+1) if li<5 else 'R':<3} {nm:<28} "
              f"{dv_tot:>7.1f}m {DV1[fi,ti]:>7.1f}m "
              f"{DV2[fi,ti]:>8.1f}m {DV_PH[fi,ti]:>6.1f}m "
              f"{T_TR[fi,ti]:>8.1f}d {T_PH[fi,ti]:>7.1f}d")

    # Tug table
    print(f"\n  {'Tug':<3} {'Debris':<28} {'ΔV m/s':>8} {'Prop kg':>8} "
          f"{'Spiral d':>9} {'Start d':>8} {'Arrive d':>9} {'Handover':>9}")
    print(f"  {'':─<3} {'':─<28} {'':─>8} {'':─>8} "
          f"{'':─>9} {'':─>8} {'':─>9} {'':─>9}")
    tug_prop_uniform = float(r['tug_mwets'][0] - TUG_DRY)
    for i, idx in enumerate(seq):
        worst_flag = '  ← sizing driver' if TUG_MPROP[idx] == TUG_MPROP[list(seq)].max() else ''
        print(f"  {i+1:<3} {NAMES[idx]:<28} "
              f"{TUG_DV[idx]:>7.1f}m {tug_prop_uniform:>7.1f}kg "
              f"{TUG_TIME[idx]:>8.1f}d {r['tug_starts'][i]:>7.1f}d "
              f"{r['tug_arrive'][i]:>8.1f}d {r['handover'][i]:>8.1f}d{worst_flag}")

    # Summary box
    feas_str = "✓ FEASIBLE" if r['feasible'] else "✗ EXCEEDS 365 d"
    tug_prop_total = float((r['tug_mwets'] - TUG_DRY).sum())
    rcs_pct = 100.0 * r['rcs_margin'] / r['rcs_allocation']
    print(f"\n  ┌{'─'*60}┐")
    print(f"  │  Mothership ΔV total        : {r['tot_dv']:>8.1f} m/s                   │")
    print(f"  │  MS orbital propellant      : {r['ms_prop_orbital']:>8.1f} kg                    │")
    print(f"  │  RCS margin (+10% budget)   : {r['rcs_allocation']:>8.1f} kg                    │")
    print(f"  │  MS total initial prop      : {r['ms_prop']:>8.1f} kg                    │")
    print(f"  │  Tug propellant total       : {tug_prop_total:>8.1f} kg                    │")
    print(f"  │  Combined propellant        : {r['ms_prop']+tug_prop_total:>8.1f} kg                    │")
    print(f"  │  MS returns to RH day       : {r['ms_return']:>8.1f} d                     │")
    print(f"  │  Total debris mass removed  : {MASS[list(seq)].sum():>8.0f} kg                    │")
    print(f"  │  RCS margin remaining       : {r['rcs_margin']:>8.1f} kg  ({rcs_pct:5.1f}% of budget)  │")
    print(f"  │  Mission completion         : {r['mission_day']:>8.1f} d  {feas_str:<14}     │")
    print(f"  └{'─'*60}┘")

# =============================================================================
# MASS TIMELINE PRINT
# =============================================================================

def print_mass_timeline(r):
    seq = r['seq']
    nf_ = [RH_IDX] + list(seq); nt_ = list(seq) + [RH_IDX]
    m0  = r['ms_wet0']
    tug_mwets = r['tug_mwets']
    W = 64

    print("\n" + "="*98)
    print("  MOTHERSHIP MASS TIMELINE")
    print("="*98)
    print("  " + "  ->  ".join(NAMES[k].split('(')[0].strip()[:14] for k in seq))
    print(f"\n  Initial: {MS_DRY:.0f} kg dry  +  {r['ms_prop_orbital']:.1f} kg orbital prop"
          f"  +  {r['rcs_allocation']:.1f} kg RCS (+10%)  +  {tug_mwets.sum():.1f} kg tugs  =  {m0:.1f} kg")
    print(f"\n  {'Day':>7}  {'Event':<{W}}  {'D Mass kg':>10}  {'Mass kg':>9}")
    print(f"  {'':->7}  {'':─<{W}}  {'':->10}  {'':->9}")

    m = m0; t = 0.0

    def prow(day, event, delta=None):
        nonlocal m
        ds = "—"
        if delta is not None:
            m += delta; ds = f"{delta:+.1f}"
        print(f"  {day:>7.1f}  {event:<{W}}  {ds:>10}  {m:>9.1f}")

    prow(t, "Depart Recycling Hub")
    for i in range(6):
        fi=nf_[i]; ti=nt_[i]
        dest = NAMES[ti].split('(')[0].strip()[:28] if ti < N_DEB else 'Recycling Hub'
        asc  = SMA_ALL[ti] >= SMA_ALL[fi]
        d1,d2,dph = DV1[fi,ti], DV2[fi,ti], DV_PH[fi,ti]

        m_b1 = m
        dm = -m*(1 - np.exp(-d1/MS_VEX))
        nf1, no1 = finite_burn_info(d1, m_b1)
        b1 = (f"  Burn 1  Hohmann dep       ({d1:6.1f} m/s | {nf1} fires, {no1} orbits)" if asc else
              f"  Burn 1  dep + plane chg   ({d1:6.1f} m/s | {nf1} fires, {no1} orbits)")
        prow(t, b1, dm)
        t += finite_burn_time(d1, m_b1, SMA_ALL[fi])

        t += T_TR[fi,ti]

        m_b2 = m
        dm = -m*(1 - np.exp(-d2/MS_VEX))
        nf2, no2 = finite_burn_info(d2, m_b2)
        b2 = (f"  Burn 2  circ + plane chg  ({d2:6.1f} m/s | {nf2} fires, {no2} orbits)" if asc else
              f"  Burn 2  circularise        ({d2:6.1f} m/s | {nf2} fires, {no2} orbits)")
        prow(t, b2, dm)
        t += finite_burn_time(d2, m_b2, SMA_ALL[ti])

        m_b3 = m
        dm_ph = -m*(1 - np.exp(-dph/MS_VEX))
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
            arrive_order_rh = np.argsort(r['tug_arrive'])
            t_rh = t
            for j_rank in range(5):
                j = int(arrive_order_rh[j_rank])
                ta = float(r['tug_arrive'][j])
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
    rcs_pct = 100.0 * rcs_remaining / r['rcs_allocation']
    print(f"\n  Orbital prop    : {r['ms_prop_orbital']:.1f} kg"
          f"  |  RCS budget (+10%) : {r['rcs_allocation']:.1f} kg"
          f"  |  Total init prop : {r['ms_prop']:.1f} kg")
    print(f"  RCS remaining   : {rcs_remaining:.1f} kg  ({rcs_pct:.1f}% of budget)"
          f"  |  Mission day : {r['mission_day']:.1f} d")

# =============================================================================
# TUG PROPELLANT MARGINS
# =============================================================================

def print_tug_margins(r):
    """For each tug in the sequence, show propellant required vs uniform budget."""
    seq              = r['seq']
    tug_prop_uniform = float(r['tug_mwets'][0] - TUG_DRY)

    print("\n" + "="*76)
    print(f"  TUG PROPELLANT MARGINS  (uniform budget: {tug_prop_uniform:.1f} kg/tug)")
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
# STANDALONE FIGURES
# =============================================================================

def make_plot(r, save_dir=None):
    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'figures', 'nominal_mission')
    os.makedirs(save_dir, exist_ok=True)

    seq  = r['seq']
    nf_  = [RH_IDX] + list(seq)
    nt_  = list(seq) + [RH_IDX]
    cyc  = plt.rcParams['axes.prop_cycle'].by_key()['color']
    saved = []

    # ── 1. Full mission timeline ──────────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(14, 4))
    events = [(0, 'Depart RH', 'gray')]
    cum = 0.0
    for i, idx in enumerate(seq):
        cum += T_TR[nf_[i], nt_[i]] + T_PH[nf_[i], nt_[i]]
        events.append((cum, f'Arr.{i+1}:{NAMES[idx].split("(")[0][:9]}', cyc[i % len(cyc)]))
        cum += T_OPS
        events.append((cum, f'Handoff Tug{i+1}', cyc[i % len(cyc)]))
    cum += T_TR[seq[-1], RH_IDX] + T_PH[seq[-1], RH_IDX]
    events.append((cum, 'MS@RH', 'green'))
    for i in range(5):
        events.append((r['handover'][i], f'HO:{NAMES[seq[i]].split("(")[0][:9]}', 'purple'))
    events.append((r['mission_day'], 'Mission Complete', 'red'))
    events.sort(key=lambda x: x[0])
    for j, (day, label, col) in enumerate(events):
        ax1.axvline(day, color=col, alpha=0.5, lw=0.9)
        yo = 0.65 if j % 2 == 0 else 0.22
        ax1.text(day, yo, f'{label}\nd{day:.0f}', ha='center', va='center',
                 fontsize=7, color=col,
                 bbox=dict(fc='white', ec=col, boxstyle='round,pad=0.25', alpha=0.9))
    ax1.axvline(MAX_DAYS, color='red', lw=1.5, ls='--', alpha=0.85)
    ax1.text(MAX_DAYS + 1.5, 0.42, '365d', color='red', fontsize=8, va='center')
    ax1.set_xlim(-8, max(MAX_DAYS + 28, r['mission_day'] + 20))
    ax1.set_ylim(0, 1)
    ax1.set_yticks([])
    ax1.set_xlabel('Mission day')
    ax1.set_title('Full Mission Timeline')
    ax1.grid(True, alpha=0.4)
    fig1.tight_layout()
    p1 = os.path.join(save_dir, '01_mission_timeline.png')
    fig1.savefig(p1, dpi=150, bbox_inches='tight')
    plt.close(fig1)
    saved.append(p1)

    # ── 2. Tug spiral Gantt ───────────────────────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(10, 5))
    for i, idx in enumerate(seq):
        s = r['tug_starts'][i]; e = r['tug_arrive'][i]
        ax2.barh(i, e - s, left=s, alpha=0.85, height=0.55, color=cyc[i % len(cyc)])
        ax2.text(e + 1.5, i, f'd{e:.0f}', va='center', fontsize=8)
    ax2.axvline(r['ms_return'], color='gray', lw=1.2, ls=':', alpha=0.8,
                label=f"MS@RH d{r['ms_return']:.0f}")
    ax2.axvline(MAX_DAYS, color='red', lw=1.3, ls='--', alpha=0.85, label='365d limit')
    ax2.set_yticks(range(5))
    ax2.set_yticklabels([NAMES[i].split('(')[0].strip()[:18] for i in seq], fontsize=9)
    ax2.set_xlabel('Mission day')
    ax2.set_title('Tug Spiral Gantt')
    ax2.legend()
    ax2.grid(True, alpha=0.4, axis='x')
    fig2.tight_layout()
    p2 = os.path.join(save_dir, '02_tug_spiral_gantt.png')
    fig2.savefig(p2, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    saved.append(p2)

    # ── 3. Mothership mass vs time ────────────────────────────────────────────
    fig3, ax3 = plt.subplots(figsize=(10, 5))
    tug_mwets = r['tug_mwets']
    m = r['ms_wet0']; t = 0.0
    t_pts = [t]; m_pts = [m]; ms_t = [t]; ms_m = [m]; ms_lbl = ['RH']
    for i in range(6):
        fi = nf_[i]; tj = nt_[i]
        for dv_b in [DV1[fi, tj], DV2[fi, tj], DV_PH[fi, tj]]:
            m_pre = m; m *= np.exp(-dv_b / MS_VEX)
            t_pts += [t, t]; m_pts += [m_pre, m]
        t += T_TR[fi, tj] + T_PH[fi, tj]
        t_pts.append(t); m_pts.append(m)
        lbl = NAMES[tj].split('(')[0][:12] if tj < N_DEB else 'RH'
        ms_t.append(t); ms_m.append(m); ms_lbl.append(lbl)
        if i < 5:
            for dv_rpo in [DV_RPO_DEBRIS, DV_RPO_DETUMBLE]:
                m_pre = m; m *= np.exp(-dv_rpo / MS_VEX)
                t_pts += [t, t]; m_pts += [m_pre, m]
            t += T_OPS; t_pts.append(t); m_pts.append(m)
            m_pre = m; m -= tug_mwets[i]
            t_pts += [t, t]; m_pts += [m_pre, m]
        else:
            for _ in range(5):
                m_pre = m; m *= np.exp(-DV_RPO_RH / MS_VEX)
                t_pts += [t, t]; m_pts += [m_pre, m]
    ax3.plot(t_pts, m_pts, lw=1.5, zorder=3, label='MS mass')
    ax3.scatter(ms_t, ms_m, s=50, zorder=5, label='Key event')
    for k, (tt, mm, lb) in enumerate(zip(ms_t, ms_m, ms_lbl)):
        yo = 8 if k % 2 == 0 else -12
        ax3.annotate(lb, (tt, mm), xytext=(0, yo), textcoords='offset points',
                     fontsize=7, ha='center')
    ax3.axhline(MS_DRY, color='red', lw=1.2, ls='--', alpha=0.85,
                label=f'Dry mass {MS_DRY:.0f} kg')
    ax3.axhline(MS_DRY + r['rcs_allocation'], color='orange', lw=1.0, ls=':',
                alpha=0.75, label=f'Dry+RCS budget {MS_DRY + r["rcs_allocation"]:.0f} kg')
    ax3.set_xlabel('Mission day')
    ax3.set_ylabel('Mothership mass [kg]')
    ax3.set_title('Mothership Mass vs Time')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.4)
    fig3.tight_layout()
    p3 = os.path.join(save_dir, '03_mothership_mass_vs_time.png')
    fig3.savefig(p3, dpi=150, bbox_inches='tight')
    plt.close(fig3)
    saved.append(p3)

    # ── 4. Per-leg ΔV stacked bar ─────────────────────────────────────────────
    fig4, ax4 = plt.subplots(figsize=(9, 5))
    components = ['Hohmann dep', 'Circ+plane', 'Phasing']
    dvs = np.array([[DV1[nf_[li], nt_[li]], DV2[nf_[li], nt_[li]], DV_PH[nf_[li], nt_[li]]]
                    for li in range(6)])
    bottoms = np.zeros(6)
    bar_colors = [cyc[0], cyc[3 % len(cyc)], cyc[4 % len(cyc)]]
    for ci, (comp, col) in enumerate(zip(components, bar_colors)):
        vals = dvs[:, ci]
        ax4.bar(range(6), vals, bottom=bottoms, alpha=0.85, label=comp, color=col)
        for li in range(6):
            if vals[li] > 15:
                ax4.text(li, bottoms[li] + vals[li] / 2, f'{vals[li]:.0f}',
                         ha='center', va='center', fontsize=7, color='white', fontweight='bold')
        bottoms += vals
    xlbls = [f'L{i+1}\n{NAMES[seq[i]].split("(")[0][:8]}' for i in range(5)] + ['Ret\nRH']
    ax4.set_xticks(range(6))
    ax4.set_xticklabels(xlbls, fontsize=8)
    ax4.set_ylabel('ΔV [m/s]')
    ax4.set_title('ΔV Breakdown per Leg')
    total_per_leg = dvs.sum(axis=1)
    for li in range(6):
        ax4.text(li, bottoms[li] + 5, f'{total_per_leg[li]:.0f}',
                 ha='center', va='bottom', fontsize=8)
    ax4.legend()
    ax4.grid(True, alpha=0.4, axis='y')
    fig4.tight_layout()
    p4 = os.path.join(save_dir, '04_dv_breakdown_per_leg.png')
    fig4.savefig(p4, dpi=150, bbox_inches='tight')
    plt.close(fig4)
    saved.append(p4)

    for p in saved:
        print(f"  Plot saved → {p}")
    return saved

# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    print("="*70)
    print("  REAVER — Nominal Mission Evaluator")
    print("="*70)
    print(f"  MS dry {MS_DRY} kg | Isp {MS_ISP} s  |  "
          f"Tug dry {TUG_DRY} kg | Isp {TUG_ISP} s | Thrust {TUG_THR} N")
    print(f"  RH: SMA {RH_SMA/1e3:.0f} km | Inc {RH_INC}° | RAAN {RH_RAAN}°"
          f"  |  T_ops {T_OPS} d  |  Limit {MAX_DAYS:.0f} d")

    debris  = get_debris_set()
    seq     = optimise_ordering(debris)
    r       = evaluate_sequence(seq)

    print_results(r)
    print_mass_timeline(r)
    print_tug_margins(r)
    make_plot(r)

    if not r['feasible']:
        print(f"\n  ⚠  Mission exceeds {MAX_DAYS:.0f}-day constraint "
              f"({r['mission_day']:.1f} d).")
    print("\n  Done.")