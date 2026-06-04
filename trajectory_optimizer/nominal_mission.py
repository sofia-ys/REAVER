"""
REAVER — Nominal Mission Evaluator
====================================
Interactive: enter 5 spacecraft IDs in visit order.
Produces the full mission profile (ΔV breakdown, mass timeline, dashboard plot).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from itertools import permutations
import time, warnings
warnings.filterwarnings('ignore')

from config import *
from reaver_core import (DV1, DV2, DV_PH, T_TR, T_PH, DV_LEG, T_LEG,
                         TUG_DV, TUG_TIME, TUG_MPROP, TUG_MWET, TUG_WET_LOADED,
                         T_PH_RH)

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

    # Backward pass: propellant sized to end at MS_DRY exactly
    m = MS_DRY
    m = m * np.exp(DV_LEG[seq[-1], RH_IDX] / MS_VEX)
    for i in range(3, -1, -1):
        m += tug_mwets[i+1]
        m  = m * np.exp(DV_LEG[seq[i], seq[i+1]] / MS_VEX)
    m += tug_mwets[0]
    m  = m * np.exp(DV_LEG[RH_IDX, seq[0]] / MS_VEX)
    ms_prop = m - MS_DRY - tug_mwets.sum()
    ms_wet0 = m

    # Forward pass: mass + timing
    dv_legs = np.array([DV_LEG[nf[i], nt[i]] for i in range(6)])
    t_legs  = np.array([T_LEG[nf[i],  nt[i]] for i in range(6)])

    mass = ms_wet0
    mass_vec = []
    t = 0.0
    tug_starts = []
    for i in range(6):
        mass = mass * np.exp(-dv_legs[i] / MS_VEX)
        mass_vec.append(mass)
        t += t_legs[i]
        if i < 5:
            t += T_OPS
            tug_starts.append(t)
            mass -= tug_mwets[i]
    ms_return = t

    tug_arrive  = np.array([tug_starts[i] + TUG_TIME[seq[i]] for i in range(5)])
    handover    = np.array([max(ms_return, tug_arrive[i]) + T_PH_RH for i in range(5)])
    mission_day = handover.max()
    tot_dv      = dv_legs.sum()

    return dict(
        seq         = seq,
        ms_prop     = ms_prop,
        ms_wet0     = ms_wet0,
        tug_mwets   = tug_mwets,
        dv_legs     = dv_legs,
        t_legs      = t_legs,
        mass_vec    = np.array(mass_vec),
        tug_starts  = np.array(tug_starts),
        tug_arrive  = tug_arrive,
        handover    = handover,
        ms_return   = ms_return,
        mission_day = mission_day,
        tot_dv      = tot_dv,
        feasible    = mission_day <= MAX_DAYS,
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
    print(f"\n  ┌{'─'*58}┐")
    print(f"  │  Mothership ΔV total      : {r['tot_dv']:>8.1f} m/s                 │")
    print(f"  │  MS propellant used       : {r['ms_prop']:>8.1f} kg                  │")
    print(f"  │  Tug propellant total     : {tug_prop_total:>8.1f} kg                  │")
    print(f"  │  Combined propellant      : {r['ms_prop']+tug_prop_total:>8.1f} kg                  │")
    print(f"  │  MS returns to RH day     : {r['ms_return']:>8.1f} d                   │")
    print(f"  │  Total debris mass removed: {MASS[list(seq)].sum():>8.0f} kg                  │")
    print(f"  │  Mission completion       : {r['mission_day']:>8.1f} d  {feas_str:<14}   │")
    print(f"  └{'─'*58}┘")

# =============================================================================
# MASS TIMELINE PRINT
# =============================================================================

def print_mass_timeline(r):
    seq = r['seq']
    nf_ = [RH_IDX] + list(seq); nt_ = list(seq) + [RH_IDX]
    m0  = r['ms_wet0']
    tug_mwets = r['tug_mwets']
    W = 42

    print("\n" + "="*74)
    print("  MOTHERSHIP MASS TIMELINE")
    print("="*74)
    print("  " + "  ->  ".join(NAMES[k].split('(')[0].strip()[:14] for k in seq))
    print(f"\n  Initial: {MS_DRY:.0f} kg dry  +  {r['ms_prop']:.1f} kg prop"
          f"  +  {tug_mwets.sum():.1f} kg tugs  =  {m0:.1f} kg")
    print(f"\n  {'Day':>7}  {'Event':<{W}}  {'D Mass kg':>10}  {'Mass kg':>9}")
    print(f"  {'':->7}  {'':−<{W}}  {'':->10}  {'':->9}")

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
        dest = NAMES[ti].split('(')[0].strip()[:26] if ti < N_DEB else 'Recycling Hub'
        asc  = SMA_ALL[ti] >= SMA_ALL[fi]
        d1,d2,dph = DV1[fi,ti], DV2[fi,ti], DV_PH[fi,ti]
        dm = -m*(1 - np.exp(-d1/MS_VEX))
        b1 = (f"  Burn 1  Hohmann dep      ({d1:6.1f} m/s)" if asc else
              f"  Burn 1  dep + plane chg  ({d1:6.1f} m/s)")
        prow(t, b1, dm)
        t += T_TR[fi,ti]
        dm = -m*(1 - np.exp(-d2/MS_VEX))
        b2 = (f"  Burn 2  circ + plane chg ({d2:6.1f} m/s)" if asc else
              f"  Burn 2  circularise      ({d2:6.1f} m/s)")
        prow(t, b2, dm)
        dm = -m*(1 - np.exp(-dph/MS_VEX))
        prow(t, f"  Burn 3  phasing          ({dph:6.1f} m/s)", dm)
        t += T_PH[fi,ti]
        if i < 5:
            prow(t, f"  Rendezvous: {dest}", None)
            t += T_OPS
            prow(t, f"  Tug {i+1} released ({tug_mwets[i]:.1f} kg)", -tug_mwets[i])
        else:
            prow(t, "  Arrive Recycling Hub", None)

    print(f"  {'':->7}  {'':−<{W}}  {'':->10}  {'':->9}")
    print(f"\n  Prop required : {r['ms_prop']:.1f} kg"
          f"  |  Final margin : {m - MS_DRY:.1f} kg"
          f"  |  Mission day : {r['mission_day']:.1f} d")

# =============================================================================
# DASHBOARD PLOT
# =============================================================================

def make_plot(r, save_path=r'C:\Projects\DSE\REAVER\trajectory_optimizer\nominal_mission.png'):
    BG='#0d1117'; CB='#161b22'; TC='#c9d1d9'; MU_='#8b949e'; GR='#21262d'
    A1='#58a6ff'; A2='#3fb950'; A3='#f78166'; A4='#d2a8ff'; A5='#ffa657'
    LC=[A1,A2,A4,A3,A5,'#79c0ff']; TC_=[A1,A2,A4,A3,A5]

    seq = r['seq']
    fig = plt.figure(figsize=(20,11)); fig.patch.set_facecolor(BG)
    gs  = GridSpec(2,3,figure=fig,hspace=0.52,wspace=0.36,height_ratios=[1,1.3])

    def sax(ax, title=''):
        ax.set_facecolor(CB)
        [s.set_edgecolor(GR) for s in ax.spines.values()]
        ax.tick_params(colors=MU_,labelsize=8)
        ax.xaxis.label.set_color(MU_); ax.yaxis.label.set_color(MU_)
        ax.grid(True,color=GR,lw=0.5,alpha=0.7)
        if title: ax.set_title(title,color=TC,fontsize=9,fontweight='bold',pad=8)

    seq_str = "  →  ".join(f"[{IDS[k]}] {NAMES[k].split('(')[0].strip()[:12]}" for k in seq)
    feas = "✓ FEASIBLE" if r['feasible'] else "✗ INFEASIBLE"
    fig.text(0.5,0.980,'REAVER — Nominal Mission Dashboard',
             ha='center',color=TC,fontsize=15,fontweight='bold')
    fig.text(0.5,0.958, seq_str, ha='center',color=MU_,fontsize=8)
    fig.text(0.5,0.940,
             f"MS prop {r['ms_prop']:.0f} kg  |  Total ΔV {r['tot_dv']:.0f} m/s  |  "
             f"Mission day {r['mission_day']:.1f}  |  {feas}",
             ha='center',color=A2 if r['feasible'] else A3,fontsize=8)

    # ── 1. Full mission timeline ──────────────────────────────────────────────
    ax5 = fig.add_subplot(gs[0,:])
    sax(ax5,'Full mission timeline')
    events = [(0,'Depart RH',MU_)]
    cum = 0.0
    nf_ = [RH_IDX]+list(seq); nt_ = list(seq)+[RH_IDX]
    for i,idx in enumerate(seq):
        cum += T_TR[nf_[i],nt_[i]] + T_PH[nf_[i],nt_[i]]
        events.append((cum,f'Arr.{i+1}:{NAMES[idx].split("(")[0][:9]}',LC[i]))
        cum += T_OPS
        events.append((cum,f'Handoff Tug{i+1}',TC_[i]))
    cum += T_TR[seq[-1],RH_IDX] + T_PH[seq[-1],RH_IDX]
    events.append((cum,'MS@RH',A2))
    for i in range(5):
        events.append((r['handover'][i],f'HO:{NAMES[seq[i]].split("(")[0][:9]}',A4))
    events.append((r['mission_day'],'✓ DONE',A3))
    events.sort(key=lambda x:x[0])
    for j,(day,label,col) in enumerate(events):
        ax5.axvline(day,color=col,alpha=0.5,lw=0.9)
        yo = 0.65 if j%2==0 else 0.22
        ax5.text(day,yo,f'{label}\nd{day:.0f}',ha='center',va='center',
                 fontsize=6.5,color=col,
                 bbox=dict(fc=CB,ec=col,boxstyle='round,pad=0.25',alpha=0.92))
    ax5.axvline(MAX_DAYS,color=A3,lw=1.5,ls='--',alpha=0.85)
    ax5.text(MAX_DAYS+1.5,0.42,'365d',color=A3,fontsize=7.5,va='center')
    ax5.set_xlim(-8,max(MAX_DAYS+28, r['mission_day']+20))
    ax5.set_ylim(0,1); ax5.set_yticks([]); ax5.set_xlabel('Mission day')

    # ── 2. Tug Gantt ─────────────────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[1,0])
    sax(ax4,'Tug spiral Gantt')
    for i,idx in enumerate(seq):
        s=r['tug_starts'][i]; e=r['tug_arrive'][i]
        ax4.barh(i,e-s,left=s,color=TC_[i],alpha=0.8,ec='none',height=0.55)
        ax4.text(e+1.5,i,f'd{e:.0f}',va='center',fontsize=7,color=TC_[i])
    ax4.axvline(r['ms_return'],color='white',lw=1.2,ls=':',alpha=0.7,
                label=f"MS@RH d{r['ms_return']:.0f}")
    ax4.axvline(MAX_DAYS,color=A3,lw=1.3,ls='--',alpha=0.85,label='365d')
    ax4.set_yticks(range(5))
    ax4.set_yticklabels([NAMES[i].split('(')[0].strip()[:18] for i in seq],fontsize=7)
    ax4.set_xlabel('Mission day')
    ax4.legend(fontsize=7,facecolor=CB,labelcolor=TC,edgecolor=GR)

    # ── 3. Mothership mass vs time ────────────────────────────────────────────
    ax_m = fig.add_subplot(gs[1,1])
    sax(ax_m,'Mothership mass vs time')
    tug_mwets = r['tug_mwets']
    m = r['ms_wet0']; t = 0.0
    t_pts=[t]; m_pts=[m]; ms_t=[t]; ms_m=[m]; ms_lbl=['RH']
    for i in range(6):
        fi=nf_[i]; tj=nt_[i]
        for dv_b in [DV1[fi,tj],DV2[fi,tj],DV_PH[fi,tj]]:
            m_pre=m; m*=np.exp(-dv_b/MS_VEX)
            t_pts+=[t,t]; m_pts+=[m_pre,m]
        t += T_TR[fi,tj]+T_PH[fi,tj]
        t_pts.append(t); m_pts.append(m)
        lbl = NAMES[tj].split('(')[0][:12] if tj<N_DEB else 'RH'
        ms_t.append(t); ms_m.append(m); ms_lbl.append(lbl)
        if i < 5:
            t += T_OPS; t_pts.append(t); m_pts.append(m)
            m_pre=m; m -= tug_mwets[i]
            t_pts+=[t,t]; m_pts+=[m_pre,m]
    ax_m.plot(t_pts,m_pts,color=A1,lw=1.5,zorder=3)
    ax_m.scatter(ms_t,ms_m,color=A2,s=40,zorder=5,lw=0)
    for k,(tt,mm,lb) in enumerate(zip(ms_t,ms_m,ms_lbl)):
        yo = 8 if k%2==0 else -12
        ax_m.annotate(lb,(tt,mm),xytext=(0,yo),textcoords='offset points',
                      fontsize=6.5,color=A2,ha='center')
    ax_m.axhline(MS_DRY,color=A3,lw=1.2,ls='--',alpha=0.85,
                 label=f'Dry mass {MS_DRY:.0f} kg')
    ax_m.set_xlabel('Mission day'); ax_m.set_ylabel('Mothership mass [kg]')
    ax_m.legend(fontsize=7,facecolor=CB,labelcolor=TC,edgecolor=GR)

    # ── 4. Per-leg ΔV stacked bar ─────────────────────────────────────────────
    ax6 = fig.add_subplot(gs[1,2])
    sax(ax6,'ΔV breakdown per leg')
    components = ['Hohmann dep', 'Circ+plane', 'Phasing']
    colors_dv  = [A1, A4, A5]
    bottoms = np.zeros(6)
    dvs = np.array([[DV1[nf_[li],nt_[li]], DV2[nf_[li],nt_[li]], DV_PH[nf_[li],nt_[li]]]
                    for li in range(6)])
    for ci, (comp, col) in enumerate(zip(components, colors_dv)):
        vals = dvs[:, ci]
        bars = ax6.bar(range(6), vals, bottom=bottoms, color=col, alpha=0.85,
                       ec='none', label=comp)
        for li in range(6):
            if vals[li] > 15:
                ax6.text(li, bottoms[li]+vals[li]/2, f'{vals[li]:.0f}',
                         ha='center',va='center',fontsize=6.5,color='white',fontweight='bold')
        bottoms += vals
    xlbls = [f'L{i+1}\n{NAMES[seq[i]].split("(")[0][:8]}' for i in range(5)] + ['Ret\nRH']
    ax6.set_xticks(range(6)); ax6.set_xticklabels(xlbls, fontsize=6.5)
    ax6.set_ylabel('ΔV [m/s]')
    total_per_leg = dvs.sum(axis=1)
    for li in range(6):
        ax6.text(li, bottoms[li]+5, f'{total_per_leg[li]:.0f}',
                 ha='center',va='bottom',fontsize=7,color=TC)
    ax6.legend(fontsize=7,facecolor=CB,labelcolor=TC,edgecolor=GR)

    plt.savefig(save_path,dpi=150,bbox_inches='tight',facecolor=BG,edgecolor='none')
    print(f"\n  Plot saved → {save_path}")

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
    make_plot(r)

    if not r['feasible']:
        print(f"\n  ⚠  Mission exceeds {MAX_DAYS:.0f}-day constraint "
              f"({r['mission_day']:.1f} d).")
    print("\n  Done.")