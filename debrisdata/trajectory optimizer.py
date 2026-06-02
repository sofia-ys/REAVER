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
  Burn 3: Drift phasing (90deg avg, delta_a/a=0.001)

Tug model: Edelbaum low-thrust optimal transfer (iterative mass sizing)
"""

import numpy as np
from itertools import combinations, permutations
import matplotlib

from target_filter import parameter_table

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from collections import Counter
import sys, time, warnings
warnings.filterwarnings('ignore')

# ── Constants ──────────────────────────────────────────────────────────────
MU   = 3.986004418e14   # m³/s²
G0   = 9.80665          # m/s²
DAY  = 86400.0
D2R  = np.pi / 180.0

# ── Spacecraft ─────────────────────────────────────────────────────────────
MS_DRY   = 1300.0       # kg
MS_ISP   = 253.0        # s  (monoprop chemical)
MS_VEX   = MS_ISP * G0  # m/s

TUG_DRY  = 200.0        # kg each
TUG_ISP  = 1850.0       # s  (E-propulsion)
TUG_VEX  = TUG_ISP * G0
TUG_THR  = 0.2          # N

# ── Mission ─────────────────────────────────────────────────────────────────
T_OPS       = 5.0       # proximity ops per debris [days]
DA_RATIO    = 0.02     # drift phasing: delta_a / a
PHASE_RAD   = np.pi/2   # 90 deg average phase gap
MAX_DAYS    = 365.0
SOFT_MASS   = 2000.0    # kg flag threshold

# ── Recycling Hub ───────────────────────────────────────────────────────────
RH_SMA  = 42164.2e3     # m  (nominal GEO)
RH_INC  = 7.0           # deg
RH_RAAN = 10.0          # deg

# ── Debris catalogue ────────────────────────────────────────────────────────
_RAW = parameter_table
# [
#     #  id    name                         mass_kg  sma_km      inc       raan
#     (489,  'NSS 5 (Intelsat 803)',        2060.46, 42494.423, 10.4813,  48.1496),
#     (495,  'Sirius 2 (GE-1E)',           1762.14, 42410.084, 11.9727,  34.1101),
#     (518,  'EUTE W2',                    1793.86, 42451.133, 11.6171,  38.4600),
#     (530,  'Skynet 4E',                  1490.00, 42516.912, 12.0550,   9.2220),
#     (556,  'EUTE 16C (SESAT 1)',         2600.00, 42537.778, 10.8393,  46.2495),
#     (593,  'Skynet 4F',                  1489.00, 42482.488, 12.1030,  18.6123),
#     (628,  'EUTE 12 West A (AB 1)',      2700.00, 42728.199,  7.2278,  68.7938),
#     (653,  'EUTE 33A (EB 3)',            1552.00, 42558.620,  9.6938,  55.8840),
#     (660,  'SESAT 2 (Express AM-22)',    2542.00, 42420.904,  8.2239,  64.3413),
#     (705,  'Meteosat 9 (MSG 2)',         2054.00, 42163.863,  9.3136,  54.3880),
#     (741,  'EUTE 8A (Sinosat 3)',        2320.00, 42716.964,  9.5162,  48.3713),
#     (804,  'COMSATBW-1',               2440.00, 42164.245,  0.0696,  89.6878),
#     (819,  'COMSATBW-2',               2440.00, 42164.548,  0.0344, 111.5455),
#     (834,  'HYLAS 1',                  2570.00, 42164.346,  5.5899,  73.2864),
#     (884,  'Meteosat 10 (MSG 3)',       2035.00, 42164.226,  4.6216,  61.1229),
#     (967,  'Meteosat 11 (MSG 4)',       2043.00, 42164.271,  3.1471,  71.2262),
# ]
N_DEB = len(_RAW)
IDS   = np.array([r[0] for r in _RAW])
NAMES = [r[1] for r in _RAW]
MASS  = np.array([r[2] for r in _RAW])          # kg
SMA   = np.array([r[3] for r in _RAW]) * 1e3    # m
INC   = np.array([r[4] for r in _RAW])          # deg
RAAN  = np.array([r[5] for r in _RAW])          # deg

# Augment with RH as index 16
SMA_ALL  = np.append(SMA,  RH_SMA)
INC_ALL  = np.append(INC,  RH_INC)
RAAN_ALL = np.append(RAAN, RH_RAAN)
MASS_ALL = np.append(MASS, 0.0)
RH_IDX   = N_DEB   # 16

# =============================================================================
# PRE-COMPUTE PAIRWISE TRANSFER TABLE  (17x17 objects incl. RH)
# =============================================================================

def build_transfer_table():
    """
    Pre-compute all pairwise transfer quantities between all 17 nodes
    (16 debris + RH).  Returns dicts of arrays indexed [from, to].
    """
    N = N_DEB + 1
    # Allocate
    dv1   = np.zeros((N, N))   # Hohmann departure burn [m/s]
    dv2   = np.zeros((N, N))   # Combined circ+plane at apoapsis [m/s]
    dv_ph = np.zeros((N, N))   # Drift phasing [m/s]
    t_tr  = np.zeros((N, N))   # Hohmann transfer time [days]
    t_ph  = np.zeros((N, N))   # Phasing time [days]

    for i in range(N):
        sa = SMA_ALL[i]
        for j in range(N):
            if i == j:
                continue
            sb = SMA_ALL[j]

            # ── Hohmann ──────────────────────────────────────────────────
            at   = (sa + sb) / 2.0
            va   = np.sqrt(MU / sa)
            vb   = np.sqrt(MU / sb)
            v_p  = np.sqrt(MU * (2.0/sa - 1.0/at))
            v_ap = np.sqrt(MU * (2.0/sb - 1.0/at))
            dv1[i,j]  = abs(v_p - va)
            t_tr[i,j] = np.pi * np.sqrt(at**3 / MU) / DAY

            # ── Combined circularisation + plane change at apoapsis ──────
            i1 = INC_ALL[i]*D2R;  o1 = RAAN_ALL[i]*D2R
            i2 = INC_ALL[j]*D2R;  o2 = RAAN_ALL[j]*D2R
            cos_dth = (np.cos(i1)*np.cos(i2) +
                       np.sin(i1)*np.sin(i2)*np.cos(o2-o1))
            dth = np.arccos(np.clip(cos_dth, -1.0, 1.0))
            dv2[i,j] = np.sqrt(v_ap**2 + vb**2 - 2*v_ap*vb*np.cos(dth))

            # ── Drift phasing ────────────────────────────────────────────
            vc_b  = vb
            n_b   = 2*np.pi / (2*np.pi*np.sqrt(sb**3/MU))
            rate  = 1.5 * n_b * DA_RATIO
            t_ph[i,j]  = (PHASE_RAD / rate) / DAY
            dv_ph[i,j] = DA_RATIO * vc_b    # 2 * 0.5 * DA_RATIO * vc_b

    return dv1, dv2, dv_ph, t_tr, t_ph

print("  Pre-computing transfer table...", end=' ', flush=True)
t0 = time.time()
DV1, DV2, DV_PH, T_TR, T_PH = build_transfer_table()
DV_LEG = DV1 + DV2 + DV_PH    # total ΔV per leg [m/s]
T_LEG  = T_TR + T_PH           # total time per leg [days]
print(f"done ({time.time()-t0:.2f}s)")

# Pre-compute tug spiral data for each debris object
# Iterative Edelbaum sizing: m_wet = TUG_DRY + m_debris + m_prop
print("  Pre-computing tug spirals...", end=' ', flush=True)
TUG_DV   = np.zeros(N_DEB)
TUG_TIME = np.zeros(N_DEB)
TUG_MPROP = np.zeros(N_DEB)
TUG_MWET  = np.zeros(N_DEB)

for k in range(N_DEB):
    m_pl = TUG_DRY + MASS[k]
    m_wet = m_pl * 1.35
    i1,o1 = INC[k], RAAN[k]
    i2,o2 = RH_INC, RH_RAAN
    s1,s2 = SMA[k], RH_SMA
    v1 = np.sqrt(MU/s1); v2 = np.sqrt(MU/s2)
    cos_dth = (np.cos(i1*D2R)*np.cos(i2*D2R) +
               np.sin(i1*D2R)*np.sin(i2*D2R)*np.cos((o2-o1)*D2R))
    dth = np.arccos(np.clip(cos_dth,-1,1))
    dv_e = np.sqrt(v1**2 + v2**2 - 2*v1*v2*np.cos(np.pi/2*dth))
    for _ in range(60):
        mf   = m_wet * np.exp(-dv_e/TUG_VEX)
        mp   = m_wet - mf
        mnew = m_pl + mp
        if abs(mnew - m_wet) < 0.05: break
        m_wet = 0.6*m_wet + 0.4*mnew
    t_s = (m_wet * TUG_VEX / TUG_THR) * (1 - np.exp(-dv_e/TUG_VEX))
    TUG_DV[k]    = dv_e
    TUG_TIME[k]  = t_s / DAY
    TUG_MPROP[k] = m_wet - m_pl
    TUG_MWET[k]  = m_wet

print(f"done ({time.time()-t0:.2f}s)")

# Phasing time for mothership to meet tug at RH (one leg)
def drift_phase_scalar(sma):
    n = np.sqrt(MU/sma**3)
    rate = 1.5 * n * DA_RATIO
    return DA_RATIO * np.sqrt(MU/sma), PHASE_RAD/rate/DAY

_, T_PH_RH = drift_phase_scalar(RH_SMA)

# =============================================================================
# VECTORISED SEQUENCE EVALUATOR
# =============================================================================

def evaluate_all_sequences(ms_prop):
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

    ms_wet0 = MS_DRY + ms_prop

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
    t_legs  = T_LEG[from_nodes,  to_nodes]   # (n_seq, 6)

    # Sequential rocket equation for mass
    # After each burn the mass drops. Apply log sum trick.
    # m_after_n = m_wet0 * exp(-sum(dv[0..n]) / MS_VEX)
    cum_dv     = np.cumsum(dv_legs, axis=1)            # (n_seq, 6)
    mass_after = ms_wet0 * np.exp(-cum_dv / MS_VEX)   # (n_seq, 6)

    # Check propellant not exhausted at any point
    prop_ok = np.all(mass_after >= MS_DRY, axis=1)    # (n_seq,)

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
    tot_dv  = cum_dv[:, 5]        # (n_seq,)
    tot_prop = ms_wet0 - mass_after[:, 5]

    print(f"done ({time.time()-t1:.2f}s)")

    # ── Phase 2: Tug spirals (parallel) ────────────────────────────────────
    print("  Evaluating tug spirals...", end=' ', flush=True)
    t2 = time.time()

    # Tug arrival = tug_start + TUG_TIME[debris_index]
    tug_spiral = TUG_TIME[sequences]             # (n_seq, 5)
    tug_arrive = tug_start + tug_spiral          # (n_seq, 5)

    # Handover day = max(ms_return, tug_arrive) + T_PH_RH (MS phases to tug)
    handover = np.maximum(ms_return[:, None], tug_arrive) + T_PH_RH  # (n_seq,5)

    # Mission completion = last handover
    mission_day = handover.max(axis=1)           # (n_seq,)

    print(f"done ({time.time()-t2:.2f}s)")

    # ── Feasibility & scoring ───────────────────────────────────────────────
    time_ok  = mission_day <= MAX_DAYS
    feasible = prop_ok & time_ok
    f_idx    = np.where(feasible)[0]

    print(f"  Feasible: {len(f_idx):,}  |  Infeasible: {n_seq-len(f_idx):,}")

    if len(f_idx) == 0:
        return None

    # Extract feasible subset
    f_seq    = sequences[f_idx]           # (nf, 5)
    f_dv     = tot_dv[f_idx]
    f_prop   = tot_prop[f_idx]
    f_mrem   = ms_wet0 - mass_after[f_idx, 5] # same as f_prop
    f_day    = mission_day[f_idx]
    f_mass   = MASS[f_seq].sum(axis=1)
    f_mret   = ms_return[f_idx]
    f_tstart = tug_start[f_idx]
    f_tarr   = tug_arrive[f_idx]
    f_ho     = handover[f_idx]
    f_heavy  = (MASS[f_seq] > SOFT_MASS).any(axis=1)
    f_mafter = mass_after[f_idx, 5]

    # Pareto score: 40% ΔV + 40% time + 20% mass (normalised)
    dv_n   = (f_dv  - f_dv.min())  / (f_dv.max()  - f_dv.min()  + 1e-9)
    time_n = (f_day - f_day.min()) / (f_day.max() - f_day.min() + 1e-9)
    mass_n = (f_mass.max() - f_mass) / (f_mass.max()-f_mass.min()+1e-9)
    score  = 0.40*dv_n + 0.40*time_n + 0.20*mass_n

    order  = np.argsort(score)

    results = {
        'sequences':   f_seq[order],
        'score':       score[order],
        'total_dv':    f_dv[order],
        'prop_used':   f_prop[order],
        'prop_margin': f_mafter[order] - MS_DRY,
        'mission_day': f_day[order],
        'ms_return':   f_mret[order],
        'mass_removed':f_mass[order],
        'tug_start':   f_tstart[order],
        'tug_arrive':  f_tarr[order],
        'handover':    f_ho[order],
        'has_heavy':   f_heavy[order],
        'n_feasible':  len(f_idx),
        'n_total':     n_seq,
        'ms_prop':     ms_prop,
        # Per-leg ΔV breakdown for top results
        'dv_legs':     dv_legs[f_idx][order],
        't_legs':      t_legs[f_idx][order],
        'cum_time':    cum_time[f_idx][order],
    }
    return results

# =============================================================================
# RESULTS PRINTER
# =============================================================================

def print_top(res, n=10):
    print("\n"+"="*70)
    print(f"  TOP {n} SEQUENCES  (40% ΔV + 40% time + 20% mass removed)")
    print("="*70)

    for rank in range(min(n, res['n_feasible'])):
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
        for i, idx in enumerate(seq):
            print(f"  {i+1:<3} {NAMES[idx]:<28} "
                  f"{TUG_DV[idx]:>7.1f}m "
                  f"{TUG_MPROP[idx]:>7.1f}kg "
                  f"{TUG_TIME[idx]:>8.1f}d "
                  f"{res['tug_arrive'][rank,i]:>8.1f}d "
                  f"{res['handover'][rank,i]:>8.1f}d")

        print(f"\n  ┌────────────────────────────────────────────────────────┐")
        print(f"  │  Mothership ΔV total     : {res['total_dv'][rank]:>8.1f} m/s             │")
        print(f"  │  Mothership prop used    : {res['prop_used'][rank]:>8.1f} kg              │")
        print(f"  │  Mothership prop margin  : {res['prop_margin'][rank]:>8.1f} kg              │")
        print(f"  │  MS returns to RH day    : {res['ms_return'][rank]:>8.1f}                │")
        print(f"  │  Total debris mass       : {res['mass_removed'][rank]:>8.0f} kg              │")
        print(f"  │  Mission completion      : {res['mission_day'][rank]:>8.1f} days  ✓       │")
        print(f"  └────────────────────────────────────────────────────────┘")

# =============================================================================
# SENSITIVITY
# =============================================================================

def sensitivity(res):
    global T_OPS, T_PH_RH, T_LEG, T_PH, DV_LEG
    best_seq = res['sequences'][0]
    ms_prop  = res['ms_prop']

    print("\n"+"="*70)
    print("  SENSITIVITY ANALYSIS  —  best sequence")
    print("="*70)
    print("  "+"  ".join(f"{i+1}.{NAMES[best_seq[i]].split('(')[0][:12]}"
                          for i in range(5)))

    # A) Propellant sweep
    print(f"\n  A) Propellant budget (T_ops={T_OPS}d fixed):")
    print(f"  {'Prop':>8} {'Feas':>5} {'Compl d':>9} {'Used kg':>9} {'Margin':>8}")
    print(f"  {'':─>8} {'':─>5} {'':─>9} {'':─>9} {'':─>8}")
    for p in np.linspace(500, 4000, 15):
        mw = MS_DRY + p
        cum_dv_seq = np.cumsum([DV_LEG[([RH_IDX]+list(best_seq))[i],
                                        (list(best_seq)+[RH_IDX])[i]]
                                 for i in range(6)])
        ma = mw * np.exp(-cum_dv_seq / MS_VEX)
        ok = (ma >= MS_DRY).all()
        if ok:
            prop_used = mw - ma[-1]
            margin    = ma[-1] - MS_DRY
            # Recompute mission day for this prop
            r2 = _eval_single(best_seq, p)
            print(f"  {p:>8.0f} {'✓':>5} {r2['day']:>9.1f} {prop_used:>9.1f} {margin:>8.1f}")
        else:
            print(f"  {p:>8.0f} {'✗':>5} {'—':>9} {'—':>9} {'—':>8}  propellant exhausted")

    # B) Ops time sweep
    orig_tops = T_OPS
    print(f"\n  B) Ops time per debris (prop={ms_prop}kg fixed):")
    print(f"  {'T_ops d':>8} {'Feas':>5} {'Compl d':>9}")
    print(f"  {'':─>8} {'':─>5} {'':─>9}")
    for t in np.linspace(1, 14, 14):
        T_OPS = t
        r2 = _eval_single(best_seq, ms_prop)
        if r2['feasible']:
            print(f"  {t:>8.1f} {'✓':>5} {r2['day']:>9.1f}")
        else:
            print(f"  {t:>8.1f} {'✗':>5} {'—':>9}  {r2['reason']}")
    T_OPS = orig_tops

def _eval_single(seq, ms_prop):
    """Scalar evaluation of one sequence for sensitivity sweeps."""
    nodes_from = [RH_IDX] + list(seq)
    nodes_to   = list(seq) + [RH_IDX]
    mw = MS_DRY + ms_prop
    day = 0.0
    for i in range(6):
        fi, ti = nodes_from[i], nodes_to[i]
        dv = DV_LEG[fi, ti]
        mw = mw * np.exp(-dv / MS_VEX)
        if mw < MS_DRY:
            return {'feasible': False, 'reason': 'propellant exhausted'}
        day += T_LEG[fi, ti]
        if i < 5:
            day += T_OPS
    ms_ret = day
    tug_starts = []
    d = 0.0
    for i in range(6):
        fi, ti = nodes_from[i], nodes_to[i]
        d += T_LEG[fi, ti]
        if i < 5:
            d += T_OPS
            tug_starts.append(d)
    handovers = [max(ms_ret, tug_starts[i] + TUG_TIME[seq[i]]) + T_PH_RH
                 for i in range(5)]
    mission_day = max(handovers)
    feasible = mission_day <= MAX_DAYS
    return {'feasible': feasible, 'day': mission_day,
            'reason': f'{mission_day:.1f}d > {MAX_DAYS}d'}

# =============================================================================
# PLOTTING
# =============================================================================

def make_plots(res, save_path='/mnt/user-data/outputs/reaver_optimizer_results.png'):
    BG='#0d1117'; CB='#161b22'; TC='#c9d1d9'; MU_='#8b949e'; GR='#21262d'
    A1='#58a6ff'; A2='#3fb950'; A3='#f78166'; A4='#d2a8ff'; A5='#ffa657'
    LC=[A1,A2,A4,A3,A5,'#79c0ff']; TC_=[A1,A2,A4,A3,A5]

    fig = plt.figure(figsize=(20,15)); fig.patch.set_facecolor(BG)
    gs  = GridSpec(3,3,figure=fig,hspace=0.48,wspace=0.36)

    def sax(ax,title=''):
        ax.set_facecolor(CB)
        [s.set_edgecolor(GR) for s in ax.spines.values()]
        ax.tick_params(colors=MU_,labelsize=8)
        ax.xaxis.label.set_color(MU_); ax.yaxis.label.set_color(MU_)
        ax.grid(True,color=GR,lw=0.5,alpha=0.7)
        if title: ax.set_title(title,color=TC,fontsize=9,fontweight='bold',pad=8)

    fig.text(0.5,0.976,'REAVER Mission Optimizer — Results Dashboard',
             ha='center',color=TC,fontsize=15,fontweight='bold')
    fig.text(0.5,0.953,
             f'MS: dry {MS_DRY}kg, Isp={MS_ISP}s  |  '
             f'Tug: dry {TUG_DRY}kg, Isp={TUG_ISP}s, T={TUG_THR}N  |  '
             f'T_ops={T_OPS}d/debris  |  90° avg phasing  |  ≤{MAX_DAYS:.0f}d',
             ha='center',color=MU_,fontsize=8)

    dv   = res['total_dv']
    tday = res['mission_day']
    mass = res['mass_removed']

    # ── 1. Pareto scatter ─────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0,:2])
    sax(ax1,'Pareto space — ΔV vs mission duration  (colour = mass removed)')
    sc = ax1.scatter(tday,dv,c=mass,cmap='plasma',alpha=0.3,s=4,lw=0)
    ax1.scatter(tday[:10],dv[:10],color=A2,s=70,zorder=5,
                edgecolors='white',lw=0.6,label='Top 10')
    for i,(t,d) in enumerate(zip(tday[:5],dv[:5])):
        ax1.annotate(f'#{i+1}',(t,d),xytext=(6,4),
                     textcoords='offset points',fontsize=8,color=A2,fontweight='bold')
    ax1.axvline(MAX_DAYS,color=A3,lw=1.3,ls='--',alpha=0.9,label=f'{MAX_DAYS:.0f}d limit')
    cb=plt.colorbar(sc,ax=ax1,pad=0.02)
    cb.set_label('Mass removed [kg]',color=MU_,fontsize=8)
    cb.ax.yaxis.set_tick_params(color=MU_,labelsize=7)
    plt.setp(cb.ax.yaxis.get_ticklabels(),color=MU_)
    ax1.set_xlabel('Mission completion [days]'); ax1.set_ylabel('Mothership ΔV [m/s]')
    ax1.legend(fontsize=7.5,facecolor=CB,labelcolor=TC,edgecolor=GR)

    # ── 2. Duration histogram ─────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[0,2])
    sax(ax2,'Mission duration distribution')
    ax2.hist(tday,bins=45,color=A1,alpha=0.75,ec='none')
    ax2.axvline(MAX_DAYS,color=A3,lw=1.5,ls='--',label=f'{MAX_DAYS:.0f}d')
    ax2.axvline(tday[:10].mean(),color=A2,lw=1.2,ls=':',
                label=f'Top10 avg {tday[:10].mean():.0f}d')
    ax2.set_xlabel('Completion day'); ax2.set_ylabel('Count')
    ax2.legend(fontsize=7,facecolor=CB,labelcolor=TC,edgecolor=GR)

    # ── 3. ΔV stacked bar top-10 ──────────────────────────────────────────
    ax3 = fig.add_subplot(gs[1,:2])
    sax(ax3,'Top 10 — mothership ΔV breakdown per leg')
    for ri in range(min(10,res['n_feasible'])):
        bot=0.0
        seq_r = res['sequences'][ri]
        nodes_f=[RH_IDX]+list(seq_r); nodes_t=list(seq_r)+[RH_IDX]
        for li in range(6):
            dv_l=DV_LEG[nodes_f[li],nodes_t[li]]
            col=LC[li] if li<5 else '#484f58'
            ax3.bar(ri,dv_l,bottom=bot,width=0.6,color=col,alpha=0.85,ec='none')
            if dv_l>25:
                ax3.text(ri,bot+dv_l/2,f'{dv_l:.0f}',
                         ha='center',va='center',fontsize=6,color='white',fontweight='bold')
            bot+=dv_l
    ax3.set_xticks(range(min(10,res['n_feasible'])))
    ax3.set_xticklabels([f'#{i+1}' for i in range(min(10,res['n_feasible']))])
    ax3.set_xlabel('Rank'); ax3.set_ylabel('ΔV [m/s]')
    pats=([mpatches.Patch(color=LC[i],label=f'Leg {i+1}') for i in range(5)]
          +[mpatches.Patch(color='#484f58',label='Return')])
    ax3.legend(handles=pats,fontsize=7,facecolor=CB,labelcolor=TC,edgecolor=GR,ncol=3)

    # ── 4. Best solution tug Gantt ────────────────────────────────────────
    ax4 = fig.add_subplot(gs[1,2])
    sax(ax4,'#1 solution — tug spiral Gantt')
    seq0=res['sequences'][0]
    for i,idx in enumerate(seq0):
        s=res['tug_start'][0,i]; e=res['tug_arrive'][0,i]
        nm=NAMES[idx].split('(')[0].strip()[:20]
        ax4.barh(i,e-s,left=s,color=TC_[i],alpha=0.8,ec='none',height=0.55)
        ax4.text(e+1.5,i,f'd{e:.0f}',va='center',fontsize=7,color=TC_[i])
    ax4.axvline(res['ms_return'][0],color='white',lw=1.2,ls=':',alpha=0.7,
                label=f"MS@RH d{res['ms_return'][0]:.0f}")
    ax4.axvline(MAX_DAYS,color=A3,lw=1.3,ls='--',alpha=0.85,label='365d')
    ax4.set_yticks(range(5))
    ax4.set_yticklabels([NAMES[i].split('(')[0].strip()[:18] for i in seq0],fontsize=7)
    ax4.set_xlabel('Mission day')
    ax4.legend(fontsize=7,facecolor=CB,labelcolor=TC,edgecolor=GR)

    # ── 5. Mission timeline ribbon ────────────────────────────────────────
    ax5 = fig.add_subplot(gs[2,:2])
    sax(ax5,'#1 solution — full mission timeline')
    events=[(0,'Depart RH',MU_)]
    cum=0.0
    for i,idx in enumerate(seq0):
        fi=[RH_IDX]+list(seq0)
        ti=list(seq0)+[RH_IDX]
        cum+=T_TR[fi[i],ti[i]]+T_PH[fi[i],ti[i]]
        events.append((cum,f'Arr.{i+1}:{NAMES[idx].split("(")[0][:9]}',LC[i]))
        cum+=T_OPS
        events.append((cum,f'Handoff Tug{i+1}',TC_[i]))
    cum+=T_TR[seq0[-1],RH_IDX]+T_PH[seq0[-1],RH_IDX]
    events.append((cum,'MS@RH',A2))
    for i in range(5):
        events.append((res['handover'][0,i],
                       f'HO:{NAMES[seq0[i]].split("(")[0][:9]}',A4))
    events.append((res['mission_day'][0],'✓ DONE',A3))
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

    # ── 6. Debris frequency top-50 ───────────────────────────────────────
    ax6 = fig.add_subplot(gs[2,2])
    sax(ax6,'Debris selection frequency (top 50 solutions)')
    freq=Counter()
    for r in range(min(50,res['n_feasible'])):
        for idx in res['sequences'][r]:
            freq[NAMES[idx].split('(')[0].strip()[:18]]+=1
    ns=sorted(freq,key=freq.get,reverse=True)[:10]
    cs=[freq[n] for n in ns]
    bars=ax6.barh(range(len(ns)),cs,color=A1,alpha=0.8,ec='none')
    ax6.set_yticks(range(len(ns))); ax6.set_yticklabels(ns,fontsize=7)
    ax6.set_xlabel('Appearances in top 50')
    for bar,cnt in zip(bars,cs):
        ax6.text(bar.get_width()+0.3,bar.get_y()+bar.get_height()/2,
                 str(cnt),va='center',fontsize=7,color=A1)

    plt.savefig(save_path,dpi=150,bbox_inches='tight',
                facecolor=BG,edgecolor='none')
    print(f"\n  Plot saved → {save_path}")

# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    MS_PROP = 2000.0

    print("="*70)
    print("  REAVER TRAJECTORY OPTIMIZER  —  Group 7")
    print("="*70)
    print(f"  MS dry {MS_DRY}kg | MS Isp {MS_ISP}s | Prop budget {MS_PROP}kg")
    print(f"  Tug dry {TUG_DRY}kg | Tug Isp {TUG_ISP}s | Thrust {TUG_THR}N")
    print(f"  T_ops {T_OPS}d | Constraint ≤{MAX_DAYS}d | 90° avg phasing")

    t0=time.time()
    res = evaluate_all_sequences(MS_PROP)
    print(f"  Total runtime: {time.time()-t0:.2f}s")

    if res is None:
        print("\n  ⚠  No feasible solutions. Trying 3500kg propellant...")
        MS_PROP=3500.0
        res = evaluate_all_sequences(MS_PROP)

    if res:
        print_top(res, n=10)
        sensitivity(res)
        make_plots(res)

        rpt='/mnt/user-data/outputs/reaver_top_sequences.txt'
        orig=sys.stdout
        with open(rpt,'w') as f:
            sys.stdout=f
            print_top(res,n=10)
            sensitivity(res)
        sys.stdout=orig
        print(f"  Report saved → {rpt}")

    print("\n  Done.")
