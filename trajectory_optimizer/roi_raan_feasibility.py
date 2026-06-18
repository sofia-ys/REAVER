"""
REAVER — ROI Propellant-Feasibility of the Extended Debris Pool
===============================================================
Selling-pitch analysis for mission extension (ROI Section 19.4).

Question
--------
With the Recycling Hub held at its design RAAN of 64 deg (station-kept, no
re-optimisation) and the mothership built to its DESIGN propellant budget,
how many 5-target campaigns over the *extended* debris catalogue (33 objects)
are actually feasible — i.e. fit inside the as-built tank AND the 365-day limit?

Design propellant budget (from the nominal mission 804->884->697->488->443):
    MS orbital propellant      = 1019.3 kg   (design point, RCS reserve intact)
    MS orbital + 15% RCS margin = 1172.2 kg  (stretch: burn into the RCS reserve)

Both caps act on the per-campaign MS *orbital* propellant (transfers only),
which is the quantity the 1019.3 kg figure represents (RPO dV is excluded from
that number in nominal_mission.py / reaver_optimizer.py).

Capacity metric
---------------
A 5-target campaign = one operational year (+EUR290M net cash flow).  Because
every object is removed exactly once, the headline number is the count of
*disjoint* feasible campaigns (operational years), found by greedy packing,
NOT the raw count of feasible combinations.

This script is self-contained: it does not import the 21-object catalogue and
does not modify config.py, reaver_core.py, or rh_raan_sweep.py.
"""

import os
import time
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Scalars only — NOT the 21-object catalogue — from the single source of truth.
from config import (MU, DAY, D2R, MS_DRY, MS_VEX, TUG_DRY, TUG_VEX, TUG_THR,
                    TUG_PROP_MARGIN, N_PHASE_REV, T_OPS, MAX_DAYS,
                    RH_SMA, RH_INC, MS_THR, MS_BURN_S, MS_WET_ESTIMATE)

# =============================================================================
# CONFIGURATION FOR THIS ANALYSIS
# =============================================================================

RH_RAAN_FIXED  = 64.0      # deg   hub held here by station-keeping (NOT optimised)
MS_PROP_DESIGN = 1019.3    # kg    design MS orbital propellant (RCS reserve intact)
MS_PROP_MARGIN = 1172.2    # kg    orbital + 15% RCS margin (stretch cap)

# =============================================================================
# EXTENDED DEBRIS CATALOGUE  (33 objects = original 21 + 12 lower-RAAN additions)
#   columns used: norad_id, name, mass_kg, sma_km, inc_deg, raan_deg
# =============================================================================

_RAW = [
    (443,  'ECHOSTAR 1',              1902.87, 42544.807,  8.5674, 63.1596),
    (464,  'AMC-1 (GE-1)',            1585.73, 42380.233,  9.0706, 59.3739),
    (488,  'AMC-3 (GE-3)',            1585.73, 42165.071,  8.0920, 64.3488),
    (490,  'ECHOSTAR 3',              2125.87, 42560.558,  9.4764, 55.9812),
    (498,  'ASTRA 1G',                2279.48, 42514.794,  9.4583, 55.0755),
    (505,  'NSS 806 (INTELSAT 806)',  2135.78, 42502.032,  7.0712, 69.9995),
    (513,  'INTELSAT 805',            1932.61, 42388.082,  7.8332, 66.0575),
    (517,  'INTELSAT 7 (PAS 7)',      2099.11, 42752.715, 10.3941, 52.3204),
    (520,  'EUTE 25A (HB 5)',         1744.30, 42681.392, 10.4738, 50.9476),
    (525,  'EUTE 115 WEST A(SATMEX 5',2229.93, 42604.474, 10.2863, 52.1720),
    (554,  'ASIASTAR',                2777.00, 42164.171,  9.0818, 58.3898),
    (579,  'AMC-7 (GE-7)',            1935.00, 42486.971,  6.8592, 70.5256),
    (590,  'AMC-8 (GE-8)',            2015.00, 42482.105,  6.0647, 72.8264),
    (628,  'EUTE 12 WEST A (AB 1)',   2700.00, 42728.199,  7.2278, 68.7938),
    (640,  'GALAXY 12',               1760.00, 42528.577,  6.9338, 69.5405),
    (653,  'EUTE 33A (EB 3)',         1552.00, 42558.620,  9.6938, 55.8840),
    (663,  'AMC-10 (GE-10)',          2315.00, 42490.063,  6.6355, 69.2760),
    (672,  'AMC-11 (GE-11)',          2315.00, 42499.585,  0.2335, 80.2545),
    (697,  'GALAXY 14',               2087.00, 42549.012,  4.9754, 75.3425),
    (699,  'GALAXY 15',               2033.00, 42542.743,  3.4810, 77.4994),
    (705,  'METEOSAT 9 (MSG 2)',      2054.00, 42163.863,  9.3136, 54.3880),
    (731,  'AMC-18',                  2081.00, 42565.455,  2.5582, 84.3917),
    (748,  'INTELSAT 11 (PAS 11)',    2450.00, 42510.914,  3.5105, 78.9272),
    (755,  'HORIZONS 2',              2304.00, 42164.580,  2.5479, 80.8867),
    (774,  'AMC-21',                  2473.00, 42522.009,  2.8956, 80.3233),
    (786,  'NSS 9',                   2290.00, 42164.923,  1.2073, 82.8539),
    (804,  'COMSATBW-1',              2440.00, 42164.245,  0.0696, 89.6878),
    (834,  'HYLAS 1',                 2570.00, 42164.346,  5.5899, 73.2864),
    (884,  'METEOSAT 10 (MSG 3)',     2035.00, 42164.226,  4.6216, 61.1229),
    (942,  'DELTA 4 R/B',             2850.00, 42549.844,  9.1489, 57.5692),
    (967,  'METEOSAT 11 (MSG 4)',     2043.00, 42164.271,  3.1471, 71.2262),
    (1009, 'DELTA 4 R/B',             2850.00, 42527.541,  8.6610, 58.1750),
    (1169, 'SES 20',                  1499.00, 42165.239,  0.0464, 98.6035),
]

# NORAD ids of the 12 objects added beyond the original 21-object design pool.
_NEW12_IDS = {464, 490, 498, 517, 520, 525, 554, 653, 705, 942, 1009, 1169}

N_DEB = len(_RAW)
IDS   = np.array([r[0] for r in _RAW])
NAMES = [r[1] for r in _RAW]
MASS  = np.array([r[2] for r in _RAW])
SMA   = np.array([r[3] for r in _RAW]) * 1e3   # m
INC   = np.array([r[4] for r in _RAW])         # deg
RAAN  = np.array([r[5] for r in _RAW])         # deg

# Augment with the RH as node index N_DEB.
RH_IDX   = N_DEB
SMA_ALL  = np.append(SMA,  RH_SMA)
INC_ALL  = np.append(INC,  RH_INC)
RAAN_ALL = np.append(RAAN, RH_RAAN_FIXED)
IS_NEW12 = np.array([i in _NEW12_IDS for i in IDS])   # (N_DEB,) bool

# =============================================================================
# PHYSICS  (replicated from reaver_core for the 33-object set, RH at 64 deg)
# =============================================================================

def build_tables(rh_raan_deg):
    """Pairwise impulsive transfer dV and component times for all 34 nodes."""
    raan_all = RAAN_ALL.copy()
    raan_all[RH_IDX] = rh_raan_deg
    N = N_DEB + 1
    dv1 = np.zeros((N, N)); dv2 = np.zeros((N, N)); dv_ph = np.zeros((N, N))
    t_tr = np.zeros((N, N)); t_ph = np.zeros((N, N))
    for i in range(N):
        sa = SMA_ALL[i]; ia = INC_ALL[i]*D2R; oa = raan_all[i]*D2R
        va = np.sqrt(MU/sa)
        for j in range(N):
            if i == j:
                continue
            sb = SMA_ALL[j]; ib = INC_ALL[j]*D2R; ob = raan_all[j]*D2R
            vb = np.sqrt(MU/sb)
            at   = (sa + sb)/2.0
            v_sa = np.sqrt(MU*(2.0/sa - 1.0/at))
            v_sb = np.sqrt(MU*(2.0/sb - 1.0/at))
            t_tr[i, j] = np.pi*np.sqrt(at**3/MU)/DAY
            cos_dth = (np.cos(ia)*np.cos(ib) +
                       np.sin(ia)*np.sin(ib)*np.cos(ob - oa))
            dth = np.arccos(np.clip(cos_dth, -1.0, 1.0))
            if sb > sa:
                dv1[i, j] = abs(v_sa - va)
                dv2[i, j] = np.sqrt(v_sb**2 + vb**2 - 2*v_sb*vb*np.cos(dth))
            else:
                dv1[i, j] = np.sqrt(v_sa**2 + va**2 - 2*v_sa*va*np.cos(dth))
                dv2[i, j] = abs(v_sb - vb)
            T_tgt      = 2*np.pi*np.sqrt(sb**3/MU)
            T_ph_orb   = T_tgt*(1.0 - 1.0/(4.0*N_PHASE_REV))
            a_ph       = (MU*(T_ph_orb/(2*np.pi))**2)**(1.0/3.0)
            dv_ph[i, j]= 2*abs(np.sqrt(MU/sb) - np.sqrt(MU*(2/sb - 1/a_ph)))
            t_ph[i, j] = N_PHASE_REV*T_ph_orb/DAY
    return dv1, dv2, dv_ph, t_tr, t_ph


def finite_burn_time(dv_req, mass_kg, sma_m):
    """Wall-clock days to deliver dv_req: 64 N thrust, 38-min fires, 2 fires/orbit."""
    if dv_req <= 0.0:
        return 0.0
    dv_per_fire = (MS_THR/mass_kg)*MS_BURN_S
    n_fires  = int(np.ceil(dv_req/dv_per_fire))
    n_orbits = int(np.ceil(n_fires/2))
    T_orb    = 2*np.pi*np.sqrt(sma_m**3/MU)
    return n_orbits*T_orb/DAY


def build_finite_leg_table(dv1, dv2, dv_ph, t_tr, t_ph, mass_kg):
    """Finite-burn-corrected leg times (matches reaver_optimizer design model)."""
    N = N_DEB + 1
    t_fin = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            m1 = mass_kg
            m2 = m1*np.exp(-dv1[i, j]/MS_VEX)
            m3 = m2*np.exp(-dv2[i, j]/MS_VEX)
            t_fin[i, j] = (finite_burn_time(dv1[i, j], m1, SMA_ALL[i])
                         + t_tr[i, j]
                         + finite_burn_time(dv2[i, j], m2, SMA_ALL[j])
                         + finite_burn_time(dv_ph[i, j], m3, SMA_ALL[j])
                         + t_ph[i, j])
    return t_fin


def compute_tug_spirals(rh_raan_deg):
    """Edelbaum tug spiral sizing (debris -> RH) for all 33 objects."""
    tug_time  = np.zeros(N_DEB)
    tug_mprop = np.zeros(N_DEB)
    v2 = np.sqrt(MU/RH_SMA)
    i2, o2 = RH_INC, rh_raan_deg
    for k in range(N_DEB):
        m_pl  = TUG_DRY + MASS[k]
        m_wet = m_pl*1.35
        i1, o1 = INC[k], RAAN[k]
        v1 = np.sqrt(MU/SMA[k])
        cos_d = (np.cos(i1*D2R)*np.cos(i2*D2R) +
                 np.sin(i1*D2R)*np.sin(i2*D2R)*np.cos((o2 - o1)*D2R))
        dth  = np.arccos(np.clip(cos_d, -1.0, 1.0))
        dv_e = np.sqrt(v1**2 + v2**2 - 2*v1*v2*np.cos(np.pi/2*dth))
        for _ in range(60):
            mf   = m_wet*np.exp(-dv_e/TUG_VEX)
            mnew = m_pl + (m_wet - mf)
            if abs(mnew - m_wet) < 0.05:
                break
            m_wet = 0.6*m_wet + 0.4*mnew
        t_s = (m_wet*TUG_VEX/TUG_THR)*(1 - np.exp(-dv_e/TUG_VEX))
        tug_time[k]  = t_s/DAY
        tug_mprop[k] = m_wet - m_pl
    return tug_time, tug_mprop


def tug_loads(tug_props):
    """REAVER loading rule: 2 tugs at (max req)*margin, 3 at (2nd req)*margin."""
    p = np.atleast_2d(np.asarray(tug_props, dtype=np.float64))
    order  = np.argsort(-p, axis=1)
    desc   = np.take_along_axis(p, order, axis=1)
    factor = 1.0 + TUG_PROP_MARGIN
    load_by_rank = np.empty_like(p)
    load_by_rank[:, :2] = desc[:, 0:1]*factor
    load_by_rank[:, 2:] = desc[:, 1:2]*factor
    loads = np.empty_like(p)
    np.put_along_axis(loads, order, load_by_rank, axis=1)
    return loads

# =============================================================================
# PRE-COMPUTE  (RH fixed at 64 deg)
# =============================================================================

print("  Building transfer & finite-burn tables (RH @ 64 deg)...", end=' ', flush=True)
t0 = time.time()
DV1, DV2, DV_PH, T_TR, T_PH = build_tables(RH_RAAN_FIXED)
DV_LEG = DV1 + DV2 + DV_PH
T_LEG_FINITE = build_finite_leg_table(DV1, DV2, DV_PH, T_TR, T_PH, MS_WET_ESTIMATE)
TUG_TIME, TUG_MPROP = compute_tug_spirals(RH_RAAN_FIXED)
print(f"done ({time.time()-t0:.2f}s)")

# =============================================================================
# EVALUATE EVERY 5-TARGET COMBINATION (best time-feasible ordering per combo)
# =============================================================================

_combos = np.array(list(itertools.combinations(range(N_DEB), 5)), dtype=np.int32)
_perms  = np.array(list(itertools.permutations(range(5))),        dtype=np.int32)
N_COMBOS = len(_combos)
N_PERMS  = len(_perms)
print(f"  Combinations C({N_DEB},5) = {N_COMBOS:,}  x {N_PERMS} orderings "
      f"= {N_COMBOS*N_PERMS:,} sequences")

T_OPS_ARR = np.array([T_OPS]*5 + [0.0])

# Per-combo result: the minimum orbital propellant achievable by any ordering
# that also satisfies the 365-day limit (inf if no ordering is time-feasible).
min_prop = np.full(N_COMBOS, np.inf)
best_seq = np.full((N_COMBOS, 5), -1, dtype=np.int32)

CHUNK = 3000   # combos per block -> ~360k sequences/block (memory-safe)
print("  Evaluating campaigns...", end=' ', flush=True)
t0 = time.time()
for s in range(0, N_COMBOS, CHUNK):
    cc = _combos[s:s+CHUNK]                       # (c, 5)
    c  = len(cc)
    seqs = cc[:, _perms].reshape(c*N_PERMS, 5)    # (c*120, 5) all orderings
    ns = len(seqs)

    frm = np.empty((ns, 6), dtype=np.int32)
    to  = np.empty((ns, 6), dtype=np.int32)
    frm[:, 0] = RH_IDX; frm[:, 1:5] = seqs[:, :4]; frm[:, 5] = seqs[:, 4]
    to[:, :5] = seqs;   to[:, 5]    = RH_IDX

    dv_legs = DV_LEG[frm, to]                     # (ns, 6)
    t_legs  = T_LEG_FINITE[frm, to]
    tug_mwet = (TUG_DRY + tug_loads(TUG_MPROP[seqs])).astype(np.float64)  # (ns, 5)

    # Backward pass -> MS orbital propellant (transfers only; final mass = MS_DRY)
    m_req = np.full(ns, MS_DRY)
    m_req = m_req*np.exp(dv_legs[:, 5]/MS_VEX)
    for leg in range(4, -1, -1):
        m_req += tug_mwet[:, leg]
        m_req  = m_req*np.exp(dv_legs[:, leg]/MS_VEX)
    ms_prop = m_req - MS_DRY - tug_mwet.sum(axis=1)

    # Forward pass -> mission duration (incl. tug spiral handovers)
    cum_time    = np.cumsum(t_legs + T_OPS_ARR[None, :], axis=1)
    tug_arrive  = cum_time[:, :5] + TUG_TIME[seqs]
    handover    = np.maximum(cum_time[:, 5:6], tug_arrive) + T_OPS
    mission_day = handover.max(axis=1)

    ms_prop_tf = np.where(mission_day <= MAX_DAYS, ms_prop, np.inf)
    ms_prop_tf = ms_prop_tf.reshape(c, N_PERMS)
    bperm  = ms_prop_tf.argmin(axis=1)
    bval   = ms_prop_tf[np.arange(c), bperm]
    min_prop[s:s+c] = bval
    best_seq[s:s+c] = seqs.reshape(c, N_PERMS, 5)[np.arange(c), bperm]
print(f"done ({time.time()-t0:.1f}s)")

# =============================================================================
# DISJOINT-COHORT PACKING  ->  operational years
# =============================================================================

def pack_cohorts(cap):
    """Greedily select disjoint feasible campaigns (lowest propellant first)."""
    mask  = min_prop <= cap
    idx   = np.where(mask)[0]
    order = idx[np.argsort(min_prop[idx])]
    used  = np.zeros(N_DEB, dtype=bool)
    chosen = []
    for ci in order:
        s = best_seq[ci]
        if not used[s].any():
            used[s] = True
            chosen.append(ci)
    return chosen

time_feasible = np.isfinite(min_prop)
caps = [('Design (1019.3 kg)', MS_PROP_DESIGN),
        ('Margin (1172.2 kg)', MS_PROP_MARGIN)]

print("\n" + "="*72)
print(f"  EXTENDED-POOL FEASIBILITY  |  RH RAAN = {RH_RAAN_FIXED:.0f} deg  |  "
      f"{N_DEB} objects ({int(IS_NEW12.sum())} new)")
print("="*72)
print(f"  Time-feasible combinations (<= {MAX_DAYS:.0f} d) : "
      f"{int(time_feasible.sum()):,} / {N_COMBOS:,}")
print(f"  Theoretical ceiling                    : {N_DEB//5} operational years "
      f"({N_DEB} objects / 5)")
print("-"*72)
print(f"  {'Propellant cap':<22}{'Feas. combos':>14}{'Op. years':>12}"
      f"{'New-12 serviced':>18}")
print("-"*72)
results = {}
for label, cap in caps:
    n_feas = int((min_prop <= cap).sum())
    cohorts = pack_cohorts(cap)
    serviced = set()
    for ci in cohorts:
        for d in best_seq[ci]:
            if IS_NEW12[d]:
                serviced.add(int(IDS[d]))
    results[label] = dict(cap=cap, n_feas=n_feas, years=len(cohorts),
                          cohorts=cohorts, new12=len(serviced))
    print(f"  {label:<22}{n_feas:>14,}{len(cohorts):>12}{len(serviced):>14} / 12")
print("="*72)

# Detail: the packed operational-year schedule at the design cap
design_cohorts = results['Design (1019.3 kg)']['cohorts']
print(f"\n  Packed schedule at design cap (1019.3 kg) — "
      f"{len(design_cohorts)} operational years:")
for yr, ci in enumerate(design_cohorts, 1):
    seq = best_seq[ci]
    names = ' -> '.join(f"{IDS[d]}" for d in seq)
    tag = '  (incl. new)' if IS_NEW12[seq].any() else ''
    print(f"    Year {yr:>2}: {names:<34} {min_prop[ci]:>7.1f} kg{tag}")

# =============================================================================
# DIAGNOSTIC — which debris drive infeasibility
#   For each object, look at every combination that CONTAINS it and measure the
#   share whose best ordering still exceeds the cap.  Objects well above the
#   fleet baseline are the ones that make campaigns infeasible.
# =============================================================================

def infeasibility_by_debris(cap):
    base = float((min_prop > cap).mean())
    rows = []
    for d in range(N_DEB):
        in_combo = (_combos == d).any(axis=1)
        mp   = min_prop[in_combo]
        frac = float((mp > cap).mean())
        rows.append((d, int(in_combo.sum()), frac, float(mp.mean())))
    rows.sort(key=lambda r: -r[2])   # worst offenders first
    return base, rows

_cap = MS_PROP_DESIGN
base_rate, diag_rows = infeasibility_by_debris(_cap)

print(f"\n  Per-debris infeasibility at design cap ({_cap:.0f} kg) — "
      f"fleet baseline {base_rate*100:.1f}% of combinations infeasible")
print(f"  {'NORAD':>6}  {'name':<24}{'RAAN':>7}{'inc':>6}{'infeas%':>9}"
      f"{'mean kg':>9}  flag")
print("  " + "-"*72)
for d, n, frac, mean_p in diag_rows:
    new = 'NEW' if IS_NEW12[d] else ''
    hot = '  <== driver' if frac > base_rate*1.5 else ''
    print(f"  {IDS[d]:>6}  {NAMES[d][:23]:<24}{RAAN[d]:>7.1f}{INC[d]:>6.1f}"
          f"{frac*100:>8.1f}%{mean_p:>9.0f}  {new:<3}{hot}")

# =============================================================================
# PLOT 2 — per-debris infeasibility share (the "what blocks campaigns" chart)
# =============================================================================

from matplotlib.patches import Patch
order  = [r[0] for r in diag_rows]
fracs  = np.array([r[2] for r in diag_rows]) * 100.0
labels = [f"{IDS[d]}  {NAMES[d][:18]}" for d in order]
colors = ['tab:orange' if IS_NEW12[d] else 'tab:blue' for d in order]

fig2, ax2 = plt.subplots(figsize=(10, 9))
fig2.patch.set_facecolor('white'); ax2.set_facecolor('white')
ypos = np.arange(N_DEB)
ax2.barh(ypos, fracs, color=colors, edgecolor='white', linewidth=0.4)
ax2.axvline(base_rate*100, color='tab:red', ls='--', lw=1.5)
ax2.set_yticks(ypos); ax2.set_yticklabels(labels, fontsize=8)
ax2.invert_yaxis()
ax2.set_xlabel('Share of campaigns over design cap (1019.3 kg) [%]')
ax2.set_xlim(0, max(fracs.max()*1.05, base_rate*100*1.2))
ax2.legend(handles=[Patch(color='tab:blue',  label='Original 21'),
                    Patch(color='tab:orange', label='New 12'),
                    plt.Line2D([0],[0], color='tab:red', ls='--',
                               label=f'Fleet baseline {base_rate*100:.1f}%')],
           fontsize=9, loc='lower right')
ax2.grid(True, axis='x', lw=0.5, alpha=0.5)

SAVE2 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'figures', 'roi_raan_feasibility', 'roi_debris_infeasibility.png')
plt.tight_layout()
plt.savefig(SAVE2, dpi=150, bbox_inches='tight')
print(f"\n  Diagnostic plot saved -> {SAVE2}")

# =============================================================================
# PLOT  — propellant distribution of campaigns vs the two design caps
# =============================================================================

plt.style.use('default')
fig, ax = plt.subplots(figsize=(11, 6))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

prop_tf = min_prop[time_feasible]
hi = max(MS_PROP_MARGIN*1.05, prop_tf.max()*1.02)
bins = np.linspace(prop_tf.min(), hi, 70)

ax.hist(prop_tf, bins=bins, color='tab:blue', alpha=0.75,
        edgecolor='white', linewidth=0.3,
        label=f'Time-feasible campaigns (n = {prop_tf.size:,})')

n_design = int((min_prop <= MS_PROP_DESIGN).sum())
n_margin = int((min_prop <= MS_PROP_MARGIN).sum())
ax.axvline(MS_PROP_DESIGN, color='tab:green', lw=2.0, ls='--',
           label=f'Design cap 1019.3 kg  ->  {n_design:,} feasible  '
                 f'({results["Design (1019.3 kg)"]["years"]} op. years)')
ax.axvline(MS_PROP_MARGIN, color='tab:orange', lw=2.0, ls='--',
           label=f'Margin cap 1172.2 kg  ->  {n_margin:,} feasible  '
                 f'({results["Margin (1172.2 kg)"]["years"]} op. years)')

ax.set_xlabel('MS orbital propellant per campaign [kg]')
ax.set_ylabel('Number of 5-target campaigns')
ax.set_xlim(prop_tf.min(), hi)
ax.legend(fontsize=9)
ax.grid(True, lw=0.5, alpha=0.5)

SAVE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'figures', 'roi_raan_feasibility', 'roi_raan_feasibility.png')
os.makedirs(os.path.dirname(SAVE), exist_ok=True)
plt.tight_layout()
plt.savefig(SAVE, dpi=150, bbox_inches='tight')
print(f"\n  Plot saved -> {SAVE}")
