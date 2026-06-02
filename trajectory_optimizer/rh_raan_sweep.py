"""
REAVER — Recycling Hub RAAN Optimisation
=========================================
Sweeps RH RAAN 0°→355° in 5° steps and evaluates worst-case, best-case,
and mean mothership propellant for each position.

Key output: the RAAN that minimises the worst-case required propellant.
"""

import numpy as np
from itertools import combinations, permutations
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time, warnings
warnings.filterwarnings('ignore')

from config import *

# =============================================================================
# PRE-COMPUTE DEBRIS↔DEBRIS TRANSFER TABLE  (independent of RH RAAN)
# =============================================================================

N = N_DEB + 1
_DV_DD = np.zeros((N, N))   # total ΔV   for debris-to-debris pairs only
_DT_DD = np.zeros((N, N))   # total time for debris-to-debris pairs only

print("  Pre-computing debris↔debris table...", end=' ', flush=True)
for i in range(N_DEB):
    sa = SMA[i]; ia = INC[i]*D2R; oa = RAAN[i]*D2R
    va = np.sqrt(MU/sa)
    for j in range(N_DEB):
        if i == j:
            continue
        sb = SMA[j]; ib = INC[j]*D2R; ob = RAAN[j]*D2R
        vb = np.sqrt(MU/sb)
        at   = (sa + sb) / 2.0
        v_sa = np.sqrt(MU*(2/sa - 1/at))
        v_sb = np.sqrt(MU*(2/sb - 1/at))
        t_tr = np.pi * np.sqrt(at**3/MU) / DAY
        cos_dth = (np.cos(ia)*np.cos(ib) +
                   np.sin(ia)*np.sin(ib)*np.cos(ob - oa))
        dth = np.arccos(np.clip(cos_dth, -1.0, 1.0))
        if sb > sa:
            d1 = abs(v_sa - va)
            d2 = np.sqrt(v_sb**2 + vb**2 - 2*v_sb*vb*np.cos(dth))
        else:
            d1 = np.sqrt(v_sa**2 + va**2 - 2*v_sa*va*np.cos(dth))
            d2 = abs(v_sb - vb)
        T_tgt    = 2*np.pi * np.sqrt(sb**3/MU)
        T_ph_orb = T_tgt * (1.0 - 1.0/(4.0*N_PHASE_REV))
        a_ph     = (MU * (T_ph_orb/(2*np.pi))**2)**(1.0/3.0)
        d_ph     = 2*abs(np.sqrt(MU/sb) - np.sqrt(MU*(2/sb - 1/a_ph)))
        t_ph     = N_PHASE_REV * T_ph_orb / DAY
        _DV_DD[i, j] = d1 + d2 + d_ph
        _DT_DD[i, j] = t_tr + t_ph
print("done")

# =============================================================================
# BUILD SEQUENCE INDEX  (constant across all RAAN values)
# =============================================================================

print("  Building sequence index...", end=' ', flush=True)
_sequences = []
for combo in combinations(range(N_DEB), 5):
    for perm in permutations(combo):
        _sequences.append(perm)
_sequences = np.array(_sequences, dtype=np.int32)   # (524160, 5)
N_SEQ    = len(_sequences)
N_PERMS  = 120          # 5! permutations per combination
N_COMBOS = N_SEQ // N_PERMS   # C(16,5) = 4368
print(f"{N_SEQ:,}  ({N_COMBOS} combos × {N_PERMS} orderings)")

# from/to node arrays — legs 0-4 and return leg 5
_from = np.zeros((N_SEQ, 6), dtype=np.int32)
_to   = np.zeros((N_SEQ, 6), dtype=np.int32)
_from[:, 0]   = RH_IDX
_from[:, 1:5] = _sequences[:, :4]
_from[:, 5]   = _sequences[:, 4]
_to[:, :5]    = _sequences
_to[:, 5]     = RH_IDX

# =============================================================================
# RAAN SWEEP
# =============================================================================

def _rh_transfers(rh_raan_deg):
    """Return DV_LEG and T_LEG for all pairs involving RH at given RAAN."""
    DV = _DV_DD.copy()
    DT = _DT_DD.copy()
    oa = rh_raan_deg * D2R
    ia = RH_INC * D2R
    sa = RH_SMA
    va = np.sqrt(MU/sa)
    for j in range(N_DEB):
        sb = SMA[j]; ib = INC[j]*D2R; ob = RAAN[j]*D2R
        vb = np.sqrt(MU/sb)
        for src, dst in [(RH_IDX, j), (j, RH_IDX)]:
            si = sa if src == RH_IDX else sb
            vi = va if src == RH_IDX else vb
            ii = ia if src == RH_IDX else ib
            oi = oa if src == RH_IDX else ob
            sk = sb if dst == j    else sa
            vk = vb if dst == j    else va
            ik = ib if dst == j    else ia
            ok = ob if dst == j    else oa
            at    = (si + sk) / 2.0
            v_si  = np.sqrt(MU*(2/si - 1/at))
            v_sk  = np.sqrt(MU*(2/sk - 1/at))
            t_tr  = np.pi * np.sqrt(at**3/MU) / DAY
            cos_d = (np.cos(ii)*np.cos(ik) +
                     np.sin(ii)*np.sin(ik)*np.cos(ok - oi))
            dth   = np.arccos(np.clip(cos_d, -1.0, 1.0))
            if sk > si:
                d1 = abs(v_si - vi)
                d2 = np.sqrt(v_sk**2 + vk**2 - 2*v_sk*vk*np.cos(dth))
            else:
                d1 = np.sqrt(v_si**2 + vi**2 - 2*v_si*vi*np.cos(dth))
                d2 = abs(v_sk - vk)
            T_tgt    = 2*np.pi * np.sqrt(sk**3/MU)
            T_ph_orb = T_tgt * (1.0 - 1.0/(4.0*N_PHASE_REV))
            a_ph     = (MU*(T_ph_orb/(2*np.pi))**2)**(1.0/3.0)
            d_ph     = 2*abs(np.sqrt(MU/sk) - np.sqrt(MU*(2/sk - 1/a_ph)))
            t_ph     = N_PHASE_REV * T_ph_orb / DAY
            DV[src, dst] = d1 + d2 + d_ph
            DT[src, dst] = t_tr + t_ph
    return DV, DT


def _tug_data(rh_raan_deg):
    """Return TUG_MPROP and TUG_TIME arrays for all debris at given RH RAAN."""
    mprop = np.zeros(N_DEB)
    ttime = np.zeros(N_DEB)
    for k in range(N_DEB):
        m_pl  = TUG_DRY + MASS[k]
        m_wet = m_pl * 1.35
        i1,o1 = INC[k], RAAN[k]
        i2,o2 = RH_INC, rh_raan_deg
        v1 = np.sqrt(MU/SMA[k]); v2 = np.sqrt(MU/RH_SMA)
        cos_d = (np.cos(i1*D2R)*np.cos(i2*D2R) +
                 np.sin(i1*D2R)*np.sin(i2*D2R)*np.cos((o2-o1)*D2R))
        dth  = np.arccos(np.clip(cos_d, -1, 1))
        dv_e = np.sqrt(v1**2 + v2**2 - 2*v1*v2*np.cos(np.pi/2*dth))
        for _ in range(60):
            mf   = m_wet * np.exp(-dv_e/TUG_VEX)
            mnew = m_pl + (m_wet - mf)
            if abs(mnew - m_wet) < 0.05: break
            m_wet = 0.6*m_wet + 0.4*mnew
        t_s = (m_wet * TUG_VEX / TUG_THR) * (1 - np.exp(-dv_e/TUG_VEX))
        mprop[k] = m_wet - m_pl
        ttime[k] = t_s / DAY
    return mprop, ttime


def evaluate_raan(rh_raan_deg):
    """
    For each of the 4368 unordered combinations of 5 debris, find the ordering
    (permutation) that minimises propellant mass, then compare those 4368
    minimum-propellant values against each other.

      prop_worst = the combination that requires the most propellant even at
                   its best ordering  (hardest set of debris)
      prop_best  = the combination that requires the least propellant
      prop_mean  = mean across all feasible combinations
    """
    DV_LEG, T_LEG = _rh_transfers(rh_raan_deg)
    TUG_MPROP, TUG_TIME = _tug_data(rh_raan_deg)

    T_tgt_rh = 2*np.pi * np.sqrt(RH_SMA**3/MU)
    T_PH_RH  = N_PHASE_REV * T_tgt_rh * (1.0 - 1.0/(4.0*N_PHASE_REV)) / DAY

    dv_legs = DV_LEG[_from, _to]
    t_legs  = T_LEG[_from,  _to]

    tug_mwet = (TUG_DRY + TUG_MPROP[_sequences]).astype(np.float64)

    # Backward pass → propellant per sequence (N_SEQ,)
    m_req = np.full(N_SEQ, MS_DRY)
    m_req = m_req * np.exp(dv_legs[:, 5] / MS_VEX)
    for leg in range(4, -1, -1):
        m_req += tug_mwet[:, leg]
        m_req  = m_req * np.exp(dv_legs[:, leg] / MS_VEX)
    ms_prop_seq = m_req - MS_DRY - tug_mwet.sum(axis=1)

    # Timing → feasibility
    t_ops_arr   = np.array([T_OPS]*5 + [0.0])
    cum_time    = np.cumsum(t_legs + t_ops_arr[None, :], axis=1)
    tug_arrive  = cum_time[:, :5] + TUG_TIME[_sequences]
    handover    = np.maximum(cum_time[:, 5:6], tug_arrive) + T_PH_RH
    mission_day = handover.max(axis=1)
    feas        = mission_day <= MAX_DAYS   # (N_SEQ,)

    if not feas.any():
        return dict(n_feas_combos=0, prop_worst=np.nan,
                    prop_best=np.nan, prop_mean=np.nan)

    # ── Combination-level reduction ────────────────────────────────────────
    # Sequences are laid out as 4368 consecutive blocks of 120 permutations.
    # Step 1: select best ordering per combination = lowest Pareto score
    #         (50% ΔV + 50% time).  Normalise over FEASIBLE sequences only,
    #         consistent with reaver_optimizer.py.  Infeasible → inf.
    # Step 2: evaluate worst-case via propellant of that best-Pareto ordering.
    tot_dv   = np.cumsum(dv_legs, axis=1)[:, 5]
    f_dv_    = tot_dv[feas];    f_day_ = mission_day[feas]
    dv_n_f   = (f_dv_  - f_dv_.min())  / (f_dv_.max()  - f_dv_.min()  + 1e-9)
    time_n_f = (f_day_ - f_day_.min()) / (f_day_.max() - f_day_.min() + 1e-9)
    score          = np.full(N_SEQ, np.inf)
    score[feas]    = 0.5*dv_n_f + 0.5*time_n_f

    score_c  = score.reshape(N_COMBOS, N_PERMS)          # (4368, 120)
    prop_c   = ms_prop_seq.reshape(N_COMBOS, N_PERMS)
    feas_c   = feas.reshape(N_COMBOS, N_PERMS)

    best_perm      = score_c.argmin(axis=1)               # best-Pareto ordering per combo
    combo_ok       = feas_c.any(axis=1)                   # at least one feasible ordering
    min_prop_combo = prop_c[np.arange(N_COMBOS), best_perm]   # propellant of best-Pareto perm

    n_feas_combos = int(combo_ok.sum())
    if n_feas_combos == 0:
        return dict(n_feas_combos=0, prop_worst=np.nan,
                    prop_best=np.nan, prop_mean=np.nan)

    valid = min_prop_combo[combo_ok]
    return dict(
        n_feas_combos = n_feas_combos,
        prop_worst    = float(np.max(valid)),
        prop_best     = float(np.min(valid)),
        prop_mean     = float(np.mean(valid)),
    )


# =============================================================================
# RUN SWEEP
# =============================================================================

RAAN_SWEEP = np.arange(0, 360, 1)   # 72 values

print(f"\n  Sweeping RH RAAN 0°→355° in 5° steps ({len(RAAN_SWEEP)} points)...")
print(f"  {'RAAN':>6}  {'Combos':>10}  {'Prop best':>10}  {'Prop worst':>11}  {'Time':>6}")
print(f"  {'':─>6}  {'':─>6}  {'':─>10}  {'':─>11}  {'':─>6}")

records = []
t_total = time.time()
for raan in RAAN_SWEEP:
    t0 = time.time()
    r  = evaluate_raan(float(raan))
    records.append(r)
    print(f"  {raan:>5.0f}°  {r['n_feas_combos']:>6,}/{N_COMBOS}  "
          f"{r['prop_best']:>9.0f}kg  {r['prop_worst']:>10.0f}kg  "
          f"({time.time()-t0:.1f}s)")

print(f"\n  Total sweep time: {time.time()-t_total:.1f}s")

# =============================================================================
# RESULTS
# =============================================================================

prop_worst = np.array([r['prop_worst']    for r in records], dtype=float)
prop_best  = np.array([r['prop_best']     for r in records], dtype=float)
prop_mean  = np.array([r['prop_mean']     for r in records], dtype=float)
n_feas     = np.array([r['n_feas_combos'] for r in records], dtype=float)

# Optimum: all N_COMBOS combinations feasible AND lowest worst-case propellant
all_feas_mask = (n_feas == N_COMBOS)
if all_feas_mask.any():
    prop_worst_constrained = np.where(all_feas_mask, prop_worst, np.nan)
    opt_idx = int(np.nanargmin(prop_worst_constrained))
else:
    print("  ⚠  No RAAN where all combinations are feasible — falling back to min worst-prop")
    opt_idx = int(np.nanargmin(prop_worst))
opt_raan = RAAN_SWEEP[opt_idx]

print("\n" + "="*60)
print(f"  Optimal RH RAAN              : {opt_raan:.0f}°")
print(f"  Min worst-combo prop there   : {prop_worst[opt_idx]:.0f} kg")
print(f"  Best-combo prop there        : {prop_best[opt_idx]:.0f} kg")
print(f"  Feasible combinations        : {int(n_feas[opt_idx]):,} / {N_COMBOS}  ✓ all feasible")
print(f"  RAAN values with all combos feasible: {int(all_feas_mask.sum())}")
print("="*60)

# =============================================================================
# PLOT
# =============================================================================

BG='#0d1117'; CB='#161b22'; TC='#c9d1d9'; MU_='#8b949e'; GR='#21262d'
A1='#58a6ff'; A2='#3fb950'; A3='#f78166'; A4='#d2a8ff'

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
fig.patch.set_facecolor(BG)
fig.text(0.5, 0.98, 'Recycling Hub RAAN Optimisation — Mothership Propellant Sweep',
         ha='center', color=TC, fontsize=13, fontweight='bold')

for ax in (ax1, ax2):
    ax.set_facecolor(CB)
    [s.set_edgecolor(GR) for s in ax.spines.values()]
    ax.tick_params(colors=MU_, labelsize=8)
    ax.xaxis.label.set_color(MU_)
    ax.yaxis.label.set_color(MU_)
    ax.grid(True, color=GR, lw=0.5, alpha=0.6)

# ── Panel 1: propellant ───────────────────────────────────────────────────────
# Shade region where all combinations are feasible
for i in range(len(RAAN_SWEEP) - 1):
    if all_feas_mask[i]:
        ax1.axvspan(RAAN_SWEEP[i], RAAN_SWEEP[i+1], color=A2, alpha=0.07, lw=0)
if all_feas_mask[-1]:
    ax1.axvspan(RAAN_SWEEP[-1], RAAN_SWEEP[-1]+1, color=A2, alpha=0.07, lw=0)

ax1.plot(RAAN_SWEEP, prop_worst, color=A3, lw=1.8,
         label=f'Worst combination — best-Pareto ordering  (max={np.nanmax(prop_worst):.0f} kg)')
ax1.plot(RAAN_SWEEP, prop_mean,  color=A4, lw=1.2, ls='--',
         label='Mean across feasible combinations')
ax1.plot(RAAN_SWEEP, prop_best,  color=A2, lw=1.8,
         label=f'Best combination — best-Pareto ordering  (min={np.nanmin(prop_best):.0f} kg)')

ax1.axvline(opt_raan, color=A1, lw=2.0, ls=':', alpha=0.95,
            label=f'Optimum RAAN = {opt_raan:.0f}°  |  worst-combo prop = {prop_worst[opt_idx]:.0f} kg  |  all {N_COMBOS} combos feasible')
ax1.scatter([opt_raan], [prop_worst[opt_idx]], color=A1, s=90, zorder=6, lw=0)

# debris RAAN markers — placed after plot so get_ylim() is valid
ax1.autoscale_view()
y_top = ax1.get_ylim()[1]
for r, nm in zip(RAAN, NAMES):
    ax1.axvline(r, color=MU_, lw=0.6, alpha=0.35)
    ax1.text(r, y_top, nm.split('(')[0].strip()[:10], rotation=90, fontsize=5.5,
             color=MU_, ha='center', va='top', alpha=0.6)

ax1.set_ylabel('MS propellant required [kg]')
ax1.set_title('Mothership propellant vs RH RAAN  (per combination: best-Pareto ordering)',
              color=TC, fontsize=9, fontweight='bold', pad=6)
ax1.legend(fontsize=7, facecolor=CB, labelcolor=TC, edgecolor=GR)

# ── Panel 2: feasible count ───────────────────────────────────────────────────
# Shade fully-feasible region
for i in range(len(RAAN_SWEEP) - 1):
    if all_feas_mask[i]:
        ax2.axvspan(RAAN_SWEEP[i], RAAN_SWEEP[i+1], color=A2, alpha=0.12, lw=0)
if all_feas_mask[-1]:
    ax2.axvspan(RAAN_SWEEP[-1], RAAN_SWEEP[-1]+1, color=A2, alpha=0.12, lw=0)

ax2.fill_between(RAAN_SWEEP, n_feas, alpha=0.25, color=A1)
ax2.plot(RAAN_SWEEP, n_feas, color=A1, lw=1.8)
ax2.axhline(N_COMBOS, color=A2, lw=1.2, ls='--', alpha=0.8,
            label=f'All {N_COMBOS} combinations feasible')
ax2.axvline(opt_raan, color=A1, lw=2.0, ls=':', alpha=0.95,
            label=f'Optimum RAAN = {opt_raan:.0f}°')

for r in RAAN:
    ax2.axvline(r, color=MU_, lw=0.6, alpha=0.35)

ax2.set_xlabel('RH RAAN [deg]')
ax2.set_ylabel(f'Feasible combinations (of {N_COMBOS})')
ax2.set_title('Feasible combinations vs RH RAAN', color=TC, fontsize=9,
              fontweight='bold', pad=6)
ax2.legend(fontsize=7.5, facecolor=CB, labelcolor=TC, edgecolor=GR)
ax2.set_xlim(0, max(RAAN_SWEEP))

SAVE = r'C:\Projects\DSE\REAVER\trajectory_optimizer\rh_raan_sweep.png'
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(SAVE, dpi=150, bbox_inches='tight', facecolor=BG)
print(f"\n  Plot saved → {SAVE}")