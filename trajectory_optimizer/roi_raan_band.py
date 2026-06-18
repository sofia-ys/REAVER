"""
REAVER — ROI Analysis 1: Hub Relocation to a Fresh RAAN Band
============================================================
Scenario (ROI Section 19.4 mission-extension case):
    The Recycling Hub relocates (with the mothership onboard) to a new debris
    cluster in the 0-30 deg RAAN band, re-optimises its parked RAAN inside that
    band, and operates there.  How many 5-target campaigns are feasible, and how
    many operational years does that unlock?

Two things this script makes explicit, because the pitch needs them:
  1. A side-by-side: the same 0-30 deg debris pool serviced FROM THE CURRENT
     64 deg hub (no move) vs. from the relocated, in-band-optimised hub.  This
     proves whether relocating actually buys anything.
  2. The one-time plane-change dV the hub pays to relocate from 64 deg to its
     parked in-band RAAN (single-burn node rotation at the hub orbit).

Feasibility caps act on per-campaign MS *orbital* propellant (transfers only):
    Design 1019.3 kg (RCS reserve intact) / Margin 1172.2 kg (into RCS reserve).
A campaign = one operational year; capacity = disjoint feasible cohorts.

Self-contained: imports only scalar constants from config; modifies nothing.
"""

import os
import time
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config import (MU, DAY, D2R, MS_DRY, MS_VEX, TUG_DRY, TUG_VEX, TUG_THR,
                    TUG_PROP_MARGIN, N_PHASE_REV, T_OPS, MAX_DAYS,
                    RH_SMA, RH_INC, MS_THR, MS_BURN_S, MS_WET_ESTIMATE)

# =============================================================================
# ANALYSIS CONFIGURATION
# =============================================================================

RH_RAAN_CURRENT = 64.0               # deg  present hub position (no-move baseline)
BAND_LO, BAND_HI, BAND_STEP = 0.0, 30.0, 2.0   # deg  in-band sweep for relocated hub
MS_PROP_DESIGN  = 1019.3             # kg   design MS orbital propellant
MS_PROP_MARGIN  = 1172.2             # kg   orbital + 15% RCS margin

# =============================================================================
# DEBRIS CATALOGUE — 0-30 deg RAAN cluster (23 objects)
#   columns used: norad_id, name, mass_kg, sma_km, inc_deg, raan_deg
# =============================================================================

_RAW = [
    (270, 'INTELSAT 602',          1780.97, 42522.783, 14.3090,  7.7100),
    (278, 'INTELSAT 603',          1780.97, 42387.086, 14.0019,  9.4923),
    (283, 'INTELSAT 604',          1780.97, 42688.424, 14.6045, 10.1280),
    (293, 'SBS 6',                 1467.79, 42527.174, 12.6059, 28.9768),
    (311, 'ASTRA 1B',              1548.07, 42665.944, 13.3408, 22.3882),
    (319, 'INTELSAT 605',          2502.48, 42430.955, 13.5487, 17.7038),
    (326, 'INTELSAT 601',          2502.48, 42357.291, 13.3671, 17.0986),
    (362, 'ASTRA 1C',              1684.84, 42570.514, 12.8225, 25.7583),
    (396, 'DIRECTV 2 (DBS 2)',     1711.60, 42630.245, 12.8820, 26.6313),
    (409, 'ASTRA 1D',              1684.84, 42539.085, 12.2616, 28.8908),
    (421, 'AMSC 1',                1635.28, 42551.507, 13.3059, 16.6547),
    (433, 'TELSTAR 4 (402R)',      2078.30, 42161.446, 13.5677, 12.8568),
    (473, 'INTELSAT 26 (JCSAT 4)', 1803.77, 42551.764, 12.2481, 29.6978),
    (474, 'USA 130',               2378.59, 42517.113, 14.4280,  3.9730),
    (477, 'DIRECTV 6 (TEMPO 2)',   2418.24, 42515.928, 13.1017, 22.9139),
    (500, 'INTELSAT 804',          2060.46, 42166.224, 13.1296, 17.7970),
    (506, 'UFO 8 (USA 138)',       1529.24, 42164.201, 11.1053, 22.9055),
    (521, 'UFO 9 (USA 140)',       1529.24, 42762.858, 11.2807, 21.4172),
    (530, 'SKYNET 4E',             1490.00, 42516.912, 12.0550,  9.2220),
    (560, 'USA 149',               2386.00, 42166.218, 13.1551, 12.0135),
    (593, 'SKYNET 4F',             1489.00, 42482.488, 12.1030, 18.6123),
    (604, 'USA 159',               2386.00, 42467.074, 12.7621, 16.5987),
    (664, 'USA 176',               2361.00, 42125.014, 12.0340, 26.2386),
]

N_DEB = len(_RAW)
IDS   = np.array([r[0] for r in _RAW])
NAMES = [r[1] for r in _RAW]
MASS  = np.array([r[2] for r in _RAW])
SMA   = np.array([r[3] for r in _RAW]) * 1e3
INC   = np.array([r[4] for r in _RAW])
RAAN  = np.array([r[5] for r in _RAW])

RH_IDX   = N_DEB
SMA_ALL  = np.append(SMA, RH_SMA)
INC_ALL  = np.append(INC, RH_INC)

# =============================================================================
# PHYSICS  (identical model to roi_raan_feasibility.py / reaver_core)
# =============================================================================

def build_tables(rh_raan_deg):
    raan_all = np.append(RAAN, rh_raan_deg)
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
    if dv_req <= 0.0:
        return 0.0
    dv_per_fire = (MS_THR/mass_kg)*MS_BURN_S
    n_fires  = int(np.ceil(dv_req/dv_per_fire))
    n_orbits = int(np.ceil(n_fires/2))
    return n_orbits*2*np.pi*np.sqrt(sma_m**3/MU)/DAY


def build_finite_leg_table(dv1, dv2, dv_ph, t_tr, t_ph, mass_kg):
    N = N_DEB + 1
    t_fin = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            m1 = mass_kg
            m2 = m1*np.exp(-dv1[i, j]/MS_VEX)
            m3 = m2*np.exp(-dv2[i, j]/MS_VEX)
            t_fin[i, j] = (finite_burn_time(dv1[i, j], m1, SMA_ALL[i]) + t_tr[i, j]
                         + finite_burn_time(dv2[i, j], m2, SMA_ALL[j])
                         + finite_burn_time(dv_ph[i, j], m3, SMA_ALL[j]) + t_ph[i, j])
    return t_fin


def compute_tug_spirals(rh_raan_deg):
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


def relocation_dv(raan_from, raan_to, inc_deg=RH_INC):
    """One-time single-burn plane-change dV [m/s] to rotate the hub node."""
    i  = inc_deg*D2R
    dO = (raan_to - raan_from)*D2R
    cth = np.cos(i)**2 + np.sin(i)**2*np.cos(dO)
    th  = np.arccos(np.clip(cth, -1.0, 1.0))
    return 2.0*np.sqrt(MU/RH_SMA)*np.sin(th/2.0)

# =============================================================================
# CAMPAIGN EVALUATOR  (best time-feasible ordering per 5-target combination)
# =============================================================================

_combos = np.array(list(itertools.combinations(range(N_DEB), 5)), dtype=np.int32)
_perms  = np.array(list(itertools.permutations(range(5))),        dtype=np.int32)
N_COMBOS, N_PERMS = len(_combos), len(_perms)
T_OPS_ARR = np.array([T_OPS]*5 + [0.0])


def evaluate_pool(rh_raan_deg):
    """Return per-combo min orbital propellant (best time-feasible ordering)."""
    dv1, dv2, dv_ph, t_tr, t_ph = build_tables(rh_raan_deg)
    DV_LEG = dv1 + dv2 + dv_ph
    T_LEG  = build_finite_leg_table(dv1, dv2, dv_ph, t_tr, t_ph, MS_WET_ESTIMATE)
    TUG_TIME, TUG_MPROP = compute_tug_spirals(rh_raan_deg)

    min_prop = np.full(N_COMBOS, np.inf)
    best_seq = np.full((N_COMBOS, 5), -1, dtype=np.int32)
    CHUNK = 6000
    for s in range(0, N_COMBOS, CHUNK):
        cc = _combos[s:s+CHUNK]; c = len(cc)
        seqs = cc[:, _perms].reshape(c*N_PERMS, 5); ns = len(seqs)
        frm = np.empty((ns, 6), np.int32); to = np.empty((ns, 6), np.int32)
        frm[:, 0] = RH_IDX; frm[:, 1:5] = seqs[:, :4]; frm[:, 5] = seqs[:, 4]
        to[:, :5] = seqs;   to[:, 5]    = RH_IDX
        dv_legs = DV_LEG[frm, to]; t_legs = T_LEG[frm, to]
        tug_mwet = (TUG_DRY + tug_loads(TUG_MPROP[seqs])).astype(np.float64)
        m_req = np.full(ns, MS_DRY)*np.exp(dv_legs[:, 5]/MS_VEX)
        for leg in range(4, -1, -1):
            m_req += tug_mwet[:, leg]; m_req *= np.exp(dv_legs[:, leg]/MS_VEX)
        ms_prop = m_req - MS_DRY - tug_mwet.sum(axis=1)
        cum_time   = np.cumsum(t_legs + T_OPS_ARR[None, :], axis=1)
        tug_arrive = cum_time[:, :5] + TUG_TIME[seqs]
        mission_day = (np.maximum(cum_time[:, 5:6], tug_arrive) + T_OPS).max(axis=1)
        mp = np.where(mission_day <= MAX_DAYS, ms_prop, np.inf).reshape(c, N_PERMS)
        bp = mp.argmin(axis=1)
        min_prop[s:s+c] = mp[np.arange(c), bp]
        best_seq[s:s+c] = seqs.reshape(c, N_PERMS, 5)[np.arange(c), bp]
    return min_prop, best_seq


def pack_cohorts(min_prop, best_seq, cap):
    idx   = np.where(min_prop <= cap)[0]
    order = idx[np.argsort(min_prop[idx])]
    used  = np.zeros(N_DEB, dtype=bool); chosen = []
    for ci in order:
        s = best_seq[ci]
        if not used[s].any():
            used[s] = True; chosen.append(ci)
    return chosen

# =============================================================================
# 1) IN-BAND RAAN SWEEP  (where to park the relocated hub)
# =============================================================================

band = np.arange(BAND_LO, BAND_HI + 1e-9, BAND_STEP)
print(f"  Pool: {N_DEB} objects, RAAN {RAAN.min():.0f}-{RAAN.max():.0f} deg, "
      f"inc {INC.min():.1f}-{INC.max():.1f} deg")
print(f"  Sweeping hub RAAN {BAND_LO:.0f}-{BAND_HI:.0f} deg "
      f"({len(band)} points), C({N_DEB},5)={N_COMBOS:,} campaigns each...")

t0 = time.time()
sweep_feas, sweep_years = [], []
for raan in band:
    mp, bs = evaluate_pool(float(raan))
    sweep_feas.append(int((mp <= MS_PROP_DESIGN).sum()))
    sweep_years.append(len(pack_cohorts(mp, bs, MS_PROP_DESIGN)))
sweep_feas  = np.array(sweep_feas); sweep_years = np.array(sweep_years)
print(f"  Sweep done ({time.time()-t0:.1f}s)")

# Optimum in-band: most operational years, tie-broken by most feasible campaigns
best_i   = int(np.lexsort((sweep_feas, sweep_years))[-1])
opt_raan = float(band[best_i])

# =============================================================================
# 2) EVALUATE THE TWO SCENARIOS  (no-move @64 deg  vs  relocated @opt)
# =============================================================================

mp_cur, bs_cur = evaluate_pool(RH_RAAN_CURRENT)
mp_opt, bs_opt = evaluate_pool(opt_raan)
reloc_dv = relocation_dv(RH_RAAN_CURRENT, opt_raan)

def scen(mp, bs):
    out = {}
    for label, cap in (('design', MS_PROP_DESIGN), ('margin', MS_PROP_MARGIN)):
        coh = pack_cohorts(mp, bs, cap)
        out[label] = (int((mp <= cap).sum()), len(coh), coh)
    return out

S_cur, S_opt = scen(mp_cur, bs_cur), scen(mp_opt, bs_opt)
theo_years = N_DEB // 5

print("\n" + "="*74)
print(f"  ANALYSIS 1 — 0-30 deg RAAN POOL  ({N_DEB} objects, ceiling "
      f"{theo_years} op. years)")
print("="*74)
print(f"  {'Hub scenario':<34}{'Feas. (design)':>16}{'Op. years':>12}")
print("  " + "-"*72)
print(f"  {'Current hub @ 64 deg (NO move)':<34}"
      f"{S_cur['design'][0]:>10,}/{N_COMBOS:,}{S_cur['design'][1]:>12}")
print(f"  {f'Relocated & optimised @ {opt_raan:.0f} deg':<34}"
      f"{S_opt['design'][0]:>10,}/{N_COMBOS:,}{S_opt['design'][1]:>12}")
print("  " + "-"*72)
print(f"  One-time relocation dV (64 -> {opt_raan:.0f} deg, hub @ inc "
      f"{RH_INC:.0f} deg) : {reloc_dv:.0f} m/s")
print(f"  Margin-cap op. years  : no-move {S_cur['margin'][1]}  |  "
      f"relocated {S_opt['margin'][1]}")
print("="*74)

print(f"\n  Packed schedule at relocated hub @ {opt_raan:.0f} deg "
      f"(design cap) — {S_opt['design'][1]} operational years:")
for yr, ci in enumerate(S_opt['design'][2], 1):
    seq = bs_opt[ci]
    print(f"    Year {yr:>2}: {' -> '.join(str(IDS[d]) for d in seq):<32} "
          f"{mp_opt[ci]:>7.1f} kg")

# Per-debris infeasibility at the relocated optimum
base = float((mp_opt > MS_PROP_DESIGN).mean())
rows = sorted(((d, float((mp_opt[(_combos == d).any(1)] > MS_PROP_DESIGN).mean()))
               for d in range(N_DEB)), key=lambda r: -r[1])
print(f"\n  Objects driving infeasibility @ {opt_raan:.0f} deg "
      f"(baseline {base*100:.1f}%):")
for d, frac in rows:
    if frac > base*1.5:
        print(f"    {IDS[d]:>5}  {NAMES[d][:24]:<25} RAAN {RAAN[d]:>5.1f}  "
              f"inc {INC[d]:>5.1f}  ->  {frac*100:.1f}%")

# =============================================================================
# PLOT 1 — in-band sweep: operational years & feasible campaigns vs hub RAAN
# =============================================================================

plt.style.use('default')
fig, ax = plt.subplots(figsize=(11, 6))
fig.patch.set_facecolor('white'); ax.set_facecolor('white')

ax.plot(band, sweep_feas, color='tab:blue', lw=2.0, marker='o', ms=4,
        label='Feasible campaigns @ design cap')
ax.axvline(opt_raan, color='tab:green', ls='--', lw=1.8,
           label=f'In-band optimum = {opt_raan:.0f} deg '
                 f'({S_opt["design"][1]} op. years)')
ax.scatter([opt_raan], [sweep_feas[best_i]], color='tab:green', s=70, zorder=6)
ax.axhline(S_cur['design'][0], color='tab:red', ls=':', lw=1.8,
           label=f'Current hub @ 64 deg, no move ({S_cur["design"][1]} op. years)')

ax.set_xlabel('Parked hub RAAN in 0-30 deg band [deg]')
ax.set_ylabel('Feasible 5-target campaigns @ 1019.3 kg')
ax.set_xlim(BAND_LO, BAND_HI)
ax.legend(fontsize=9, loc='best')
ax.grid(True, lw=0.5, alpha=0.5)

SAVE1 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'figures', 'roi_raan_band', 'roi_band_sweep.png')
os.makedirs(os.path.dirname(SAVE1), exist_ok=True)
plt.tight_layout(); plt.savefig(SAVE1, dpi=150, bbox_inches='tight')
print(f"\n  Sweep plot saved -> {SAVE1}")

# =============================================================================
# PLOT 2 — campaign propellant distribution at the relocated optimum
# =============================================================================

fig2, ax2 = plt.subplots(figsize=(11, 6))
fig2.patch.set_facecolor('white'); ax2.set_facecolor('white')
prop = mp_opt[np.isfinite(mp_opt)]
hi   = max(MS_PROP_MARGIN*1.05, prop.max()*1.02)
ax2.hist(prop, bins=np.linspace(prop.min(), hi, 60), color='tab:blue',
         alpha=0.75, edgecolor='white', linewidth=0.3,
         label=f'Campaigns @ hub {opt_raan:.0f} deg (n = {prop.size:,})')
ax2.axvline(MS_PROP_DESIGN, color='tab:green', ls='--', lw=2.0,
            label=f'Design cap 1019.3 kg -> {S_opt["design"][0]:,} feasible')
ax2.axvline(MS_PROP_MARGIN, color='tab:orange', ls='--', lw=2.0,
            label=f'Margin cap 1172.2 kg -> {S_opt["margin"][0]:,} feasible')
ax2.set_xlabel('MS orbital propellant per campaign [kg]')
ax2.set_ylabel('Number of 5-target campaigns')
ax2.set_xlim(prop.min(), hi)
ax2.legend(fontsize=9); ax2.grid(True, lw=0.5, alpha=0.5)

SAVE2 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'figures', 'roi_raan_band', 'roi_band_propellant.png')
plt.tight_layout(); plt.savefig(SAVE2, dpi=150, bbox_inches='tight')
print(f"  Propellant plot saved -> {SAVE2}")
