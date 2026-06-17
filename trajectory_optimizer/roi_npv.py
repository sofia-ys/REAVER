"""
REAVER — ROI Section 19.4/19.5: Mission-Extension NPV with Both RAAN Analyses
============================================================================
Folds the two feasibility studies into the SME-SMAD discounted-cash-flow model.

Operational-year capacity (5 targets = 1 operational year = +EUR290M net):
  * Analysis 2 (roi_raan_feasibility.py): hub HELD at 64 deg, extended 33-object
    pool -> 6 operational years (1 base + 5 extension).  No relocation cost.
  * Analysis 1 (roi_raan_band.py): hub RELOCATES to the 0-30 deg band (parked at
    14 deg) -> +4 further operational years, after a one-time 314 m/s plane change.

Cash-flow convention follows ROI Section 19.3/19.5:
  - Years 1-9 (2026-2034) development: -27.78 EUR M/yr
  - Year 10 (2035) base mission: +300 (first-year ops inside the sunk EUR250M)
  - Each extension year: +290 net (= +300 revenue - 10 ops)
  - One-time hub relocation propellant cost charged in the first Phase-B year.

NPV per SME-SMAD Eq. 12-2 at i = 18% / 20% / 25%.
Self-contained; imports only G0 from config.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from config import G0

# =============================================================================
# ASSUMPTIONS
# =============================================================================

DEV_COST   = -27.78     # EUR M/yr, development outflow (Years 1-9)
DEV_YEARS  = 9
BASE_NET   = 300.0      # EUR M, Year-10 base mission net cash flow
EXT_NET    = 290.0      # EUR M, each extension year (+300 rev - 10 ops)
RATES      = [0.18, 0.20, 0.25]

PHASE_A_EXT = 5         # extension years at hub @ 64 deg   (Analysis 2: 6 op yrs total)
PHASE_B_EXT = 4         # extension years after relocation  (Analysis 1)

# Hub relocation (Analysis 1)
ISS_MASS   = 419725.0   # kg, ISS reference mass
RH_WET     = ISS_MASS / 2.0     # kg, assumed hub wet mass (user estimate)
RH_ISP     = 250.0      # s
RELOC_DV   = 314.0      # m/s, one-time 64 -> 14 deg plane change (roi_raan_band.py)
COST_PER_KG_NOMINAL = 20000.0   # EUR/kg, propellant delivered to super-GEO (assumption)

# =============================================================================
# RELOCATION PROPELLANT  (Tsiolkovsky)
# =============================================================================

ve        = RH_ISP * G0
reloc_kg  = RH_WET * (1.0 - np.exp(-RELOC_DV / ve))      # propellant burned [kg]

def reloc_cost_MEUR(cost_per_kg):
    return reloc_kg * cost_per_kg / 1e6                  # EUR M

# =============================================================================
# CASH-FLOW BUILDERS & NPV
# =============================================================================

def cashflow(n_ext_phaseA, n_ext_phaseB=0, reloc_MEUR=0.0):
    cf  = [DEV_COST]*DEV_YEARS          # Years 1-9
    cf += [BASE_NET]                    # Year 10 base
    cf += [EXT_NET]*n_ext_phaseA        # Phase-A extension years
    if n_ext_phaseB > 0:
        cf += [EXT_NET - reloc_MEUR]    # first Phase-B year carries relocation cost
        cf += [EXT_NET]*(n_ext_phaseB-1)
    return np.array(cf)

def npv(cf, i):
    t = np.arange(1, len(cf)+1)
    return float(np.sum(cf / (1.0+i)**t))

# Scenarios -------------------------------------------------------------------
SCEN = {
    'Doc baseline (5 op yrs, no studies)' : cashflow(4),
    'Analysis 2 only  (6 op yrs @64deg)'  : cashflow(PHASE_A_EXT),
    f'Analysis 1+2 (10 op yrs, relocate)' : cashflow(PHASE_A_EXT, PHASE_B_EXT,
                                                     reloc_cost_MEUR(COST_PER_KG_NOMINAL)),
}

# =============================================================================
# REPORT
# =============================================================================

print("="*78)
print("  HUB RELOCATION PROPELLANT (Analysis 1)")
print("="*78)
print(f"  Hub wet mass (half ISS)     : {RH_WET:>12,.0f} kg")
print(f"  Relocation dV (64->14 deg)  : {RELOC_DV:>12.0f} m/s")
print(f"  Exhaust velocity (Isp 250s) : {ve:>12.1f} m/s")
print(f"  Propellant burned           : {reloc_kg:>12,.0f} kg  ({reloc_kg/1000:.1f} t)")
print(f"  Cost @ {COST_PER_KG_NOMINAL:,.0f} EUR/kg      : {reloc_cost_MEUR(COST_PER_KG_NOMINAL):>12,.0f} EUR M")

print("\n" + "="*78)
print("  NET PRESENT VALUE BY SCENARIO  (EUR M, SME-SMAD Eq. 12-2)")
print("="*78)
print(f"  {'Scenario':<38}{'NPV@18%':>11}{'NPV@20%':>11}{'NPV@25%':>11}")
print("  " + "-"*74)
for name, cf in SCEN.items():
    vals = [npv(cf, i) for i in RATES]
    print(f"  {name:<38}" + "".join(f"{v:>11.1f}" for v in vals))
print("  " + "-"*74)

# Incremental NPV deltas (vs doc baseline) at 20%
base20 = npv(SCEN['Doc baseline (5 op yrs, no studies)'], 0.20)
a2_20  = npv(SCEN['Analysis 2 only  (6 op yrs @64deg)'], 0.20)
full20 = npv(list(SCEN.values())[2], 0.20)
print(f"\n  At i = 20%:")
print(f"    Analysis 2 adds  : {a2_20-base20:+.1f} EUR M  (1 free extension year @64 deg)")
print(f"    Analysis 1 adds  : {full20-a2_20:+.1f} EUR M  (4 yrs - relocation cost "
      f"@ {COST_PER_KG_NOMINAL/1000:.0f}k EUR/kg)")

# Breakeven propellant cost for the relocation (Phase-B NPV neutral)
print("\n" + "="*78)
print("  RELOCATION BREAKEVEN — max propellant cost for Analysis 1 to stay NPV-positive")
print("="*78)
for i in RATES:
    a2  = npv(cashflow(PHASE_A_EXT), i)
    # find reloc_MEUR where full NPV == Analysis-2 NPV
    lo, hi = 0.0, 5000.0
    for _ in range(80):
        mid = 0.5*(lo+hi)
        if npv(cashflow(PHASE_A_EXT, PHASE_B_EXT, mid), i) > a2:
            lo = mid
        else:
            hi = mid
    be_MEUR = 0.5*(lo+hi)
    print(f"  i={i*100:.0f}%:  breakeven relocation cost = {be_MEUR:>7.0f} EUR M  "
          f"->  {be_MEUR*1e6/reloc_kg:>8,.0f} EUR/kg")
print("="*78)

# =============================================================================
# PLOT — NPV of the full (relocate) scenario vs propellant unit cost
# =============================================================================

cpk = np.linspace(0, 45000, 200)
plt.style.use('default')
fig, ax = plt.subplots(figsize=(11, 6))
fig.patch.set_facecolor('white'); ax.set_facecolor('white')
colors = {0.18: 'tab:green', 0.20: 'tab:blue', 0.25: 'tab:red'}
for i in RATES:
    npv_full = [npv(cashflow(PHASE_A_EXT, PHASE_B_EXT, reloc_cost_MEUR(c)), i) for c in cpk]
    npv_a2   = npv(cashflow(PHASE_A_EXT), i)
    ax.plot(cpk/1000, npv_full, color=colors[i], lw=2.0,
            label=f'Relocate (10 op yrs), i={i*100:.0f}%')
    ax.axhline(npv_a2, color=colors[i], lw=1.2, ls=':')
ax.axvline(COST_PER_KG_NOMINAL/1000, color='gray', lw=1.5, ls='--',
           label=f'Assumed {COST_PER_KG_NOMINAL/1000:.0f}k EUR/kg')
ax.set_xlabel('Relocation propellant delivered cost [thousand EUR/kg]')
ax.set_ylabel('Extended-mission NPV [EUR M]')
ax.set_xlim(0, 45)
ax.legend(fontsize=8.5, loc='upper right', ncol=2)
ax.grid(True, lw=0.5, alpha=0.5)
ax.text(0.6, npv(cashflow(PHASE_A_EXT), 0.20)+1.5,
        'dotted = Analysis 2 only (no relocation)', fontsize=8, color='dimgray')

SAVE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'figures', 'roi_npv', 'roi_npv_relocation.png')
os.makedirs(os.path.dirname(SAVE), exist_ok=True)
plt.tight_layout(); plt.savefig(SAVE, dpi=150, bbox_inches='tight')
print(f"\n  NPV plot saved -> {SAVE}")
