"""
REAVER Core Physics
====================
Shared functions used by reaver_optimizer.py, nominal_mission.py, and
rh_raan_sweep.py.  Import from here instead of duplicating in each script.

Module-level pre-computation runs once at import time using the default
RH RAAN from config.py, giving the same module-level names (DV1, DV_LEG,
TUG_MPROP, etc.) that the scripts previously defined themselves.
"""

import numpy as np
import time

from config import *

# =============================================================================
# SHARED FUNCTIONS
# =============================================================================

def build_transfer_table(rh_raan_deg=RH_RAAN):
    """
    Build the full 17x17 pairwise transfer table (16 debris + RH).

    Parameters
    ----------
    rh_raan_deg : float
        RH RAAN in degrees.  Defaults to the value in config.py.

    Returns
    -------
    dv1, dv2, dv_ph, t_tr, t_ph : (17,17) arrays
    """
    raan_all = RAAN_ALL.copy()
    raan_all[RH_IDX] = rh_raan_deg

    N = N_DEB + 1
    dv1   = np.zeros((N, N))
    dv2   = np.zeros((N, N))
    dv_ph = np.zeros((N, N))
    t_tr  = np.zeros((N, N))
    t_ph  = np.zeros((N, N))

    for i in range(N):
        sa = SMA_ALL[i]
        ia = INC_ALL[i] * D2R
        oa = raan_all[i] * D2R
        va = np.sqrt(MU / sa)
        for j in range(N):
            if i == j:
                continue
            sb = SMA_ALL[j]
            ib = INC_ALL[j] * D2R
            ob = raan_all[j] * D2R
            vb = np.sqrt(MU / sb)

            at   = (sa + sb) / 2.0
            v_sa = np.sqrt(MU * (2.0/sa - 1.0/at))
            v_sb = np.sqrt(MU * (2.0/sb - 1.0/at))
            t_tr[i, j] = np.pi * np.sqrt(at**3 / MU) / DAY

            cos_dth = (np.cos(ia)*np.cos(ib) +
                       np.sin(ia)*np.sin(ib)*np.cos(ob - oa))
            dth = np.arccos(np.clip(cos_dth, -1.0, 1.0))

            if sb > sa:
                dv1[i, j] = abs(v_sa - va)
                dv2[i, j] = np.sqrt(v_sb**2 + vb**2 - 2*v_sb*vb*np.cos(dth))
            else:
                dv1[i, j] = np.sqrt(v_sa**2 + va**2 - 2*v_sa*va*np.cos(dth))
                dv2[i, j] = abs(v_sb - vb)

            T_tgt      = 2*np.pi * np.sqrt(sb**3 / MU)
            T_ph_orb   = T_tgt * (1.0 - 1.0 / (4.0 * N_PHASE_REV))
            a_ph       = (MU * (T_ph_orb / (2*np.pi))**2) ** (1.0/3.0)
            dv_ph[i,j] = 2*abs(np.sqrt(MU/sb) - np.sqrt(MU*(2/sb - 1/a_ph)))
            t_ph[i, j] = N_PHASE_REV * T_ph_orb / DAY

    return dv1, dv2, dv_ph, t_tr, t_ph


def compute_tug_spirals(rh_raan_deg=RH_RAAN):
    """
    Edelbaum iterative sizing for all N_DEB debris-to-RH tug spirals.

    Parameters
    ----------
    rh_raan_deg : float
        RH RAAN in degrees.  Defaults to the value in config.py.

    Returns
    -------
    tug_dv, tug_time, tug_mprop, tug_mwet : (N_DEB,) arrays
    """
    tug_dv    = np.zeros(N_DEB)
    tug_time  = np.zeros(N_DEB)
    tug_mprop = np.zeros(N_DEB)
    tug_mwet  = np.zeros(N_DEB)

    v2 = np.sqrt(MU / RH_SMA)
    i2, o2 = RH_INC, rh_raan_deg

    for k in range(N_DEB):
        m_pl  = TUG_DRY + MASS[k]
        m_wet = m_pl * 1.35
        i1, o1 = INC[k], RAAN[k]
        v1 = np.sqrt(MU / SMA[k])
        cos_d = (np.cos(i1*D2R)*np.cos(i2*D2R) +
                 np.sin(i1*D2R)*np.sin(i2*D2R)*np.cos((o2-o1)*D2R))
        dth  = np.arccos(np.clip(cos_d, -1.0, 1.0))
        dv_e = np.sqrt(v1**2 + v2**2 - 2*v1*v2*np.cos(np.pi/2*dth))
        for _ in range(60):
            mf   = m_wet * np.exp(-dv_e / TUG_VEX)
            mnew = m_pl + (m_wet - mf)
            if abs(mnew - m_wet) < 0.05:
                break
            m_wet = 0.6*m_wet + 0.4*mnew
        t_s = (m_wet * TUG_VEX / TUG_THR) * (1 - np.exp(-dv_e / TUG_VEX))
        tug_dv[k]    = dv_e
        tug_time[k]  = t_s / DAY
        tug_mprop[k] = m_wet - m_pl
        tug_mwet[k]  = m_wet

    return tug_dv, tug_time, tug_mprop, tug_mwet


def phasing_hohmann(sma):
    """
    Double-Hohmann phasing burn and duration at a given orbit.

    Returns
    -------
    dv : float  [m/s]
    t  : float  [days]
    """
    T_tgt    = 2*np.pi * np.sqrt(sma**3 / MU)
    T_ph_orb = T_tgt * (1.0 - 1.0 / (4.0 * N_PHASE_REV))
    a_ph     = (MU * (T_ph_orb / (2*np.pi))**2) ** (1.0/3.0)
    dv       = 2*abs(np.sqrt(MU/sma) - np.sqrt(MU*(2/sma - 1/a_ph)))
    t        = N_PHASE_REV * T_ph_orb / DAY
    return float(dv), float(t)


# =============================================================================
# MODULE-LEVEL PRE-COMPUTATION  (default RH RAAN from config)
# =============================================================================

print("  Pre-computing transfer table...", end=' ', flush=True)
t0 = time.time()
DV1, DV2, DV_PH, T_TR, T_PH = build_transfer_table()
DV_LEG = DV1 + DV2 + DV_PH
T_LEG  = T_TR + T_PH
print(f"done ({time.time()-t0:.2f}s)")

print("  Pre-computing tug spirals...", end=' ', flush=True)
TUG_DV, TUG_TIME, TUG_MPROP, TUG_MWET = compute_tug_spirals()
TUG_WET_LOADED = TUG_DRY + TUG_MPROP.max()
print(f"done ({time.time()-t0:.2f}s)")
print(f"  Tug loaded wet mass (worst-case sizing): {TUG_WET_LOADED:.1f} kg")

_, T_PH_RH = phasing_hohmann(RH_SMA)
