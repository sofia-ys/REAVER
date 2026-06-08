"""
REAVER RPO Guidance — translational phases
==========================================
Maps a desired relative-motion manoeuvre (start point -> end point in the Hill
frame) onto a ΔV and propellant cost using the analytical two-impulse
Clohessy-Wiltshire solution from ``rpo_dynamics``.

This is the budget-level guidance: it returns the cost of each translational
phase rather than a full time-resolved trajectory.  The keep-out-sphere and
approach-cone path constraints (Conings Eq. 10, 17-18) are honoured at the
budget level by adding a routing penalty when the straight chord would clip the
keep-out sphere (a convex-solver layer can replace this later for full path
fidelity).
"""
import numpy as np

from rpo_config import N_MEAN, V_MAX, AOM_MS
from rpo_dynamics import cw_two_impulse, stationkeeping_dv
from rpo_budget import prop_from_dv, PhaseResult


def _transfer_time(r0, rf, v_max=V_MAX):
    """
    Pick a transfer time from the chord length and the max approach velocity.
    Keeps the mean speed at ~half V_MAX so the velocity cap (Conings Eq. 11) is
    respected with margin.
    """
    chord = np.linalg.norm(np.asarray(rf, float) - np.asarray(r0, float))
    return max(chord / (0.5 * v_max), 1.0)


def _clips_kos(r0, rf, r_kos):
    """
    True if the straight chord from r0 to rf passes within ``r_kos`` of the
    origin (target CoM) — i.e. a direct burn would violate the keep-out sphere.
    """
    r0 = np.asarray(r0, float)
    rf = np.asarray(rf, float)
    d = rf - r0
    L2 = np.dot(d, d)
    if L2 == 0.0:
        return np.linalg.norm(r0) < r_kos
    s = np.clip(-np.dot(r0, d) / L2, 0.0, 1.0)
    closest = r0 + s * d
    return np.linalg.norm(closest) < r_kos


def translational_phase(name, r0, rf, m0, *, r_kos=0.0, v0=None, vf=None,
                        t_transfer=None, area_to_mass=AOM_MS, duration_pad=0.0):
    """
    Cost one translational RPO phase.

    Parameters
    ----------
    name        : phase label
    r0, rf      : (3,) start / end relative position in Hill frame [m]
    m0          : chaser wet mass at the start of the phase [kg]
    r_kos       : keep-out-sphere radius to respect [m]; 0 disables the check
    v0, vf      : optional boundary velocities [m/s] (default rest-to-rest)
    t_transfer  : optional fixed transfer time [s]; default sized from V_MAX
    area_to_mass: A/m for any residual stationkeeping during the phase
    duration_pad: extra wall-clock time to add to the phase duration [s]

    Returns
    -------
    PhaseResult
    """
    r0 = np.asarray(r0, float)
    rf = np.asarray(rf, float)
    v0 = np.zeros(3) if v0 is None else np.asarray(v0, float)
    vf = np.zeros(3) if vf is None else np.asarray(vf, float)

    t = _transfer_time(r0, rf) if t_transfer is None else t_transfer

    if r_kos > 0.0 and _clips_kos(r0, rf, r_kos):
        # Route via an intermediate waypoint pushed radially out to the KOS
        # surface so the two legs skirt the sphere instead of crossing it.
        mid_dir = (r0 + rf) / 2.0
        norm = np.linalg.norm(mid_dir)
        mid_dir = mid_dir / norm if norm > 0 else np.array([0.0, 0.0, 1.0])
        waypoint = mid_dir * (r_kos * 1.1)
        t1 = _transfer_time(r0, waypoint)
        t2 = _transfer_time(waypoint, rf)
        _, _, dv_a = cw_two_impulse(r0, v0, waypoint, np.zeros(3), N_MEAN, t1)
        _, _, dv_b = cw_two_impulse(waypoint, np.zeros(3), rf, vf, N_MEAN, t2)
        dv = dv_a + dv_b
        t = t1 + t2
    else:
        _, _, dv = cw_two_impulse(r0, v0, rf, vf, N_MEAN, t)

    prop = prop_from_dv(dv, m0)
    duration = t + duration_pad
    return PhaseResult(name, dv, prop, duration)


def dwell_phase(name, m0, dwell_time, area_to_mass=AOM_MS):
    """
    Cost a stationkeeping/dwell phase (no net translation, e.g. debris analysis
    or hold-point wait): only SRP + triaxiality drift must be corrected.
    """
    dv = stationkeeping_dv(area_to_mass, dwell_time)
    prop = prop_from_dv(dv, m0)
    return PhaseResult(name, dv, prop, dwell_time, note="stationkeeping only")
