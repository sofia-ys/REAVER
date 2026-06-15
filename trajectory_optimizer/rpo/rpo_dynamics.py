"""
REAVER RPO Dynamics
===================
Physics layer for the close-range RPO simulation.

Implements the three equation sets from Conings & Mooij (AIAA-2025-0085),
adapted from LEO to the REAVER GEO recycling-hub orbit:

  * Translational : Clohessy-Wiltshire / Hill relative motion (paper Eq. 4-6),
                    here as the closed-form state-transition matrix (STM) so
                    that two-impulse transfers can be solved analytically.
  * Rotational    : Euler rigid-body equation (paper Eq. 2) + quaternion
                    kinematics (paper Eq. 3) for tumbling-target propagation.
  * Perturbations : SRP + Earth triaxiality (GEO-relevant; drag & magnetic
                    torque from the paper are dropped — negligible at GEO).

All relative-motion quantities are expressed in the target Hill frame
(x = radial / R-bar, y = along-track / V-bar, z = cross-track / H-bar).
"""
import numpy as np
from scipy.integrate import solve_ivp

from rpo_config import P_SOL, C_R, A_TRIAX


# =============================================================================
# CLOHESSY-WILTSHIRE RELATIVE TRANSLATIONAL MOTION
# =============================================================================

def cw_stm(n, t):
    """
    Clohessy-Wiltshire state-transition matrix Phi(t) for a circular reference
    orbit of mean motion ``n``.

    State ordering: x = [r(3), v(3)] in the Hill frame.
    Returns the 6x6 matrix mapping x(0) -> x(t) under force-free relative motion.
    """
    s = np.sin(n * t)
    c = np.cos(n * t)

    Phi_rr = np.array([
        [4.0 - 3.0 * c, 0.0, 0.0],
        [6.0 * (s - n * t), 1.0, 0.0],
        [0.0, 0.0, c],
    ])
    Phi_rv = np.array([
        [s / n, 2.0 * (1.0 - c) / n, 0.0],
        [2.0 * (c - 1.0) / n, (4.0 * s - 3.0 * n * t) / n, 0.0],
        [0.0, 0.0, s / n],
    ])
    Phi_vr = np.array([
        [3.0 * n * s, 0.0, 0.0],
        [6.0 * n * (c - 1.0), 0.0, 0.0],
        [0.0, 0.0, -n * s],
    ])
    Phi_vv = np.array([
        [c, 2.0 * s, 0.0],
        [-2.0 * s, 4.0 * c - 3.0, 0.0],
        [0.0, 0.0, c],
    ])

    Phi = np.zeros((6, 6))
    Phi[:3, :3] = Phi_rr
    Phi[:3, 3:] = Phi_rv
    Phi[3:, :3] = Phi_vr
    Phi[3:, 3:] = Phi_vv
    return Phi


def cw_two_impulse(r0, v0, rf, vf, n, t):
    """
    Two-impulse Clohessy-Wiltshire targeting.

    Find the burn at t=0 that places the chaser at relative position ``rf``
    after time ``t``, then the burn at t=t that matches the desired final
    velocity ``vf``.

    Parameters
    ----------
    r0, v0 : (3,) initial relative position [m] and velocity [m/s]
    rf, vf : (3,) desired final relative position [m] and velocity [m/s]
    n      : mean motion [rad/s]
    t      : transfer time [s]

    Returns
    -------
    dv1, dv2 : (3,) the two impulse vectors [m/s]
    dv_total : float, |dv1| + |dv2|  [m/s]
    """
    r0 = np.asarray(r0, float)
    v0 = np.asarray(v0, float)
    rf = np.asarray(rf, float)
    vf = np.asarray(vf, float)

    Phi = cw_stm(n, t)
    Phi_rr, Phi_rv = Phi[:3, :3], Phi[:3, 3:]
    Phi_vr, Phi_vv = Phi[3:, :3], Phi[3:, 3:]

    # Velocity just after the first burn that hits rf at time t.
    v0_plus = np.linalg.solve(Phi_rv, rf - Phi_rr @ r0)
    dv1 = v0_plus - v0

    # Velocity just before the second burn.
    vf_minus = Phi_vr @ r0 + Phi_vv @ v0_plus
    dv2 = vf - vf_minus

    dv_total = float(np.linalg.norm(dv1) + np.linalg.norm(dv2))
    return dv1, dv2, dv_total


def cw_sample_path(r0, rf, n, t, n_pts=80, v0=None, vf=None):
    """
    Sample the relative-position trajectory of a two-impulse CW transfer for
    plotting.  Applies the first impulse at t=0, then propagates the resulting
    state with the STM over [0, t].

    Returns (n_pts, 3) array of Hill-frame positions [m].
    """
    v0 = np.zeros(3) if v0 is None else np.asarray(v0, float)
    vf = np.zeros(3) if vf is None else np.asarray(vf, float)
    dv1, _, _ = cw_two_impulse(r0, v0, rf, vf, n, t)
    x0 = np.concatenate([np.asarray(r0, float), v0 + dv1])
    ts = np.linspace(0.0, t, n_pts)
    return np.array([(cw_stm(n, tt) @ x0)[:3] for tt in ts])


def cw_min_dv_transfer(r0, rf, n, t, v0=None, vf=None):
    """
    Convenience wrapper for a rest-to-rest CW transfer (v0 = vf = 0 by default),
    returning only the total ΔV.  This is the dominant cost model for each
    translational RPO phase used in the budget.
    """
    v0 = np.zeros(3) if v0 is None else v0
    vf = np.zeros(3) if vf is None else vf
    _, _, dv = cw_two_impulse(r0, v0, rf, vf, n, t)
    return dv


# =============================================================================
# GEO PERTURBATIONS  (stationkeeping driver during dwell / hold phases)
# =============================================================================

def srp_accel(area_to_mass, c_r=C_R):
    """Solar-radiation-pressure acceleration magnitude [m/s^2]."""
    return c_r * area_to_mass * P_SOL


def stationkeeping_dv(area_to_mass, dwell_time, a_triax=A_TRIAX, c_r=C_R):
    """
    ΔV required to null the drift accumulated during a hold/dwell of
    ``dwell_time`` seconds against the combined SRP + triaxiality acceleration.

    First-order model: a constant disturbance acceleration integrates into a
    drift velocity  a*T  that must be cancelled  ->  ΔV ≈ |a_srp + a_triax|*T.
    """
    a_dist = srp_accel(area_to_mass, c_r) + a_triax
    return float(a_dist * dwell_time)


# =============================================================================
# RIGID-BODY ATTITUDE  (Euler rigid body + quaternion kinematics)
# =============================================================================

def box_inertia(mass, dims):
    """
    Principal inertia tensor of a uniform rectangular box.

    Parameters
    ----------
    mass : float [kg]
    dims : (a, b, c) edge lengths along x, y, z [m]

    Returns
    -------
    (3,3) diagonal inertia tensor [kg m^2]
    """
    a, b, c = dims
    Ixx = mass / 12.0 * (b**2 + c**2)
    Iyy = mass / 12.0 * (a**2 + c**2)
    Izz = mass / 12.0 * (a**2 + b**2)
    return np.diag([Ixx, Iyy, Izz])


def debris_inertia(mass, bus_dims, panel_span, panel_frac):
    """
    Asymmetric inertia tensor of a GEO comsat: central bus plus two deployed
    solar panels modelled as point masses at the panel half-span (the panels
    dominate the off-axis moments, per the checklist).

    Panels are assumed to deploy along the body y-axis; their large lever arm
    inflates Ixx and Izz while leaving Iyy near the bus value.
    """
    m_bus = mass * (1.0 - panel_frac)
    m_pan = mass * panel_frac                      # both panels combined
    J = box_inertia(m_bus, bus_dims)

    half = panel_span / 2.0
    # two panels, each m_pan/2 at +/- half along y: contributes to Ixx and Izz
    J[0, 0] += m_pan * half**2                     # about x: y-offset matters
    J[2, 2] += m_pan * half**2                     # about z: y-offset matters
    return J


def parallel_axis(J_cm, mass, offset):
    """
    Shift an inertia tensor from a body's own CoM to a reference point using the
    parallel-axis theorem.  ``offset`` is the vector from the reference point to
    the body CoM [m].
    """
    d = np.asarray(offset, float)
    return J_cm + mass * (np.dot(d, d) * np.eye(3) - np.outer(d, d))


def combined_inertia(bodies):
    """
    Inertia tensor of an assembly about its common centre of mass.

    Parameters
    ----------
    bodies : list of (mass, J_cm, position) tuples, where ``position`` is the
             body CoM location in a shared reference frame [m].

    Returns
    -------
    J_total : (3,3) inertia tensor about the assembly CoM [kg m^2]
    r_cm    : (3,) assembly centre of mass [m]
    m_total : float total mass [kg]
    """
    m_total = sum(b[0] for b in bodies)
    r_cm = sum(b[0] * np.asarray(b[2], float) for b in bodies) / m_total

    J_total = np.zeros((3, 3))
    for m, J_cm, pos in bodies:
        offset = np.asarray(pos, float) - r_cm
        J_total += parallel_axis(J_cm, m, offset)
    return J_total, r_cm, m_total


def quat_kinematics(q, omega):
    """
    Quaternion derivative qdot = 1/2 * Xi(q) * omega   (Conings Eq. 3).
    Quaternion convention: q = [q1, q2, q3, q4] with q4 the scalar part.
    """
    q1, q2, q3, q4 = q
    Xi = np.array([
        [q4, -q3, q2],
        [q3, q4, -q1],
        [-q2, q1, q4],
        [-q1, -q2, -q3],
    ])
    return 0.5 * Xi @ omega


def euler_rigid_body(omega, J, M_net=None):
    """
    Euler's rigid-body equation  omega_dot = J^-1 (M_net - omega x J omega)
    (Conings Eq. 2).  ``M_net`` defaults to zero (free tumble).
    """
    M_net = np.zeros(3) if M_net is None else np.asarray(M_net, float)
    Jw = J @ omega
    return np.linalg.solve(J, M_net - np.cross(omega, Jw))


def propagate_tumble(q0, omega0, J, t_span, n_eval=200):
    """
    Propagate a free-tumbling rigid body (Euler + quaternion) over ``t_span``.

    Returns
    -------
    t       : (n_eval,) time grid [s]
    q_hist  : (n_eval, 4) quaternion history
    w_hist  : (n_eval, 3) body-rate history [rad/s]
    """
    def rhs(_t, y):
        q = y[:4]
        w = y[4:]
        qd = quat_kinematics(q, w)
        wd = euler_rigid_body(w, J)
        return np.concatenate([qd, wd])

    y0 = np.concatenate([q0, omega0])
    t_eval = np.linspace(0.0, t_span, n_eval)
    sol = solve_ivp(rhs, (0.0, t_span), y0, t_eval=t_eval,
                    rtol=1e-9, atol=1e-12, method='RK45')
    y = sol.y.T
    q_hist = y[:, :4]
    q_hist /= np.linalg.norm(q_hist, axis=1, keepdims=True)   # renormalise
    return sol.t, q_hist, y[:, 4:]
