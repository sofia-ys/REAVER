"""
REAVER RPO -- Plots
===================
Four separate report-ready figures for the close-range RPO simulation,
using standard matplotlib style (no custom colours).

  Figure 1 -- rpo_fig1_debris_traj.png
      Debris-RPO trajectory in the target Hill frame, with keep-out
      spheres and the final approach cone.

  Figure 2 -- rpo_fig2_rh_traj.png
      RH docking-cycle trajectory (cooperative, fixed docking axis).

  Figure 3 -- rpo_fig3_dv_breakdown.png
      Phase-by-phase dV breakdown, debris RPO vs one RH cycle.

  Figure 4 -- rpo_fig4_tumble_sweep.png
      Tumble-rate sensitivity of DV_RPO_DEBRIS and DV_RPO_DETUMBLE.

    cd c:\\Projects\\DSE\\REAVER
    python trajectory_optimizer/rpo_plots.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

from rpo_config import (N_MEAN, INSPECT_STANDOFF, H_BAR_STANDOFF, R_KOS2, R_KOS2_MEET, R_CUTOFF,
                        R_KOS1_DOCK, R_KOS2_DOCK, DOCK_PORT_HEIGHT,
                        APPROACH_CONE_DEG, V_MAX, OMEGA_REQ_DPS,
                        DEBRIS_PANEL_SPAN_WC, KOS1_MARGIN)
R_KOS1 = (DEBRIS_PANEL_SPAN_WC / 2.0) * KOS1_MARGIN   # half-span + 20% safety margin
from rpo_dynamics import cw_sample_path
from rpo_debris_sim import simulate_debris_rpo
from rpo_rh_sim import (simulate_rh_cycle, R_MEET_START,
                        R_CAPTURE, R_DOCK_CLOSE, R_PORT)

SAVE_DIR = r'C:\Projects\DSE\REAVER\trajectory_optimizer'


# =============================================================================
# SHARED HELPERS
# =============================================================================

def _transfer_time(r0, rf):
    chord = np.linalg.norm(np.asarray(rf, float) - np.asarray(r0, float))
    return max(chord / (0.5 * V_MAX), 1.0)


def _path_through(waypoints):
    """CW-sampled path through a list of Hill-frame waypoints."""
    segs = []
    for a, b in zip(waypoints[:-1], waypoints[1:]):
        a = np.asarray(a, float)
        b = np.asarray(b, float)
        segs.append(cw_sample_path(a, b, N_MEAN, _transfer_time(a, b), n_pts=80))
    return np.vstack(segs)



def _add_shadows(ax, path, color='black', alpha=0.40):
    """Project path onto the three bounding planes (floor, back wall, side wall)."""
    X, Y, Z = path[:, 0], path[:, 1], path[:, 2]
    xl, yl, zl = ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()
    kw = dict(color=color, lw=1.8, alpha=alpha, zorder=1)
    ax.plot(X, Y, np.full_like(Z, zl[0]),   **kw)  # floor  (H-bar = min)
    ax.plot(X, np.full_like(Y, yl[1]), Z,   **kw)  # back wall (V-bar = max)
    ax.plot(np.full_like(X, xl[0]), Y, Z,   **kw)  # side wall (R-bar = min)


def _wire_sphere(ax, r, color, alpha=0.50, n_lat=10, n_lon=16, full=False,
                 ypos=False, pole_axis='z'):
    if ypos:
        u = np.linspace(0, np.pi, n_lon)   # sin(u) >= 0  →  y >= 0 hemisphere
        v = np.linspace(0, np.pi, n_lat)
    else:
        u = np.linspace(0, 2 * np.pi, n_lon)
        v = np.linspace(0, np.pi if full else np.pi / 2, n_lat)
    s = r * np.outer(np.cos(u), np.sin(v))   # "x-like"
    t = r * np.outer(np.sin(u), np.sin(v))   # "y-like"
    p = r * np.outer(np.ones_like(u), np.cos(v))  # pole axis
    if pole_axis == 'y':
        x, y, z = s, p, t   # poles along V-bar (y)
    else:
        x, y, z = s, t, p   # poles along H-bar (z), default
    ax.plot_wireframe(x, y, z, color=color, lw=0.8, alpha=alpha)


def _approach_cone(ax, axis_dir, half_deg, length, color, alpha=0.22):
    axis_dir = np.asarray(axis_dir, float)
    axis_dir /= np.linalg.norm(axis_dir)
    tmp = np.array([1., 0., 0.]) if abs(axis_dir[0]) < 0.9 else np.array([0., 1., 0.])
    e1 = np.cross(axis_dir, tmp); e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis_dir, e1)
    a = np.radians(half_deg)

    h   = np.linspace(0, length, 35)
    phi = np.linspace(0, 2 * np.pi, 60)
    H, Phi = np.meshgrid(h, phi)
    r = H * np.tan(a)
    X = H * axis_dir[0] + r * (np.cos(Phi) * e1[0] + np.sin(Phi) * e2[0])
    Y = H * axis_dir[1] + r * (np.cos(Phi) * e1[1] + np.sin(Phi) * e2[1])
    Z = H * axis_dir[2] + r * (np.cos(Phi) * e1[2] + np.sin(Phi) * e2[2])
    ax.plot_surface(X, Y, Z, color=color, alpha=alpha, linewidth=0, antialiased=True)


# =============================================================================
# FIGURE 1 -- DEBRIS RPO TRAJECTORY
# =============================================================================

def fig_debris_traj(save=True):
    """
    3-D Hill-frame trajectory for the debris RPO.  Two-stage approach:
    (1) hold + inspect on V-bar (the only naturally fixed CW hold point) to
    identify the spin rate + tumble axis, then (2) transition to the tumble
    axis (modelled as +H-bar) and run the capture corridor KOS1 -> KOS2.
    """
    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    C_INSP   = prop_cycle[0]    # blue   — V-bar inspection approach
    C_TRANS  = prop_cycle[2]    # green  — transition to tumble axis
    C_CAP    = prop_cycle[1]    # orange — tumble-axis capture corridor
    C_KOS1   = 'darkorange'
    C_KOS2   = 'teal'
    C_CONE   = 'gold'
    C_TARGET = prop_cycle[4]    # purple

    # Stage 1: arrive along -V-bar to the inspection hold.
    insp = _path_through([[0.0, -(INSPECT_STANDOFF + 8.0), 0.0],
                          [0.0, -INSPECT_STANDOFF, 0.0]])
    # Stage 1->2: transition arc from the V-bar hold to the tumble-axis standoff.
    trans = _path_through([[0.0, -INSPECT_STANDOFF, 0.0],
                           [0.0, 0.0, H_BAR_STANDOFF]])
    # Stage 2: capture corridor down the tumble axis.
    cap = _path_through([[0.0, 0.0, H_BAR_STANDOFF],
                         [0.0, 0.0, R_KOS1],
                         [0.0, 0.0, R_KOS2]])

    fig = plt.figure(figsize=(7.5, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-12, 12)
    ax.set_ylim(-30, 10)
    ax.set_zlim(0, 23)
    # Equal data-unit scaling so the keep-out spheres render round (not tilted).
    ax.set_box_aspect((24, 40, 23))

    _add_shadows(ax, np.vstack([insp, trans, cap]))

    _wire_sphere(ax, R_KOS1, C_KOS1)
    _wire_sphere(ax, R_KOS2, C_KOS2)
    _approach_cone(ax, [0, 0, 1], APPROACH_CONE_DEG, R_KOS1, C_CONE, alpha=0.45)

    ax.plot(insp[:, 0], insp[:, 1], insp[:, 2], color=C_INSP, lw=2.0, zorder=5,
            label='P1 approach (along -V-bar)')
    ax.plot(trans[:, 0], trans[:, 1], trans[:, 2], color=C_TRANS, lw=2.0, ls='--',
            zorder=5, label='P3 transition (V-bar → H-bar)')
    ax.plot(cap[:, 0], cap[:, 1], cap[:, 2], color=C_CAP, lw=2.0, zorder=5,
            label='P6 capture corridor (H-bar)')

    _halo = [pe.withStroke(linewidth=3, foreground='white')]
    ax.scatter(0.0, -INSPECT_STANDOFF, 0.0, color='k', s=55, marker='o',
               depthshade=False, zorder=6, label='V-bar inspection hold (20 m)',
               path_effects=_halo)
    ax.scatter(0.0, 0.0, H_BAR_STANDOFF, color='k', s=45, marker='^',
               depthshade=False, zorder=6, label='H-bar standoff (20 m)',
               path_effects=_halo)
    ax.scatter(0.0, 0.0, R_KOS1, color='k', s=50, marker='D', depthshade=False,
               zorder=6, label='KOS1 hold point', path_effects=_halo)
    ax.scatter(0.0, 0.0, R_KOS2, color='k', s=50, marker='s', depthshade=False,
               zorder=6, label='KOS2 (MEV lock)', path_effects=_halo)
    ax.scatter(0, 0, 0, color=C_TARGET, s=80, marker='*', depthshade=False,
               zorder=6, label='Target CoM', path_effects=_halo)

    ax.set_xlabel('R-bar [m]')
    ax.set_ylabel('V-bar [m]')
    ax.set_zlabel('H-bar [m]', labelpad=10)

    _kos_handles = [
        Line2D([0], [0], color=C_KOS1, lw=1.5, label=f'KOS1 — interaction ({R_KOS1:.1f} m)'),
        Line2D([0], [0], color=C_KOS2, lw=1.5, label=f'KOS2 — MEV lock ({R_KOS2:.2f} m)'),
    ]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + _kos_handles, labels=labels + [h.get_label() for h in _kos_handles],
              fontsize=7.5, loc='upper left')
    fig.tight_layout()

    path_out = f'{SAVE_DIR}\\rpo_fig1_debris_traj.png'
    if save:
        fig.savefig(path_out, dpi=150, bbox_inches='tight', pad_inches=0.3)
        print(f'Saved -> {path_out}')
    return fig


# =============================================================================
# FIGURE 2a -- RH PHASE A: MEET (close approach to tug+debris, V-bar axis)
# =============================================================================

def fig_rh_meet(save=True):
    """
    Close-range cooperative approach of the MS to the tug+debris assembly.

    Differences from the debris approach (Fig 1):
      - V-bar approach axis (direct, no H-bar detour — target is cooperative)
      - No 8 h dwell: tug relays state telemetry, LiDAR estimates relative pose,
        go/no-go from ground takes ~30 min
      - KOS1 still sized on debris panel span (panels present, not spinning)
      - Approach cone same 5 deg terminal corridor

    Frame origin = tug+debris CoM.
    MS approaches from +V-bar (from the RH side) toward the tug at the origin.
    """
    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    C_PATH   = prop_cycle[0]   # blue
    C_KOS1   = 'darkorange'
    C_KOS2   = 'teal'
    C_CONE   = 'gold'
    C_TARGET = prop_cycle[4]   # purple — tug+debris CoM

    # MS approaches from -V-bar; tug+debris CoM at origin.
    wp = [
        [0.0, -30.0,        0.0],   # close-approach start
        [0.0, -R_KOS1,      0.0],   # KOS1 hold (synch + comm relay, go/no-go)
        [0.0, -R_KOS2_MEET, 0.0],   # KOS2 — bare-arm capture range (arm length only)
    ]
    path = _path_through(wp)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlim(-12, 12)
    ax.set_ylim(-35, 15)
    ax.set_zlim(-12, 12)

    _add_shadows(ax, path)

    _wire_sphere(ax, R_KOS1, C_KOS1, full=True, pole_axis='y')
    _wire_sphere(ax, R_KOS2_MEET, C_KOS2, full=True, pole_axis='y')
    # Cone opens toward -V-bar (toward approaching MS)
    _approach_cone(ax, [0, -1, 0], APPROACH_CONE_DEG, R_KOS1, C_CONE, alpha=0.45)

    ax.plot(path[:, 0], path[:, 1], path[:, 2],
            color=C_PATH, lw=2.0, label='MS approach path (V-bar, cooperative)', zorder=5)
    _halo = [pe.withStroke(linewidth=3, foreground='white')]
    ax.scatter(*wp[0], color='k', s=50, marker='o', depthshade=False,
               zorder=6, label='30 m standoff', path_effects=_halo)
    ax.scatter(*wp[1], color='k', s=50, marker='D', depthshade=False,
               zorder=6, label='KOS1 hold point', path_effects=_halo)
    ax.scatter(*wp[2], color='k', s=50, marker='s', depthshade=False,
               zorder=6, label='KOS2 hold point', path_effects=_halo)
    ax.scatter(0, 0, 0, color=C_TARGET, s=80, marker='^',
               depthshade=False, zorder=6, label='Tug+debris CoM',
               path_effects=_halo)

    ax.set_xlabel('R-bar [m]')
    ax.set_ylabel('V-bar [m]')
    ax.set_zlabel('H-bar [m]', labelpad=10)

    _kos_handles = [
        Line2D([0], [0], color=C_KOS1, lw=1.5,
               label=f'KOS1 — panel clearance ({R_KOS1:.1f} m)'),
        Line2D([0], [0], color=C_KOS2, lw=1.5,
               label=f'KOS2 — bare-arm capture ({R_KOS2_MEET:.2f} m)'),
    ]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + _kos_handles, labels=labels + [h.get_label() for h in _kos_handles],
              fontsize=8, loc='upper left')
    fig.tight_layout()

    path_out = f'{SAVE_DIR}\\rpo_fig2a_rh_meet.png'
    if save:
        fig.savefig(path_out, dpi=150, bbox_inches='tight', pad_inches=0.5)
        print(f'Saved -> {path_out}')
    return fig


# =============================================================================
# FIGURE 2b -- RH PHASE B: DOCK (combined body 500 m -> RH port)
# =============================================================================

def fig_rh_dock(save=True):
    """
    Phase B — DOCK.  The combined body (MS + EXTENDED arm + tug + debris) runs the
    close-approach corridor to the RH and hard-docks via a dedicated port on TOP
    of the MS bus — the arm stays extended, holding the payload out to the side.
    Same close-approach treatment as Fig 1 / Fig 2a, in the RH-port frame:
    30 m standoff -> KOS1_DOCK alignment hold -> KOS2_DOCK top-port contact.
    """
    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    C_PATH = prop_cycle[1]    # orange — combined body
    C_KOS1 = 'darkorange'
    C_KOS2 = 'teal'
    C_CONE = 'gold'
    C_PORT = prop_cycle[4]    # purple — RH docking port

    # MS approaches from -V-bar; RH docking port at origin.
    wp = [
        [0.0, -30.0,        0.0],   # close-approach start (end of the 500 m carry)
        [0.0, -R_KOS1_DOCK, 0.0],   # KOS1 alignment hold (go/no-go)
        [0.0, -R_KOS2_DOCK, 0.0],   # KOS2 top-port hard-dock contact
    ]
    path = _path_through(wp)

    fig = plt.figure(figsize=(8, 6))
    ax  = fig.add_subplot(111, projection='3d')

    ax.set_xlim(-14, 14)
    ax.set_ylim(-35, 15)
    ax.set_zlim(-14, 14)

    _add_shadows(ax, path)

    _wire_sphere(ax, R_KOS1_DOCK, C_KOS1, full=True, pole_axis='y')
    _wire_sphere(ax, R_KOS2_DOCK, C_KOS2, full=True, pole_axis='y')
    # Cone opens toward -V-bar (toward the approaching MS)
    _approach_cone(ax, [0, -1, 0], APPROACH_CONE_DEG, R_KOS1_DOCK, C_CONE, alpha=0.45)

    ax.plot(path[:, 0], path[:, 1], path[:, 2],
            color=C_PATH, lw=2.0, label='Combined-body approach (V-bar)', zorder=5)
    _halo = [pe.withStroke(linewidth=3, foreground='white')]
    ax.scatter(*wp[0], color='k', s=50, marker='o', depthshade=False,
               zorder=6, label='30 m standoff (from 500 m carry)', path_effects=_halo)
    ax.scatter(*wp[1], color='k', s=50, marker='D', depthshade=False,
               zorder=6, label='KOS1 alignment hold', path_effects=_halo)
    ax.scatter(*wp[2], color='k', s=50, marker='s', depthshade=False,
               zorder=6, label='KOS2 top-port contact', path_effects=_halo)
    ax.scatter(0, 0, 0, color=C_PORT, s=120, marker='s',
               depthshade=False, zorder=6, label='RH docking port', path_effects=_halo)

    ax.set_xlabel('R-bar [m]')
    ax.set_ylabel('V-bar [m]')
    ax.set_zlabel('H-bar [m]', labelpad=10)

    _kos_handles = [
        Line2D([0], [0], color=C_KOS1, lw=1.5,
               label=f'KOS1 — alignment hold ({R_KOS1_DOCK:.1f} m)'),
        Line2D([0], [0], color=C_KOS2, lw=1.5,
               label=f'KOS2 — top-port dock ({R_KOS2_DOCK:.2f} m)'),
    ]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + _kos_handles, labels=labels + [h.get_label() for h in _kos_handles],
              fontsize=8, loc='upper left')
    fig.tight_layout()

    path_out = f'{SAVE_DIR}\\rpo_fig2b_rh_dock.png'
    if save:
        fig.savefig(path_out, dpi=150, bbox_inches='tight', pad_inches=0.5)
        print(f'Saved -> {path_out}')
    return fig


# =============================================================================
# FIGURE 2c -- RH PHASE B: CARRY (500 m -> 30 m far-range approach)
# =============================================================================

def fig_rh_carry(save=True):
    """
    Phase B, far-range context: the combined body (MS + extended arm + tug +
    debris) carries from the 500 m capture hold in to the 30 m standoff where the
    close-approach corridor (Fig 2b) begins.  RH-port frame; the CW two-impulse
    transfer shows the characteristic R-bar excursion of the coasting arc.
    """
    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    C_CARRY = prop_cycle[1]   # orange — combined body carry
    C_START = prop_cycle[0]   # blue   — capture / start
    C_HAND  = prop_cycle[2]   # green  — handoff to close approach
    C_PORT  = prop_cycle[4]   # purple — RH docking port

    path = cw_sample_path(R_CAPTURE, R_DOCK_CLOSE, N_MEAN,
                          _transfer_time(R_CAPTURE, R_DOCK_CLOSE), n_pts=160)

    fig = plt.figure(figsize=(8, 6))
    ax  = fig.add_subplot(111, projection='3d')

    # Auto-fit the R-bar excursion so the arc is visible.
    x_pad = max(2.0, 1.2 * np.max(np.abs(path[:, 0])))
    ax.set_xlim(-x_pad, x_pad)
    ax.set_ylim(-520, 20)
    ax.set_zlim(-1.0, 1.0)

    _add_shadows(ax, path)

    ax.plot(path[:, 0], path[:, 1], path[:, 2],
            color=C_CARRY, lw=2.0, zorder=5,
            label='Combined-body carry (B1 — 500 m → 30 m)')

    _halo = [pe.withStroke(linewidth=3, foreground='white')]
    ax.scatter(*R_CAPTURE, color=C_START, s=90, marker='^', depthshade=False,
               zorder=6, label='Arm capture — combined body start (500 m)', path_effects=_halo)
    ax.scatter(*R_DOCK_CLOSE, color=C_HAND, s=70, marker='o', depthshade=False,
               zorder=6, label='30 m standoff — handoff to close approach (Fig 2b)',
               path_effects=_halo)
    ax.scatter(0, 0, 0, color=C_PORT, s=120, marker='s', depthshade=False,
               zorder=6, label='RH docking port', path_effects=_halo)

    ax.set_xlabel('R-bar [m]')
    ax.set_ylabel('V-bar [m]')
    ax.set_zlabel('H-bar [m]', labelpad=10)
    ax.legend(fontsize=8, loc='upper left')
    fig.tight_layout()

    path_out = f'{SAVE_DIR}\\rpo_fig2c_rh_carry.png'
    if save:
        fig.savefig(path_out, dpi=150, bbox_inches='tight', pad_inches=0.3)
        print(f'Saved -> {path_out}')
    return fig


# =============================================================================
# FIGURE COMBINED -- FIG1 + FIG2a + FIG2b in one PNG
# =============================================================================

def fig_combined(save=True):
    """Fig1, Fig2a and Fig2b side-by-side in a single 1×3 figure."""
    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    _halo = [pe.withStroke(linewidth=3, foreground='white')]

    fig = plt.figure(figsize=(21, 7))
    ax1 = fig.add_subplot(131, projection='3d')
    ax2 = fig.add_subplot(132, projection='3d')
    ax3 = fig.add_subplot(133, projection='3d')

    # ── Panel (a): Debris RPO (V-bar inspect -> tumble-axis capture) ─────────
    C_INSP = prop_cycle[0]; C_TRANS = prop_cycle[2]; C_CAP = prop_cycle[1]
    C_KOS1 = 'darkorange'; C_KOS2 = 'teal'; C_CONE = 'gold'; C_TARGET = prop_cycle[4]

    insp1  = _path_through([[0.0, -(INSPECT_STANDOFF + 8.0), 0.0], [0.0, -INSPECT_STANDOFF, 0.0]])
    trans1 = _path_through([[0.0, -INSPECT_STANDOFF, 0.0], [0.0, 0.0, H_BAR_STANDOFF]])
    cap1   = _path_through([[0.0, 0.0, H_BAR_STANDOFF], [0.0, 0.0, R_KOS1], [0.0, 0.0, R_KOS2]])

    ax1.set_xlim(-12, 12); ax1.set_ylim(-26, 12); ax1.set_zlim(0, 23)
    _add_shadows(ax1, np.vstack([insp1, trans1, cap1]))
    _wire_sphere(ax1, R_KOS1, C_KOS1)
    _wire_sphere(ax1, R_KOS2, C_KOS2)
    _approach_cone(ax1, [0, 0, 1], APPROACH_CONE_DEG, R_KOS1, C_CONE, alpha=0.45)
    ax1.plot(insp1[:, 0], insp1[:, 1], insp1[:, 2], color=C_INSP, lw=2.0, label='P1 approach (-V-bar)')
    ax1.plot(trans1[:, 0], trans1[:, 1], trans1[:, 2], color=C_TRANS, lw=2.0, ls='--', label='P3 transition')
    ax1.plot(cap1[:, 0], cap1[:, 1], cap1[:, 2], color=C_CAP, lw=2.0, label='P6 capture (H-bar)')
    ax1.scatter(0.0, -INSPECT_STANDOFF, 0.0, color='k', s=40, marker='o', depthshade=False, zorder=6,
                label='V-bar inspection hold', path_effects=_halo)
    ax1.scatter(0.0, 0.0, H_BAR_STANDOFF, color='k', s=35, marker='^', depthshade=False, zorder=6,
                label='H-bar standoff', path_effects=_halo)
    ax1.scatter(0.0, 0.0, R_KOS1, color='k', s=40, marker='D', depthshade=False, zorder=6,
                label='KOS1 hold point', path_effects=_halo)
    ax1.scatter(0.0, 0.0, R_KOS2, color='k', s=40, marker='s', depthshade=False, zorder=6,
                label='KOS2 (MEV lock)', path_effects=_halo)
    ax1.scatter(0, 0, 0, color=C_TARGET, s=70, marker='*', depthshade=False,
                zorder=6, label='Target CoM', path_effects=_halo)
    ax1.set_xlabel('R-bar [m]', fontsize=8); ax1.set_ylabel('V-bar [m]', fontsize=8)
    ax1.set_zlabel('H-bar [m]', labelpad=10, fontsize=8)
    kos1_h = [Line2D([0],[0], color=C_KOS1, lw=1.2,
                     label=f'KOS1 ({R_KOS1:.1f} m)'),
              Line2D([0],[0], color=C_KOS2, lw=1.2,
                     label=f'KOS2 ({R_KOS2:.2f} m)')]
    h, l = ax1.get_legend_handles_labels()
    ax1.legend(handles=h+kos1_h, labels=l+[x.get_label() for x in kos1_h],
               fontsize=6, loc='upper left')
    ax1.set_title('(a) Debris RPO — V-bar inspect → tumble-axis capture', fontsize=9, pad=6)

    # ── Panel (b): RH Phase A (cooperative close approach) ──────────────────
    wp2 = [[0.0, -30.0, 0.0], [0.0, -R_KOS1, 0.0], [0.0, -R_KOS2_MEET, 0.0]]
    path2 = _path_through(wp2)

    ax2.set_xlim(-12, 12); ax2.set_ylim(-35, 15); ax2.set_zlim(-12, 12)
    _add_shadows(ax2, path2)
    _wire_sphere(ax2, R_KOS1, C_KOS1, full=True, pole_axis='y')
    _wire_sphere(ax2, R_KOS2_MEET, C_KOS2, full=True, pole_axis='y')
    _approach_cone(ax2, [0, -1, 0], APPROACH_CONE_DEG, R_KOS1, C_CONE, alpha=0.45)
    ax2.plot(path2[:, 0], path2[:, 1], path2[:, 2], color=C_INSP, lw=2.0,
             label='MS approach (V-bar, cooperative)')
    ax2.scatter(*wp2[0], color='k', s=40, marker='o', depthshade=False, zorder=6,
                label='30 m standoff', path_effects=_halo)
    ax2.scatter(*wp2[1], color='k', s=40, marker='D', depthshade=False, zorder=6,
                label='KOS1 hold point', path_effects=_halo)
    ax2.scatter(*wp2[2], color='k', s=40, marker='s', depthshade=False, zorder=6,
                label='KOS2 hold point', path_effects=_halo)
    ax2.scatter(0, 0, 0, color=C_TARGET, s=70, marker='^', depthshade=False,
                zorder=6, label='Tug+debris CoM', path_effects=_halo)
    ax2.set_xlabel('R-bar [m]', fontsize=8); ax2.set_ylabel('V-bar [m]', fontsize=8)
    ax2.set_zlabel('H-bar [m]', labelpad=10, fontsize=8)
    kos2_h = [Line2D([0],[0], color=C_KOS1, lw=1.2,
                     label=f'KOS1 ({R_KOS1:.1f} m)'),
              Line2D([0],[0], color=C_KOS2, lw=1.2,
                     label=f'KOS2 ({R_KOS2_MEET:.2f} m)')]
    h, l = ax2.get_legend_handles_labels()
    ax2.legend(handles=h+kos2_h, labels=l+[x.get_label() for x in kos2_h],
               fontsize=6, loc='upper left')
    ax2.set_title('(b) RH Phase A — tug meet (cooperative)', fontsize=9, pad=6)

    # ── Panel (c): RH Phase B (combined body top-port dock) ─────────────────
    C_PATH3 = prop_cycle[1]; C_PORT = prop_cycle[4]

    wp3 = [[0.0, -30.0, 0.0], [0.0, -R_KOS1_DOCK, 0.0], [0.0, -R_KOS2_DOCK, 0.0]]
    path3 = _path_through(wp3)

    ax3.set_xlim(-14, 14); ax3.set_ylim(-35, 15); ax3.set_zlim(-14, 14)
    _add_shadows(ax3, path3)
    _wire_sphere(ax3, R_KOS1_DOCK, C_KOS1, full=True, pole_axis='y')
    _wire_sphere(ax3, R_KOS2_DOCK, C_KOS2, full=True, pole_axis='y')
    _approach_cone(ax3, [0, -1, 0], APPROACH_CONE_DEG, R_KOS1_DOCK, C_CONE, alpha=0.45)
    ax3.plot(path3[:, 0], path3[:, 1], path3[:, 2], color=C_PATH3, lw=2.0,
             label='Combined-body approach (V-bar)')
    ax3.scatter(*wp3[0], color='k', s=40, marker='o', depthshade=False, zorder=6,
                label='30 m standoff', path_effects=_halo)
    ax3.scatter(*wp3[1], color='k', s=40, marker='D', depthshade=False, zorder=6,
                label='KOS1 alignment hold', path_effects=_halo)
    ax3.scatter(*wp3[2], color='k', s=40, marker='s', depthshade=False, zorder=6,
                label='KOS2 top-port contact', path_effects=_halo)
    ax3.scatter(0, 0, 0, color=C_PORT, s=90, marker='s', depthshade=False,
                zorder=6, label='RH docking port', path_effects=_halo)
    ax3.set_xlabel('R-bar [m]', fontsize=8); ax3.set_ylabel('V-bar [m]', fontsize=8)
    ax3.set_zlabel('H-bar [m]', labelpad=10, fontsize=8)
    kos3_h = [Line2D([0],[0], color=C_KOS1, lw=1.2,
                     label=f'KOS1 ({R_KOS1_DOCK:.1f} m)'),
              Line2D([0],[0], color=C_KOS2, lw=1.2,
                     label=f'KOS2 ({R_KOS2_DOCK:.2f} m)')]
    h3, l3 = ax3.get_legend_handles_labels()
    ax3.legend(handles=h3+kos3_h, labels=l3+[x.get_label() for x in kos3_h],
               fontsize=6, loc='upper left')
    ax3.set_title('(c) RH Phase B — top-port dock', fontsize=9, pad=6)

    fig.tight_layout(pad=1.5)
    path_out = f'{SAVE_DIR}\\rpo_fig_combined.png'
    if save:
        fig.savefig(path_out, dpi=150, bbox_inches='tight', pad_inches=0.3)
        print(f'Saved -> {path_out}')
    return fig


# =============================================================================
# FIGURE 3 -- PHASE dV BREAKDOWN
# =============================================================================

def fig_dv_breakdown(debris_res=None, rh_res=None, save=True):
    """Stacked bar chart: per-phase dV for debris RPO vs one RH cycle."""
    if debris_res is None:
        debris_res = simulate_debris_rpo(omega_dps=OMEGA_REQ_DPS)
    if rh_res is None:
        rh_res = simulate_rh_cycle(1884.0)

    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # ── Debris phases ──────────────────────────────────────────────────────
    d_labels = ['P1  Approach',
                'P2  Debris analysis',
                'P3  Spin sync',
                'P4  Arm extension',
                'P5a Final approach',
                'P5b Detumble',
                'P6  Retreat']
    d_vals = [p.dv for p in debris_res['phases']]
    d_cols = [prop_cycle[i % len(prop_cycle)] for i in range(len(d_labels))]

    # ── RH phases ──────────────────────────────────────────────────────────
    r_labels = ['A1  Hold',
                'A2  Approach',
                'A3  Clamp',
                'B1  Carry',
                'B1  Re-point',
                'B2  Hard-dock']
    r_vals = [p.dv for p in rh_res['phases']]
    r_cols = [prop_cycle[i % len(prop_cycle)] for i in range(len(r_labels))]

    fig, axes = plt.subplots(1, 2, figsize=(10, 6), sharey=False)

    for ax, vals, cols, labels, title in [
        (axes[0], d_vals, d_cols, d_labels, 'Debris RPO'),
        (axes[1], r_vals, r_cols, r_labels, 'RH docking cycle'),
    ]:
        bottom = 0.0
        for v, c, lab in zip(vals, cols, labels):
            ax.bar(0, v, bottom=bottom, color=c, width=0.5,
                   edgecolor='white', linewidth=0.5, label=lab)
            if v > 0.08:
                ax.text(0, bottom + v / 2, f'{v:.2f} m/s',
                        ha='center', va='center', fontsize=8,
                        color='white', fontweight='bold')
            bottom += v

        total = sum(vals)
        ax.set_title(f'{title}\nTotal: {total:.2f} m/s', fontsize=10)
        ax.set_xticks([])
        ax.set_ylabel('dV [m/s]')
        ax.set_xlim(-0.5, 0.5)
        ax.legend(fontsize=7.5, loc='upper right',
                  bbox_to_anchor=(1.0, 1.0), framealpha=0.9)
        ax.yaxis.set_major_locator(mticker.MultipleLocator(1.0))
        ax.grid(axis='y', linestyle='--', alpha=0.5)

    fig.suptitle('RPO phase dV breakdown', fontsize=12, fontweight='bold', y=1.01)
    fig.tight_layout()

    path_out = f'{SAVE_DIR}\\rpo_fig3_dv_breakdown.png'
    if save:
        fig.savefig(path_out, dpi=150, bbox_inches='tight')
        print(f'Saved -> {path_out}')
    return fig


# =============================================================================
# FIGURE 4 -- TUMBLE-RATE SENSITIVITY
# =============================================================================

def fig_tumble_sweep(save=True):
    """dV vs target tumble rate for DV_RPO_DEBRIS and DV_RPO_DETUMBLE."""
    rates = np.linspace(0, 6, 25)
    dv_deb, dv_det = [], []
    for w in rates:
        r = simulate_debris_rpo(omega_dps=float(w))
        dv_deb.append(r['dv_debris'])
        dv_det.append(r['dv_detumble'])

    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax = plt.subplots(figsize=(7, 4))

    ax.plot(rates, dv_deb, '-o', color=prop_cycle[0], ms=4, lw=1.8,
            label='DV_RPO_DEBRIS (incl. +50% abort margin)')
    ax.plot(rates, dv_det, '-s', color=prop_cycle[1], ms=4, lw=1.8,
            label='DV_RPO_DETUMBLE (momentum dumping)')

    ax.axvline(OMEGA_REQ_DPS, color='gray', ls='--', lw=1.2)
    ymax = max(max(dv_deb), max(dv_det))
    ax.text(OMEGA_REQ_DPS + 0.08, ymax * 0.55,
            f'REQ-MIS-06 ceiling\n({OMEGA_REQ_DPS:.0f} deg/s, 1 rpm)',
            fontsize=8, color='gray', va='center')

    ax.set_xlabel('Target tumble rate [deg/s]')
    ax.set_ylabel('dV [m/s]')
    ax.set_title('Debris RPO -- tumble-rate sensitivity', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(linestyle='--', alpha=0.5)
    ax.set_xlim(0, 6.3)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    path_out = f'{SAVE_DIR}\\rpo_fig4_tumble_sweep.png'
    if save:
        fig.savefig(path_out, dpi=150, bbox_inches='tight')
        print(f'Saved -> {path_out}')
    return fig


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == '__main__':
    # Pre-compute shared results so debris sim only runs once
    debris_res = simulate_debris_rpo(omega_dps=OMEGA_REQ_DPS)
    rh_res     = simulate_rh_cycle(1884.0)

    fig_debris_traj()
    fig_rh_meet()
    fig_rh_carry()
    fig_rh_dock()
    fig_combined()
    fig_dv_breakdown(debris_res, rh_res)
    fig_tumble_sweep()
    print('All figures saved.')
