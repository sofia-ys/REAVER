"""
REAVER RPO — Plots
==================
Report-ready figures for the close-range RPO simulation:

  1. Debris-RPO trajectory in the target Hill frame (with keep-out spheres and
     the final approach cone).
  2. RH docking-cycle trajectory (cooperative, fixed docking axis).
  3. Phase-by-phase ΔV breakdown — debris RPO vs one RH cycle.
  4. Tumble-rate sensitivity of DV_RPO_DEBRIS and DV_RPO_DETUMBLE.

    cd c:\\Projects\\DSE\\REAVER
    python trajectory_optimizer/rpo_plots.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from rpo_config import (N_MEAN, R_CUTOFF, H_BAR_STANDOFF, R_KOS1, R_KOS2,
                        APPROACH_CONE_DEG, V_MAX, OMEGA_REQ_DPS)
from rpo_dynamics import cw_sample_path
from rpo_debris_sim import simulate_debris_rpo
from rpo_rh_sim import simulate_rh_cycle, R_MEET_START, R_CAPTURE, R_PORT

# ── Dark theme (matches nominal_mission.py) ──────────────────────────────────
BG = '#0d1117'; CB = '#161b22'; TC = '#c9d1d9'; MU_ = '#8b949e'; GR = '#21262d'
A1 = '#58a6ff'; A2 = '#3fb950'; A3 = '#f78166'; A4 = '#d2a8ff'; A5 = '#ffa657'

SAVE_PATH = r'C:\Projects\DSE\REAVER\trajectory_optimizer\rpo_dashboard.png'


def _style(ax, title=''):
    ax.set_facecolor(CB)
    for s in ax.spines.values():
        s.set_edgecolor(GR)
    ax.tick_params(colors=MU_, labelsize=8)
    ax.xaxis.label.set_color(MU_); ax.yaxis.label.set_color(MU_)
    ax.grid(True, color=GR, lw=0.5, alpha=0.7)
    if title:
        ax.set_title(title, color=TC, fontsize=10, fontweight='bold', pad=8)


def _style3d(ax, title=''):
    ax.set_facecolor(CB)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.label.set_color(MU_)
        axis.set_pane_color((0.09, 0.11, 0.13, 1.0))
        axis._axinfo['grid']['color'] = GR
    ax.tick_params(colors=MU_, labelsize=7)
    if title:
        ax.set_title(title, color=TC, fontsize=10, fontweight='bold', pad=2)


def _transfer_time(r0, rf):
    chord = np.linalg.norm(np.asarray(rf, float) - np.asarray(r0, float))
    return max(chord / (0.5 * V_MAX), 1.0)


def _path_through(waypoints):
    """Concatenate CW-sampled paths through a list of Hill-frame waypoints."""
    segs = []
    for a, b in zip(waypoints[:-1], waypoints[1:]):
        a = np.asarray(a, float); b = np.asarray(b, float)
        segs.append(cw_sample_path(a, b, N_MEAN, _transfer_time(a, b), n_pts=60))
    return np.vstack(segs)


def _wire_sphere(ax, r, color, alpha=0.25):
    u = np.linspace(0, 2 * np.pi, 18)
    v = np.linspace(0, np.pi, 12)
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(x, y, z, color=color, lw=0.4, alpha=alpha)


def _approach_cone(ax, axis_dir, half_deg, length, color):
    """Draw a few generatrix lines of the approach cone about ``axis_dir``."""
    axis_dir = np.asarray(axis_dir, float)
    axis_dir = axis_dir / np.linalg.norm(axis_dir)
    # build two vectors perpendicular to the axis
    tmp = np.array([1.0, 0.0, 0.0]) if abs(axis_dir[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(axis_dir, tmp); e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis_dir, e1)
    a = np.radians(half_deg)
    for phi in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        d = np.cos(a) * axis_dir + np.sin(a) * (np.cos(phi) * e1 + np.sin(phi) * e2)
        pts = np.outer(np.linspace(0, length, 2), d)
        ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=color, lw=0.5, alpha=0.5)


# =============================================================================
# PANELS
# =============================================================================

def panel_debris_traj(ax):
    _style3d(ax, 'Debris RPO trajectory  (target Hill frame)')
    # Close-range approach waypoints (the long cut-off leg is annotated, not drawn,
    # so the 5 m / 1.35 m keep-out spheres stay visible).
    wp = [
        [0.0, 0.0, H_BAR_STANDOFF],   # H-bar standoff (20 m)
        [0.0, 0.0, R_KOS1],           # KOS1 (interaction distance)
        [0.0, 0.0, R_KOS2],           # KOS2 (final)
    ]
    path = _path_through(wp)
    ax.plot(path[:, 0], path[:, 1], path[:, 2], color=A1, lw=2.0, label='approach')
    ax.scatter(*np.array(wp).T, color=A5, s=30, depthshade=False)
    ax.scatter(0, 0, 0, color=A3, s=60, marker='*', depthshade=False, label='target CoM')

    _wire_sphere(ax, R_KOS1, A4, 0.30)
    _wire_sphere(ax, R_KOS2, A3, 0.40)
    _approach_cone(ax, [0, 0, 1], APPROACH_CONE_DEG, R_KOS1, A2)

    ax.set_xlabel('R-bar [m]'); ax.set_ylabel('V-bar [m]'); ax.set_zlabel('H-bar [m]')
    ax.set_xlim(-12, 12); ax.set_ylim(-12, 12); ax.set_zlim(0, 22)
    ax.text(0, 0, 21.5, 'from cut-off 150 m', color=MU_, fontsize=7)
    ax.legend(fontsize=7, facecolor=CB, labelcolor=TC, edgecolor=GR, loc='upper left')


def panel_rh_traj(ax):
    _style3d(ax, 'RH docking-cycle trajectory  (RH Hill frame)')
    wp = [R_MEET_START, R_CAPTURE, R_PORT]
    path = _path_through(wp)
    ax.plot(path[:, 0], path[:, 1], path[:, 2], color=A2, lw=2.0, label='MEET + DOCK')
    ax.scatter(*np.array(wp).T, color=A5, s=30, depthshade=False)
    ax.scatter(0, 0, 0, color=A1, s=70, marker='s', depthshade=False, label='RH port')
    ax.text(0, -500, 0, '500 m hold', color=MU_, fontsize=7)
    ax.text(0, -15, 0, 'capture ~15 m', color=MU_, fontsize=7)
    ax.set_xlabel('R-bar [m]'); ax.set_ylabel('V-bar [m]'); ax.set_zlabel('H-bar [m]')
    ax.set_ylim(-520, 40)
    ax.legend(fontsize=7, facecolor=CB, labelcolor=TC, edgecolor=GR, loc='upper left')


def panel_phase_bars(ax, debris_res, rh_res):
    _style(ax, 'Phase ΔV breakdown')
    # Debris stack (exclude detumble from the DEBRIS bar; show it stacked on top
    # in a distinct colour to highlight it is a separate constant).
    d_phases = debris_res['phases']                  # P1,P2,P3,P4,P5app,P5det
    d_labels = ['P1 appr', 'P2 dwell', 'P3 spin', 'P4 arm', 'P5 appr', 'P5 detumble']
    d_vals = [p.dv for p in d_phases]
    d_cols = [A1, MU_, A4, A5, A1, A3]

    r_phases = rh_res['phases']                       # A1,A2,A3,B1c,B1p,B2
    r_labels = ['A1 hold', 'A2 appr', 'A3 clamp', 'B1 carry', 'B1 point', 'B2 dock']
    r_vals = [p.dv for p in r_phases]
    r_cols = [MU_, A2, A4, A2, A5, A1]

    # Two stacked bars side by side.
    for x, vals, cols, labels, title in [
        (0, d_vals, d_cols, d_labels, 'Debris RPO'),
        (1, r_vals, r_cols, r_labels, 'RH cycle'),
    ]:
        bottom = 0.0
        for v, c, lab in zip(vals, cols, labels):
            ax.bar(x, v, bottom=bottom, color=c, width=0.55, ec=BG, lw=0.5)
            if v > 0.15:
                ax.text(x, bottom + v / 2, lab, ha='center', va='center',
                        fontsize=6.5, color=BG, fontweight='bold')
            bottom += v
        ax.text(x, bottom + 0.1, f'{title}\nΣ {bottom:.2f} m/s', ha='center',
                va='bottom', fontsize=8, color=TC, fontweight='bold')

    ax.set_xticks([0, 1]); ax.set_xticklabels(['Debris RPO', 'RH cycle'], color=TC)
    ax.set_ylabel('ΔV [m/s]')
    ax.set_xlim(-0.6, 1.6)


def panel_tumble_sweep(ax):
    _style(ax, 'Tumble-rate sensitivity')
    rates = np.linspace(0, 6, 13)
    dv_deb = []; dv_det = []
    for w in rates:
        r = simulate_debris_rpo(omega_dps=float(w))
        dv_deb.append(r['dv_debris']); dv_det.append(r['dv_detumble'])
    ax.plot(rates, dv_deb, '-o', color=A1, ms=4, lw=1.8, label='DV_RPO_DEBRIS')
    ax.plot(rates, dv_det, '-s', color=A3, ms=4, lw=1.8, label='DV_RPO_DETUMBLE')
    ax.axvline(OMEGA_REQ_DPS, color=A5, ls='--', lw=1.2, alpha=0.8)
    ax.text(OMEGA_REQ_DPS - 0.1, ax.get_ylim()[1] * 0.5, '1 rpm (REQ-MIS-06)',
            color=A5, fontsize=7, rotation=90, va='center', ha='right')
    ax.set_xlabel('target tumble rate [deg/s]'); ax.set_ylabel('ΔV [m/s]')
    ax.legend(fontsize=7, facecolor=CB, labelcolor=TC, edgecolor=GR)


# =============================================================================
# DASHBOARD
# =============================================================================

def make_dashboard(save_path=SAVE_PATH):
    debris_res = simulate_debris_rpo(omega_dps=OMEGA_REQ_DPS)
    rh_res = simulate_rh_cycle(1884.0)   # representative heaviest cycle

    fig = plt.figure(figsize=(18, 11))
    fig.patch.set_facecolor(BG)
    gs = GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.22)

    fig.text(0.5, 0.965, 'REAVER — Close-Range RPO Dashboard', ha='center',
             color=TC, fontsize=16, fontweight='bold')
    fig.text(0.5, 0.938,
             f"DV_RPO_DEBRIS {debris_res['dv_debris']:.2f}  |  "
             f"DV_RPO_DETUMBLE {debris_res['dv_detumble']:.2f}  |  "
             f"DV_RPO_RH {rh_res['dv_rh']:.2f}  m/s", ha='center',
             color=MU_, fontsize=9)

    panel_debris_traj(fig.add_subplot(gs[0, 0], projection='3d'))
    panel_rh_traj(fig.add_subplot(gs[0, 1], projection='3d'))
    panel_phase_bars(fig.add_subplot(gs[1, 0]), debris_res, rh_res)
    panel_tumble_sweep(fig.add_subplot(gs[1, 1]))

    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor=BG)
    print(f"Dashboard saved -> {save_path}")


if __name__ == '__main__':
    make_dashboard()
