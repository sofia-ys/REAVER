"""
REAVER RPO — Equation Reference PDF
=====================================
Generates rpo_equations.pdf: a formatted summary of every equation used in
the RPO simulation (rpo_dynamics.py, rpo_guidance.py, rpo_control.py).
Intended as a reference for the final report.

    python trajectory_optimizer/rpo_equations.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

SAVE_PATH = r'C:\Projects\DSE\REAVER\trajectory_optimizer\rpo_equations.pdf'

# ── Mathtext normaliser ───────────────────────────────────────────────────────
# matplotlib mathtext does not support every LaTeX command.
# This function normalises strings to the supported subset before rendering.
import re as _re

def _m(s):
    """Normalise a mathtext string to matplotlib-supported commands."""
    s = s.replace(r'\tfrac', r'\frac')
    s = s.replace(r'\boldsymbol', r'\mathbf')
    s = s.replace(r'^\top', r'^{T}')
    s = s.replace(r'^\intercal', r'^{T}')
    s = s.replace(r'\bigl', r'\left')
    s = s.replace(r'\bigr', r'\right')
    s = s.replace(r'\mathbb{R}', r'R')
    # \mathrm{word} → word (roman text in subscripts etc.)
    s = _re.sub(r'\\mathrm\{([^}]+)\}', r'\1', s)
    return s

# ── Style ────────────────────────────────────────────────────────────────────
TITLE_FS  = 13
HEAD_FS   = 11
EQ_FS     = 10.5
NOTE_FS   = 9
LABEL_FS  = 9.5
C_HEAD    = '#1a3a5c'     # dark navy for section headings
C_EQ      = '#1a1a1a'     # near-black for equations
C_NOTE    = '#555555'     # grey for notes / parameter definitions
C_BOX     = '#eef3f8'     # light blue box fill
C_LINE    = '#2a6099'     # rule colour

def new_fig():
    fig, ax = plt.subplots(figsize=(8.27, 11.69))   # A4
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis('off')
    return fig, ax

def hline(ax, y, lw=0.8, color=C_LINE):
    ax.axhline(y, xmin=0.05, xmax=0.95, lw=lw, color=color)

def heading(ax, y, text, fs=HEAD_FS):
    ax.text(0.05, y, text, transform=ax.transAxes,
            fontsize=fs, fontweight='bold', color=C_HEAD, va='top')

def eq(ax, y, lhs, rhs, x_lhs=0.08, x_rhs=0.22):
    """Render one equation: LHS (plain) and RHS (mathtext)."""
    ax.text(x_lhs, y, _m(lhs), transform=ax.transAxes,
            fontsize=EQ_FS, color=C_EQ, va='top', family='monospace')
    ax.text(x_rhs, y, _m(rhs), transform=ax.transAxes,
            fontsize=EQ_FS, color=C_EQ, va='top',
            usetex=False)   # matplotlib mathtext, no LaTeX needed

def note(ax, y, text, x=0.08):
    ax.text(x, y, _m(text), transform=ax.transAxes,
            fontsize=NOTE_FS, color=C_NOTE, va='top', style='italic')

def param_block(ax, y, lines, x=0.10):
    """Print a block of parameter definitions."""
    dy = 0.032
    for line in lines:
        ax.text(x, y, line, transform=ax.transAxes,
                fontsize=NOTE_FS, color=C_NOTE, va='top')
        y -= dy
    return y

def box(ax, y_top, y_bot, x0=0.05, x1=0.95):
    rect = mpatches.FancyBboxPatch((x0, y_bot), x1 - x0, y_top - y_bot,
                                   boxstyle='round,pad=0.008',
                                   linewidth=0.6, edgecolor=C_LINE,
                                   facecolor=C_BOX, transform=ax.transAxes,
                                   zorder=0)
    ax.add_patch(rect)

def footer(ax, page_num, total):
    ax.text(0.5, 0.015, f'REAVER RPO — Equation Reference     page {page_num} of {total}',
            transform=ax.transAxes, fontsize=8, color='#888888',
            ha='center', va='bottom')
    hline(ax, 0.025, lw=0.5)

# =============================================================================
# PAGE DEFINITIONS
# =============================================================================

def page_title(pdf):
    fig, ax = new_fig()
    ax.text(0.5, 0.72, 'REAVER RPO Simulation',
            transform=ax.transAxes, fontsize=20, fontweight='bold',
            color=C_HEAD, ha='center', va='center')
    ax.text(0.5, 0.65, 'Equation Reference',
            transform=ax.transAxes, fontsize=16, color=C_HEAD,
            ha='center', va='center')
    hline(ax, 0.62, lw=1.5)
    ax.text(0.5, 0.56,
            'Summary of all equations used in rpo\\_dynamics.py,\n'
            'rpo\\_guidance.py, and rpo\\_control.py.',
            transform=ax.transAxes, fontsize=10.5, color=C_NOTE,
            ha='center', va='center', linespacing=1.8)

    toc = [
        ('1', 'Orbital Mechanics — Mean Motion'),
        ('2', 'Relative Translational Dynamics (Hill / CW)'),
        ('3', 'CW State-Transition Matrix'),
        ('4', 'Two-Impulse CW Targeting & Transfer Time'),
        ('5', 'GEO Perturbations & Stationkeeping'),
        ('6', 'Rigid-Body Rotational Dynamics'),
        ('7', 'Quaternion Kinematics'),
        ('8', 'Body Geometry & Inertia Tensors'),
        ('9', 'Attitude Control Budget (Angular Impulse Method)'),
        ('10', 'Propellant & Tsiolkovsky Equation'),
    ]
    y = 0.48
    ax.text(0.5, y + 0.04, 'Contents', transform=ax.transAxes,
            fontsize=11, fontweight='bold', color=C_HEAD, ha='center')
    for num, title in toc:
        ax.text(0.22, y, f'{num}.', transform=ax.transAxes,
                fontsize=9.5, color=C_HEAD, ha='right')
        ax.text(0.24, y, title, transform=ax.transAxes,
                fontsize=9.5, color=C_EQ)
        y -= 0.038
    footer(ax, 1, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_mean_motion(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '1.  Orbital Mechanics — Mean Motion')
    hline(ax, 0.915)

    note(ax, 0.885,
         'All CW relative-motion equations are parameterised by the mean motion n of '
         'the reference (target) orbit.\n'
         'At GEO (a = 42 878 km) n is approximately 7.11 × 10⁻⁵ rad/s, '
         'roughly 15× smaller than a typical LEO value.')

    box(ax, 0.82, 0.67)
    heading(ax, 0.805, 'Mean motion of a circular orbit', fs=LABEL_FS)
    eq(ax, 0.765,
       'n  =',
       r'$n = \sqrt{\dfrac{\mu}{a^3}}$',
       x_lhs=0.08, x_rhs=0.18)
    param_block(ax, 0.695, [
        r'$\mu$  = 3.986 × 10¹⁴  m³/s²   (Earth gravitational parameter)',
        r'$a$   = semi-major axis of the target orbit  [m]',
        r'$n$   = mean motion  [rad/s]',
    ], x=0.25)

    note(ax, 0.635,
         'REAVER values:   a = 42 878 km  →  n = 7.11 × 10⁻⁵ rad/s\n'
         '                 Orbital period T = 2π/n ≈ 88 341 s ≈ 24.5 h')

    footer(ax, 2, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_cw(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '2.  Relative Translational Dynamics (Hill / CW Equations)')
    hline(ax, 0.915)

    note(ax, 0.885,
         'The Clohessy–Wiltshire (Hill) equations describe the relative motion of the chaser\n'
         'with respect to the target in the local Hill frame:\n'
         '   x  = radial (R-bar),   y  = along-track (V-bar),   z  = cross-track (H-bar).\n'
         'Assumptions: target on a circular orbit; inter-satellite distance ≪ orbital radius.')

    box(ax, 0.815, 0.57)
    heading(ax, 0.800, 'Hill / CW equations of relative motion', fs=LABEL_FS)
    eq(ax, 0.768, r'$\ddot{x}$  =',
       r'$\ddot{x} - 2n\dot{y} - 3n^2 x = f_x$', x_rhs=0.20)
    eq(ax, 0.730, r'$\ddot{y}$  =',
       r'$\ddot{y} + 2n\dot{x} = f_y$',              x_rhs=0.20)
    eq(ax, 0.692, r'$\ddot{z}$  =',
       r'$\ddot{z} + n^2 z = f_z$',                  x_rhs=0.20)

    param_block(ax, 0.648, [
        r'$x, y, z$           = relative position components [m]',
        r'$f_x, f_y, f_z$     = specific force (thrust + perturbations) [m/s²]',
        r'$n$                  = mean motion of target orbit [rad/s]',
    ], x=0.25)

    note(ax, 0.530,
         'State vector form used in the simulation:   '
         r'$\mathbf{x} = (x,\, y,\, z,\, \dot{x},\, \dot{y},\, \dot{z})^{\!\top}$')

    box(ax, 0.488, 0.345)
    heading(ax, 0.473, 'Perturbation forces at GEO (replaces LEO drag/magnetic from paper)', fs=LABEL_FS)
    note(ax, 0.442,
         'Atmospheric drag: negligible at GEO altitude (42 878 km).   '
         'Magnetic dipole torque: negligible.\n'
         'Two forces remain:')
    eq(ax, 0.402, '',
       r'$a_\mathrm{SRP} = \dfrac{C_R \cdot (A/m) \cdot P_\odot}{1}$',
       x_lhs=0.08, x_rhs=0.12)
    note(ax, 0.360,
         r'$a_\mathrm{triax}$  from Earth J₂₂ tesseral harmonic (GEO longitude-dependent); '
         'representative value ≈ 3 × 10⁻⁸ m/s²', x=0.12)
    param_block(ax, 0.320, [
        r'$C_R$      = 1.5     (reflectivity, MLI-covered bus)',
        r'$A/m$      = 0.015  m²/kg  (MS area-to-mass ratio, placeholder)',
        r'$P_\odot$  = 4.56 × 10⁻⁶  N/m²  (solar radiation pressure at 1 AU)',
    ], x=0.25)

    footer(ax, 3, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_stm(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '3.  CW State-Transition Matrix (STM)')
    hline(ax, 0.915)

    note(ax, 0.885,
         r'Under zero forcing (free drift),  $\mathbf{x}(t) = \Phi(t)\,\mathbf{x}(0)$  '
         'where Φ(t) is the 6×6 STM.\n'
         'Shorthand:  s = sin(nt),   c = cos(nt),   '
         r'$\mathbf{x} = (x,y,z,\dot{x},\dot{y},\dot{z})^\top$')

    box(ax, 0.840, 0.290)
    heading(ax, 0.825, r'Position-to-position block  $\Phi_{rr}$  (3×3)', fs=LABEL_FS)
    for i,(lhs,rhs) in enumerate([
        (r'$\Phi_{rr}[1,1]$', r'$= 4 - 3c$'),
        (r'$\Phi_{rr}[2,1]$', r'$= 6(s - nt)$'),
        (r'$\Phi_{rr}[3,3]$', r'$= c$'),
        (r'all other entries', r'$= 0$'),
    ]):
        y = 0.790 - i*0.038
        ax.text(0.12, y, _m(lhs), transform=ax.transAxes, fontsize=EQ_FS, color=C_EQ, va='top')
        ax.text(0.40, y, _m(rhs), transform=ax.transAxes, fontsize=EQ_FS, color=C_EQ, va='top')

    heading(ax, 0.630, r'Position-to-velocity block  $\Phi_{rv}$  (3x3)', fs=LABEL_FS)
    for i,(lhs,rhs) in enumerate([
        (r'$\Phi_{rv}[1,1]$', r'$= s/n$'),
        (r'$\Phi_{rv}[1,2]$', r'$= 2(1-c)/n$'),
        (r'$\Phi_{rv}[2,1]$', r'$= 2(c-1)/n$'),
        (r'$\Phi_{rv}[2,2]$', r'$= (4s - 3nt)/n$'),
        (r'$\Phi_{rv}[3,3]$', r'$= s/n$,   all other entries $= 0$'),
    ]):
        y = 0.595 - i*0.038
        ax.text(0.12, y, _m(lhs), transform=ax.transAxes, fontsize=EQ_FS, color=C_EQ, va='top')
        ax.text(0.40, y, _m(rhs), transform=ax.transAxes, fontsize=EQ_FS, color=C_EQ, va='top')

    heading(ax, 0.405, r'Velocity-to-position block  $\Phi_{vr}$  (3x3)', fs=LABEL_FS)
    for i,(lhs,rhs) in enumerate([
        (r'$\Phi_{vr}[1,1]$', r'$= 3ns$'),
        (r'$\Phi_{vr}[2,1]$', r'$= 6n(c-1)$'),
        (r'$\Phi_{vr}[3,3]$', r'$= -ns$,   all other entries $= 0$'),
    ]):
        y = 0.370 - i*0.038
        ax.text(0.12, y, _m(lhs), transform=ax.transAxes, fontsize=EQ_FS, color=C_EQ, va='top')
        ax.text(0.40, y, _m(rhs), transform=ax.transAxes, fontsize=EQ_FS, color=C_EQ, va='top')

    note(ax, 0.255,
         r'Velocity-to-velocity block  $\Phi_{vv}$:  '
         r'$\Phi_{vv}[1,1]=c,\ \Phi_{vv}[1,2]=2s,\ \Phi_{vv}[2,1]=-2s,$'
         '\n'
         r'    $\Phi_{vv}[2,2]=4c-3,\ \Phi_{vv}[3,3]=c$,  all other entries = 0')
    note(ax, 0.205,
         'Implementation: rpo_dynamics.cw_stm(n, t)  →  full 6×6 matrix.\n'
         'Used for: free-drift propagation, two-impulse targeting, path sampling for plots.')

    footer(ax, 4, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_targeting(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '4.  Two-Impulse CW Targeting & Transfer Time')
    hline(ax, 0.915)

    note(ax, 0.885,
         'For each translational phase, find the two burns that transfer the chaser\n'
         r'from rest at $\mathbf{r}_0$ to rest at $\mathbf{r}_f$ in time $t$.')

    box(ax, 0.848, 0.56)
    heading(ax, 0.833, 'Two-impulse CW targeting', fs=LABEL_FS)

    eq(ax, 0.800,
       'Step 1:',
       r'$\mathbf{v}_0^+ = \Phi_{rv}^{-1}\!\left(\mathbf{r}_f - \Phi_{rr}\,\mathbf{r}_0\right)$',
       x_rhs=0.22)
    note(ax, 0.760,
         r'Solve for the required velocity just after the first burn.  $\Phi_{rv}$ must be invertible '
         r'(non-singular transfer).', x=0.22)

    eq(ax, 0.718,
       'Step 2:',
       r'$\Delta\mathbf{v}_1 = \mathbf{v}_0^+ - \mathbf{v}_0$',
       x_rhs=0.22)
    note(ax, 0.682, r'First burn: difference between required and current velocity.', x=0.22)

    eq(ax, 0.646,
       'Step 3:',
       r'$\mathbf{v}_f^- = \Phi_{vr}\,\mathbf{r}_0 + \Phi_{vv}\,\mathbf{v}_0^+$',
       x_rhs=0.22)
    note(ax, 0.610, r'Velocity at arrival, just before the second burn.', x=0.22)

    eq(ax, 0.574,
       'Step 4:',
       r'$\Delta\mathbf{v}_2 = \mathbf{v}_f - \mathbf{v}_f^-$',
       x_rhs=0.22)

    eq(ax, 0.535, r'$\Delta V$  =',
       r'$\Delta V = |\Delta\mathbf{v}_1| + |\Delta\mathbf{v}_2|$',
       x_rhs=0.22)

    box(ax, 0.498, 0.32)
    heading(ax, 0.484, 'Transfer time sizing', fs=LABEL_FS)
    eq(ax, 0.450, 'chord =',
       r'$d = |\mathbf{r}_f - \mathbf{r}_0|$',
       x_rhs=0.24)
    eq(ax, 0.412, 't  =',
       r'$t = \dfrac{d}{0.5\; V_\mathrm{max}}$',
       x_rhs=0.24)
    note(ax, 0.368,
         r'The transfer time is sized so the mean speed equals $0.5\,V_\mathrm{max}$, '
         'ensuring the velocity constraint\n'
         r'$|\mathbf{v}(t)| \leq V_\mathrm{max} = 1.0$ m/s (low-impact docking requirement) is met with margin.',
         x=0.24)
    param_block(ax, 0.312, [
        r'$V_\mathrm{max}$ = 1.0 m/s   (maximum approach velocity)',
        r'$t_\mathrm{min}$ = 1.0 s      (hard lower bound to avoid STM singularity)',
    ], x=0.25)

    footer(ax, 5, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_perturbations(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '5.  GEO Perturbations & Stationkeeping ΔV')
    hline(ax, 0.915)

    note(ax, 0.885,
         'During hold / dwell phases (P2, RH A1) the chaser must null the accumulated\n'
         'drift from two persistent GEO perturbations.  No net translation occurs.')

    box(ax, 0.848, 0.645)
    heading(ax, 0.833, 'Solar Radiation Pressure (SRP)', fs=LABEL_FS)
    eq(ax, 0.800, r'$a_\mathrm{SRP}$  =',
       r'$a_\mathrm{SRP} = C_R \cdot \dfrac{A}{m} \cdot P_\odot$',
       x_rhs=0.30)
    param_block(ax, 0.748, [
        r'$C_R$     = 1.5               reflectivity coefficient (MLI)',
        r'$A/m$     = 0.015 m²/kg       area-to-mass ratio (placeholder)',
        r'$P_\odot$ = 4.56 × 10⁻⁶ N/m² solar pressure at 1 AU',
    ], x=0.25)
    note(ax, 0.683, r'Result:  $a_\mathrm{SRP} \approx 1.03 \times 10^{-7}$ m/s²', x=0.25)

    box(ax, 0.640, 0.490)
    heading(ax, 0.625, 'Earth Triaxiality (J₂₂ tesseral harmonic)', fs=LABEL_FS)
    note(ax, 0.594,
         'A longitude-dependent tangential acceleration from Earth\'s equatorial\n'
         'ellipticity.  Represented as a constant bias per phase.',
         x=0.25)
    eq(ax, 0.552, r'$a_\mathrm{triax}$  ≈',
       r'$a_\mathrm{triax} \approx 3.0 \times 10^{-8}$ m/s$^2$',
       x_rhs=0.34)

    box(ax, 0.452, 0.260)
    heading(ax, 0.438, 'Stationkeeping ΔV (dwell phases)', fs=LABEL_FS)
    note(ax, 0.408,
         'First-order model: constant disturbance integrates into drift velocity\n'
         'that must be cancelled once per dwell.', x=0.25)
    eq(ax, 0.368, r'$a_\mathrm{dist}$  =',
       r'$a_\mathrm{dist} = a_\mathrm{SRP} + a_\mathrm{triax}$',
       x_rhs=0.34)
    eq(ax, 0.325, r'$\Delta V_\mathrm{SK}$  =',
       r'$\Delta V_\mathrm{SK} = a_\mathrm{dist} \cdot t_\mathrm{dwell}$',
       x_rhs=0.34)
    note(ax, 0.278,
         r'$t_\mathrm{dwell}$ = 28 800 s (8 h, confirmed)  →  '
         r'$\Delta V_\mathrm{SK} \approx 0.004$ m/s  (negligible for ΔV budget,\n'
         '   significant for T_OPS update)', x=0.25)

    footer(ax, 6, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_euler(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '6.  Rigid-Body Rotational Dynamics (Euler\'s Equation)')
    hline(ax, 0.915)

    note(ax, 0.885,
         'Used to propagate the tumbling state of the target debris.\n'
         'The chaser\'s own attitude is controlled, not freely tumbling.')

    box(ax, 0.848, 0.620)
    heading(ax, 0.833, 'Euler\'s rigid-body equation', fs=LABEL_FS)
    eq(ax, 0.798, r'$\dot{\boldsymbol{\omega}}$  =',
       r'$\dot{\boldsymbol{\omega}} = \mathbf{J}^{-1}\!\left(\mathbf{M}_\mathrm{net} - \boldsymbol{\omega} \times \mathbf{J}\boldsymbol{\omega}\right)$',
       x_rhs=0.26)
    param_block(ax, 0.748, [
        r'$\boldsymbol{\omega}$     = angular velocity vector of the body  [rad/s]',
        r'$\mathbf{J}$              = instantaneous inertia tensor  [kg·m²]',
        r'$\mathbf{M}_\mathrm{net}$ = net external torque vector  [N·m]',
        r'                           = 0 for free tumble (target debris)',
    ], x=0.25)
    note(ax, 0.645,
         'The target is modelled as free-tumbling ('
         r'$\mathbf{M}_\mathrm{net} = \mathbf{0}$'
         ') at the requirement\nceiling of 1 rpm (6 deg/s) about an arbitrary axis.',
         x=0.25)

    note(ax, 0.590,
         'Integration: scipy.integrate.solve_ivp with RK45,\n'
         '    rtol = 10⁻⁹,  atol = 10⁻¹²,  quaternion renormalised every step.')

    box(ax, 0.553, 0.380)
    heading(ax, 0.538, 'Conservation check — rotational kinetic energy', fs=LABEL_FS)
    note(ax, 0.508,
         'Under free tumble with no external torques, rotational KE is conserved.\n'
         'Used as a numerical verification of the integrator.', x=0.25)
    eq(ax, 0.466, r'$T_\mathrm{rot}$  =',
       r'$T_\mathrm{rot} = \tfrac{1}{2}\,\boldsymbol{\omega}^\top \mathbf{J}\,\boldsymbol{\omega} = \mathrm{const}$',
       x_rhs=0.28)
    note(ax, 0.420,
         r'Test passes if $\sigma(T_\mathrm{rot}) / \mu(T_\mathrm{rot}) < 10^{-3}$  '
         'over the propagation interval.', x=0.25)

    footer(ax, 7, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_quat(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '7.  Quaternion Kinematics')
    hline(ax, 0.915)

    note(ax, 0.885,
         r'Attitude is represented by the unit quaternion  $\mathbf{q} = (q_1, q_2, q_3, q_4)^\top$'
         '\nwhere  q₄  is the scalar part.  Propagated alongside the Euler equation.')

    box(ax, 0.848, 0.570)
    heading(ax, 0.833, r'Quaternion kinematic equation  (Conings & Mooij Eq. 3)', fs=LABEL_FS)
    eq(ax, 0.798, r'$\dot{\mathbf{q}}$  =',
       r'$\dot{\mathbf{q}} = \dfrac{1}{2}\,\Xi(\mathbf{q})\,\boldsymbol{\omega}$',
       x_rhs=0.22)

    heading(ax, 0.740, r'$\Xi(\mathbf{q})$  — the 4×3 kinematic matrix', fs=LABEL_FS)
    note(ax, 0.710,
         r'Row 1:  $[\ q_4,\ -q_3,\  q_2\ ]$'
         '\n'
         r'Row 2:  $[\ q_3,\  q_4,\ -q_1\ ]$'
         '\n'
         r'Row 3:  $[\ -q_2,\ q_1,\  q_4\ ]$'
         '\n'
         r'Row 4:  $[\ -q_1,\ -q_2,\ -q_3\ ]$',
         x=0.18)

    note(ax, 0.578,
         r'$\Xi(\mathbf{q})$ is 4×3.  The product $\Xi\,\boldsymbol{\omega}$ yields '
         r'a 4-vector $\dot{\mathbf{q}}$.')

    note(ax, 0.536,
         'After each integration step the quaternion is renormalised:\n'
         r'    $\mathbf{q} \leftarrow \mathbf{q} / |\mathbf{q}|$'
         '   to prevent numerical drift from the unit-norm constraint.')

    box(ax, 0.498, 0.300)
    heading(ax, 0.483, 'Combined state vector for target tumble propagation', fs=LABEL_FS)
    note(ax, 0.452,
         r'State:  $\mathbf{y} = (\mathbf{q},\, \boldsymbol{\omega})^\top \in \mathbb{R}^7$',
         x=0.20)
    eq(ax, 0.412,
       r'$\dot{q}_1,\dot{q}_2,\dot{q}_3,\dot{q}_4$  =',
       r'$\tfrac{1}{2}\,\Xi(\mathbf{q})\,\boldsymbol{\omega}$',
       x_lhs=0.10, x_rhs=0.50)
    eq(ax, 0.370,
       r'$\dot{\boldsymbol{\omega}}$  =',
       r'$\mathbf{J}^{-1}\!\left(\mathbf{M}_\mathrm{net} - \boldsymbol{\omega}\times\mathbf{J}\boldsymbol{\omega}\right)$',
       x_lhs=0.10, x_rhs=0.50)

    footer(ax, 8, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_inertia(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '8.  Body Geometry & Inertia Tensors')
    hline(ax, 0.915)

    box(ax, 0.900, 0.728)
    heading(ax, 0.885, 'Uniform rectangular box (used for MS, tug, debris bus)', fs=LABEL_FS)
    eq(ax, 0.852, r'$I_{xx}$  =',
       r'$I_{xx} = \dfrac{m}{12}(b^2+c^2),\quad I_{yy} = \dfrac{m}{12}(a^2+c^2),\quad I_{zz} = \dfrac{m}{12}(a^2+b^2)$',
       x_rhs=0.24)
    param_block(ax, 0.802, [
        r'$a, b, c$  = edge lengths along x, y, z  [m]',
        r'MS:  3.0 × 3.0 × 3.5 m,   Tug: 0.8 × 0.8 × 0.8 m  (both confirmed)',
        r'Debris bus: 2.5 × 2.5 × 3.0 m  (placeholder)',
    ], x=0.25)

    box(ax, 0.695, 0.490)
    heading(ax, 0.680, 'Debris inertia with deployed solar panels', fs=LABEL_FS)
    note(ax, 0.650,
         'Panels modelled as point masses at ±half-span along the body y-axis.\n'
         'They inflate Ixx and Izz while leaving Iyy near the bus value:', x=0.20)
    eq(ax, 0.608, r'$\mathbf{J}_\mathrm{debris}$  =',
       r'$\mathbf{J}_\mathrm{bus} + m_\mathrm{panel} \cdot \left(\dfrac{s}{2}\right)^{\!2} \cdot \mathrm{diag}(1,\,0,\,1)$',
       x_rhs=0.38)
    param_block(ax, 0.555, [
        r'$m_\mathrm{panel}$ = 0.15 × m_debris  (both panels combined, 15% mass fraction)',
        r'$s$ = 17.2 m  (COMSATBW-1, BSS-702, confirmed solar panel span)',
    ], x=0.25)

    box(ax, 0.455, 0.260)
    heading(ax, 0.440, 'Parallel-axis theorem (combined body assembly)', fs=LABEL_FS)
    note(ax, 0.410,
         'To shift a body\'s inertia from its own CoM to a reference point:', x=0.20)
    eq(ax, 0.368, r'$\mathbf{J}_\mathrm{ref}$  =',
       r'$\mathbf{J}_\mathrm{ref} = \mathbf{J}_\mathrm{cm} + m\left(|\mathbf{d}|^2\mathbf{I} - \mathbf{d}\otimes\mathbf{d}\right)$',
       x_rhs=0.28)
    note(ax, 0.322,
         r'$\mathbf{d}$ = vector from reference point to body CoM.  '
         'Applied to every sub-body\n'
         '   to build the combined (MS + arm + tug + debris) inertia for the detumble phase.', x=0.25)

    footer(ax, 9, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_attitude_control(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '9.  Attitude Control Budget (Angular Impulse Method)')
    hline(ax, 0.915)

    note(ax, 0.885,
         'Reaction wheels carry smooth tracking torque but cannot hold net spin indefinitely.\n'
         'The net angular momentum change is therefore delivered (or dumped) by the RCS thrusters.\n'
         'This gives a conservative (all-thruster) propellant bound — correct posture for sizing.')

    box(ax, 0.848, 0.600)
    heading(ax, 0.833, 'Angular momentum and torque', fs=LABEL_FS)
    eq(ax, 0.800, r'$H$  =',
       r'$H = |\mathbf{J}\,\boldsymbol{\omega}|$',
       x_rhs=0.22)
    eq(ax, 0.758, r'$\tau$  =',
       r'$\tau = 2\,F\,L$',
       x_rhs=0.22)
    note(ax, 0.718,
         r'Opposed thruster pair: two thrusters each at force $F$, '
         r'separated by moment arm $L$ → torque $\tau = 2FL$.', x=0.22)
    param_block(ax, 0.675, [
        r'$F$ = 64 N  (RCS thruster force, confirmed)',
        r'$L$ = 1.75 m  (effective moment arm, confirmed)',
    ], x=0.25)

    box(ax, 0.588, 0.360)
    heading(ax, 0.573, 'Propellant and equivalent ΔV', fs=LABEL_FS)
    eq(ax, 0.538, r'$m_\mathrm{prop}$  =',
       r'$m_\mathrm{prop} = \dfrac{H}{L \cdot v_\mathrm{ex}}$',
       x_rhs=0.28)
    eq(ax, 0.490, r'$\Delta V_\mathrm{att}$  =',
       r'$\Delta V_\mathrm{att} = v_\mathrm{ex}\,\ln\!\left(\dfrac{m_0}{m_0 - m_\mathrm{prop}}\right)$',
       x_rhs=0.28)
    eq(ax, 0.438, r'$t_\mathrm{burn}$  =',
       r'$t_\mathrm{burn} = \dfrac{H}{2\,F\,L}$',
       x_rhs=0.28)

    note(ax, 0.388, r'Phase 4 (arm extension) uses the momentum delta between two configurations:', x=0.20)
    box(ax, 0.360, 0.245)
    eq(ax, 0.332, r'$\Delta H$  =',
       r'$\Delta H = \bigl|\,H_\mathrm{ready} - H_\mathrm{stowed}\,\bigr|$',
       x_rhs=0.28)
    note(ax, 0.292,
         r'$H_\mathrm{stowed}$ = $|\mathbf{J}_\mathrm{MS}\,\boldsymbol{\omega}|$   '
         r'(tug not yet extended)'
         '\n'
         r'$H_\mathrm{ready}$  = $|\mathbf{J}_\mathrm{MS+tug}\,\boldsymbol{\omega}|$   '
         r'(tug on arm at full extension, debris NOT coupled)',
         x=0.25)

    footer(ax, 10, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def page_tsiolkovsky(pdf):
    fig, ax = new_fig()
    heading(ax, 0.93, '10.  Propellant Budget & Tsiolkovsky Equation')
    hline(ax, 0.915)

    box(ax, 0.898, 0.665)
    heading(ax, 0.883, 'Tsiolkovsky rocket equation', fs=LABEL_FS)
    eq(ax, 0.848, r'$\Delta V$  =',
       r'$\Delta V = v_\mathrm{ex}\,\ln\!\dfrac{m_0}{m_0 - m_\mathrm{prop}}$',
       x_rhs=0.22)
    eq(ax, 0.795, r'$m_\mathrm{prop}$  =',
       r'$m_\mathrm{prop} = m_0\!\left(1 - e^{-\Delta V / v_\mathrm{ex}}\right)$',
       x_rhs=0.22)
    eq(ax, 0.742, r'$v_\mathrm{ex}$  =',
       r'$v_\mathrm{ex} = I_\mathrm{sp} \cdot g_0 = 253 \times 9.8066 = 2481\;\mathrm{m/s}$',
       x_rhs=0.22)
    param_block(ax, 0.698, [
        r'$m_0$              = chaser wet mass at the start of the phase  [kg]',
        r'                     (decreases through the mission as prop is burned)',
        r'$I_\mathrm{sp}$    = 253 s  (MS monoprop, confirmed)',
        r'$g_0$              = 9.8066 m/s²  (standard gravity)',
    ], x=0.25)

    box(ax, 0.628, 0.445)
    heading(ax, 0.613, 'KOS1 derivation (per-target, uplinked by GO at cut-off)', fs=LABEL_FS)
    eq(ax, 0.578, r'$r_\mathrm{KOS1}$  =',
       r'$r_\mathrm{KOS1} = \dfrac{s_\mathrm{panel}}{2} \times \kappa_\mathrm{margin}$',
       x_rhs=0.28)
    param_block(ax, 0.528, [
        r'$s_\mathrm{panel}$ = solar-panel tip-to-tip span, confirmed per target from GO uplink',
        r'$\kappa_\mathrm{margin}$ = 1.20  (+20% safety margin)',
        r'COMSATBW-1 (BSS-702, confirmed):  17.2/2 × 1.20 = 10.32 m',
    ], x=0.25)

    box(ax, 0.410, 0.260)
    heading(ax, 0.395, 'KOS2 derivation (fixed geometry)', fs=LABEL_FS)
    eq(ax, 0.360, r'$r_\mathrm{KOS2}$  =',
       r'$r_\mathrm{KOS2} = L_\mathrm{arm} + \dfrac{L_\mathrm{tug}}{2}$',
       x_rhs=0.28)
    param_block(ax, 0.310, [
        r'$L_\mathrm{arm}$ = 5.0 m  (arm length TBD MIRON)',
        r'$L_\mathrm{tug}$ = 0.8 m  (tug body length, confirmed)',
        r'Current value: 5.0 + 0.4 = 5.40 m  (auto-updates when MIRON provides arm length)',
    ], x=0.25)

    footer(ax, 11, 11)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

with PdfPages(SAVE_PATH) as pdf:
    page_title(pdf)
    page_mean_motion(pdf)
    page_cw(pdf)
    page_stm(pdf)
    page_targeting(pdf)
    page_perturbations(pdf)
    page_euler(pdf)
    page_quat(pdf)
    page_inertia(pdf)
    page_attitude_control(pdf)
    page_tsiolkovsky(pdf)

    d = pdf.infodict()
    d['Title']   = 'REAVER RPO Simulation — Equation Reference'
    d['Subject'] = 'RPO ΔV simulation equations for final report'

print(f"Saved -> {SAVE_PATH}")
