"""
REAVER Mission – Primary Structure Truss Analysis
==================================================
Octagonal CFRP truss, 3.5 m tall x 3.0 m wide
Falcon 9 launch environment

Analyses:
  1. Static load cases  (axial 6g, lateral 2g, combined)
  2. Modal analysis     (natural frequencies vs F9 min. requirements)
  3. Miles equation     (random vibration equivalent quasi-static loads)

Material: CFRP Quasi-Isotropic
  E   = 54  GPa
  rho = 1550 kg/m³
  su  = 600  MPa

Tube cross-section (parameterisable):
  OD  = 40 mm  (outer diameter)
  t   =  2 mm  (wall thickness)
"""

import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# 0.  PARAMETERS  (edit here)
# =============================================================================

# --- Material (CFRP Q-I) ---
E_Pa   = 54e9          # Young's modulus [Pa]
rho    = 1550.0        # Density [kg/m³]
su_Pa  = 600e6         # Ultimate tensile/compressive strength [Pa]
SF     = 1.5           # Safety factor (ESA ECSS typical)

# --- Tube cross-section ---
OD_m   = 0.040         # Outer diameter [m]
t_m    = 0.002         # Wall thickness [m]
ID_m   = OD_m - 2*t_m
A      = np.pi/4 * (OD_m**2 - ID_m**2)   # Cross-sectional area [m²]
I_sec  = np.pi/64 * (OD_m**4 - ID_m**4)  # Second moment of area [m⁴] (info only)

# --- Spacecraft ---
total_mass = 3800.0    # kg  (total S/C mass incl. payload)
g0         = 9.80665   # m/s²

# --- Falcon 9 quasi-static load factors (SPX payload user guide) ---
ax_g  =  6.0    # axial (thrust direction, worst-case liftoff/MECO)
lat_g =  2.0    # lateral

# --- Falcon 9 minimum frequency requirements ---
f_axial_req  = 25.0    # Hz  (axial)
f_lat_req    = 10.0    # Hz  (lateral)

# --- Random vibration (GRMS input) ---
GRMS       = 5.13      # g_rms  (mission requirement)
Q_factor   = 10.0      # amplification factor (typical for lightly damped structure)
# Miles eq.: a_peak = sqrt(pi/2 * f_n * Q * PSD)
# Equivalent: a_eq [g] = Q_factor * GRMS  (simplified bounding approach)
# More precisely per Miles: a_rms = sqrt(pi/2 * fn * Q * W0)
# We derive equivalent quasi-static g from modal freq below.

# --- Geometry ---
H      = 3.5           # total height [m]
R      = 1.5           # radius of octagon (centre to vertex) [m]
n_side = 8             # octagon
n_ring = 4             # number of horizontal rings (nodes per ring = n_side)

# =============================================================================
# 1.  GEOMETRY: BUILD NODE COORDINATES AND ELEMENT CONNECTIVITY
# =============================================================================

def octagon_ring(z, R, n=8):
    """Returns (n,3) array of node coords for one horizontal ring."""
    angles = np.array([2*np.pi*i/n for i in range(n)])
    x = R * np.cos(angles)
    y = R * np.sin(angles)
    return np.column_stack([x, y, np.full(n, z)])

# Ring heights: 0, H/3, 2H/3, H  (4 rings → 3 bays matching image)
z_rings = np.linspace(0, H, n_ring)

nodes = np.vstack([octagon_ring(z, R) for z in z_rings])
N = len(nodes)   # total nodes = 32

# Build connectivity: (element_start_node, element_end_node)
elements = []

def node_id(ring, vertex):
    return ring * n_side + (vertex % n_side)

for r in range(n_ring):
    # --- horizontal ring members ---
    for v in range(n_side):
        elements.append((node_id(r, v), node_id(r, v+1)))

for r in range(n_ring - 1):
    for v in range(n_side):
        # --- vertical members ---
        elements.append((node_id(r, v), node_id(r+1, v)))
        # --- diagonal bracing (alternating) ---
        if v % 2 == 0:
            elements.append((node_id(r, v),   node_id(r+1, v+1)))
        else:
            elements.append((node_id(r, v+1), node_id(r+1, v)))
        # --- cross-bracing inside each bay (radial diagonals matching image) ---
        elements.append((node_id(r, v), node_id(r+1, (v + n_side//2) % n_side)))

elements = list(set(tuple(sorted(e)) for e in elements))  # remove duplicates
n_el = len(elements)

# =============================================================================
# 2.  GLOBAL STIFFNESS AND MASS ASSEMBLY (pin-jointed truss)
# =============================================================================

DOF = 3
n_dof = N * DOF

K_global = np.zeros((n_dof, n_dof))
M_global = np.zeros((n_dof, n_dof))

elem_lengths = []

for (i, j) in elements:
    xi, yi, zi = nodes[i]
    xj, yj, zj = nodes[j]
    L = np.sqrt((xj-xi)**2 + (yj-yi)**2 + (zj-zi)**2)
    elem_lengths.append(L)

    cx = (xj-xi)/L;  cy = (yj-yi)/L;  cz = (zj-zi)/L
    T  = np.array([cx, cy, cz])

    # Local stiffness → global (axial only)
    k = (E_Pa * A / L)
    TT = np.outer(T, T)
    ke = k * np.block([[ TT, -TT],
                        [-TT,  TT]])

    # Consistent mass matrix (bar element)
    m_el = rho * A * L
    me_local = (m_el/6) * np.block([[2*np.eye(3), np.eye(3)],
                                     [np.eye(3),   2*np.eye(3)]])

    # Assemble
    dofs = [i*3, i*3+1, i*3+2, j*3, j*3+1, j*3+2]
    for a in range(6):
        for b in range(6):
            K_global[dofs[a], dofs[b]] += ke[a, b]
            M_global[dofs[a], dofs[b]] += me_local[a, b]

# Add non-structural mass (payload, systems) as lumped mass at top ring
ns_mass = total_mass - rho * A * np.sum(elem_lengths)
ns_mass = max(ns_mass, 0)
top_ring_nodes = [node_id(n_ring-1, v) for v in range(n_side)]
m_lump = ns_mass / len(top_ring_nodes)
for n_idx in top_ring_nodes:
    for d in range(3):
        M_global[n_idx*3+d, n_idx*3+d] += m_lump

# =============================================================================
# 3.  BOUNDARY CONDITIONS: fix bottom ring (all 3 DOF)
# =============================================================================

fixed_dofs = []
for v in range(n_side):
    n_idx = node_id(0, v)
    fixed_dofs += [n_idx*3, n_idx*3+1, n_idx*3+2]

free_dofs = [d for d in range(n_dof) if d not in fixed_dofs]

Kff = K_global[np.ix_(free_dofs, free_dofs)]
Mff = M_global[np.ix_(free_dofs, free_dofs)]

# =============================================================================
# 4.  STATIC ANALYSIS
# =============================================================================

def apply_load(ax_g_val, lat_g_val, mass_kg):
    """Build force vector: load distributed over free nodes (top 3 rings)."""
    F = np.zeros(n_dof)
    loaded = [node_id(r, v) for r in range(1, n_ring) for v in range(n_side)]
    m_per_node = mass_kg / len(loaded)
    for n_idx in loaded:
        F[n_idx*3+2] -= m_per_node * ax_g_val  * g0   # axial  (Z)
        F[n_idx*3+0] += m_per_node * lat_g_val * g0   # lateral (X)
    return F

def static_solve(F):
    Ff = F[free_dofs]
    U_free = np.linalg.solve(Kff, Ff)
    U = np.zeros(n_dof)
    for idx, d in enumerate(free_dofs):
        U[d] = U_free[idx]
    return U

def element_stress(U):
    """Axial stress in each element [Pa]."""
    stresses = []
    for (i, j) in elements:
        xi, yi, zi = nodes[i]; xj, yj, zj = nodes[j]
        L  = np.sqrt((xj-xi)**2 + (yj-yi)**2 + (zj-zi)**2)
        cx = (xj-xi)/L; cy = (yj-yi)/L; cz = (zj-zi)/L
        T  = np.array([cx, cy, cz])
        ui = U[i*3:i*3+3]; uj = U[j*3:j*3+3]
        delta = np.dot(T, uj - ui)
        sigma = E_Pa * delta / L
        stresses.append(sigma)
    return np.array(stresses)

load_cases = {
    'Axial (6g)':            (6.0, 0.0),
    'Lateral (2g)':          (0.0, 2.0),
    'Combined (6g + 2g)':    (6.0, 2.0),
}

print("=" * 65)
print("REAVER PRIMARY STRUCTURE – STATIC ANALYSIS")
print("=" * 65)
print(f"  Tube: OD={OD_m*1e3:.0f} mm, t={t_m*1e3:.0f} mm,  A={A*1e6:.2f} cm²")
print(f"  n_nodes={N},  n_elements={n_el}")
print(f"  Structural mass: {rho*A*np.sum(elem_lengths):.1f} kg")
print(f"  Non-structural lumped mass: {ns_mass:.1f} kg")
print(f"  Safety factor: {SF}")
print()

static_results = {}
for name, (ag, lg) in load_cases.items():
    F = apply_load(ag, lg, total_mass)
    U = static_solve(F)
    sigma = element_stress(U)

    max_tensile    = np.max(sigma)
    max_compressive= np.min(sigma)
    max_abs        = np.max(np.abs(sigma))
    margin         = su_Pa / (SF * max_abs) - 1.0

    max_disp = np.max(np.abs(U[free_dofs]))

    static_results[name] = {
        'U': U, 'sigma': sigma,
        'max_t': max_tensile, 'max_c': max_compressive,
        'max_abs': max_abs, 'margin': margin, 'max_disp': max_disp
    }

    status = "PASS ✓" if margin > 0 else "FAIL ✗"
    print(f"  Load case: {name}")
    print(f"    Max tensile stress    : {max_tensile/1e6:+.1f} MPa")
    print(f"    Max compressive stress: {max_compressive/1e6:+.1f} MPa")
    print(f"    Max nodal displacement: {max_disp*1e3:.2f} mm")
    print(f"    Margin of Safety (MoS): {margin:+.3f}  [{status}]")
    print()

# =============================================================================
# 5.  MODAL ANALYSIS
# =============================================================================

print("=" * 65)
print("MODAL ANALYSIS")
print("=" * 65)

n_modes = 10
eigenvalues, eigenvectors = eigh(Kff, Mff, subset_by_index=[0, n_modes-1])
eigenvalues = np.maximum(eigenvalues, 0)
freqs_hz = np.sqrt(eigenvalues) / (2 * np.pi)

# Identify mode character (dominant direction of motion)
def mode_character(evec, free_dofs):
    U_full = np.zeros(n_dof)
    for idx, d in enumerate(free_dofs):
        U_full[d] = evec[idx]
    ux = np.sum(U_full[0::3]**2)
    uy = np.sum(U_full[1::3]**2)
    uz = np.sum(U_full[2::3]**2)
    dominant = np.argmax([ux, uy, uz])
    labels = ['Lateral-X', 'Lateral-Y', 'Axial-Z']
    return labels[dominant]

print(f"  {'Mode':<6} {'Freq [Hz]':>10} {'Character':<14} {'F9 Req [Hz]':>12} {'Status'}")
print(f"  {'-'*6} {'-'*10} {'-'*14} {'-'*12} {'-'*10}")

modal_results = []
for m in range(n_modes):
    f  = freqs_hz[m]
    ev = eigenvectors[:, m]
    ch = mode_character(ev, free_dofs)
    req = f_axial_req if 'Axial' in ch else f_lat_req
    ok  = "PASS ✓" if f >= req else "FAIL ✗"
    modal_results.append({'f': f, 'char': ch, 'req': req, 'pass': f >= req})
    print(f"  {m+1:<6} {f:>10.2f} {ch:<14} {req:>12.1f} {ok}")

f_axial_actual = min((r['f'] for r in modal_results if 'Axial' in r['char']), default=None)
f_lat_actual   = min((r['f'] for r in modal_results if 'Lateral' in r['char']), default=None)
print()

# =============================================================================
# 6.  MILES EQUATION – RANDOM VIBRATION EQUIVALENT LOADS
# =============================================================================

print("=" * 65)
print("RANDOM VIBRATION – MILES EQUATION")
print("=" * 65)
print(f"  Input GRMS = {GRMS:.2f} g,  Q = {Q_factor:.0f}")
print()

# PSD level assumed flat (0dB/oct) — conservative bounding
# GRMS² = PSD * BW  →  PSD = GRMS² / BW
# For Miles: a_rms = sqrt(pi/2 * fn * Q * PSD)
# We use BW = fn (conservative single-peak assumption)

def miles_g(fn_hz, grms, Q):
    """Peak equivalent g from Miles equation (3-sigma = 3*a_rms)."""
    PSD = grms**2 / fn_hz     # flat PSD assumption [g²/Hz]
    a_rms = np.sqrt(np.pi/2 * fn_hz * Q * PSD)
    return 3.0 * a_rms        # 3-sigma peak

print(f"  {'Mode':<6} {'Freq [Hz]':>10} {'Character':<14} {'Miles g_eq':>12} {'MoS':>8} {'Status'}")
print(f"  {'-'*6} {'-'*10} {'-'*14} {'-'*12} {'-'*8} {'-'*8}")

for m, r in enumerate(modal_results):
    fn = r['f']
    if fn < 0.1:
        continue
    ch = r['char']
    g_eq = miles_g(fn, GRMS, Q_factor)

    ax_l = g_eq if 'Axial' in ch else 0.0
    la_l = g_eq if 'Lateral' in ch else 0.0
    F    = apply_load(ax_l, la_l, total_mass)
    U    = static_solve(F)
    sig  = element_stress(U)
    mos  = su_Pa / (SF * np.max(np.abs(sig))) - 1.0
    ok   = "PASS ✓" if mos > 0 else "FAIL ✗"
    print(f"  {m+1:<6} {fn:>10.2f} {ch:<14} {g_eq:>12.2f} {mos:>8.3f} {ok}")

print()

# =============================================================================
# 7.  SUMMARY
# =============================================================================

print("=" * 65)
print("SUMMARY")
print("=" * 65)

all_mos   = [r['margin'] for r in static_results.values()]
min_mos   = min(all_mos)
print(f"  Min static MoS (all load cases): {min_mos:+.3f}  "
      f"{'PASS ✓' if min_mos > 0 else 'FAIL ✗'}")

if f_axial_actual:
    m_ax = f_axial_actual / f_axial_req - 1
    print(f"  Axial freq:   {f_axial_actual:.2f} Hz  (req ≥ {f_axial_req} Hz)  "
          f"margin {m_ax:+.2f}  {'PASS ✓' if m_ax > 0 else 'FAIL ✗'}")
if f_lat_actual:
    m_la = f_lat_actual / f_lat_req - 1
    print(f"  Lateral freq: {f_lat_actual:.2f} Hz  (req ≥ {f_lat_req} Hz)  "
          f"margin {m_la:+.2f}  {'PASS ✓' if m_la > 0 else 'FAIL ✗'}")

print()

# =============================================================================
# 8.  PLOTS
# =============================================================================

fig = plt.figure(figsize=(18, 12))
fig.patch.set_facecolor('#1a1a2e')

# --- 8a. 3-D truss geometry ---
ax1 = fig.add_subplot(2, 3, 1, projection='3d')
ax1.set_facecolor('#1a1a2e')
for (i, j) in elements:
    xs = [nodes[i,0], nodes[j,0]]
    ys = [nodes[i,1], nodes[j,1]]
    zs = [nodes[i,2], nodes[j,2]]
    ax1.plot(xs, ys, zs, 'c-', lw=0.8, alpha=0.7)
ax1.scatter(nodes[:,0], nodes[:,1], nodes[:,2], c='white', s=20, zorder=5)
ax1.set_title('Truss Geometry', color='white', fontsize=10)
ax1.tick_params(colors='white', labelsize=7)
for pane in [ax1.xaxis.pane, ax1.yaxis.pane, ax1.zaxis.pane]:
    pane.fill = False
ax1.set_xlabel('X [m]', color='white', fontsize=7)
ax1.set_ylabel('Y [m]', color='white', fontsize=7)
ax1.set_zlabel('Z [m]', color='white', fontsize=7)

# --- 8b. Element stresses – combined load case ---
ax2 = fig.add_subplot(2, 3, 2)
ax2.set_facecolor('#1a1a2e')
sig_comb = static_results['Combined (6g + 2g)']['sigma'] / 1e6
colors = ['#ff4444' if s > 0 else '#4488ff' for s in sig_comb]
ax2.bar(range(n_el), sig_comb, color=colors, width=1.0, alpha=0.8)
ax2.axhline(su_Pa/1e6/SF,  color='yellow', lw=1.5, ls='--', label=f'+σ_allow={su_Pa/SF/1e6:.0f} MPa')
ax2.axhline(-su_Pa/1e6/SF, color='yellow', lw=1.5, ls='--', label=f'-σ_allow')
ax2.set_title('Element Stresses – Combined (6g+2g)', color='white', fontsize=10)
ax2.set_xlabel('Element index', color='white', fontsize=8)
ax2.set_ylabel('Stress [MPa]', color='white', fontsize=8)
ax2.tick_params(colors='white')
ax2.legend(fontsize=7, facecolor='#2a2a4e', labelcolor='white')
ax2.spines[:].set_color('#444')

# --- 8c. Natural frequencies vs requirements ---
ax3 = fig.add_subplot(2, 3, 3)
ax3.set_facecolor('#1a1a2e')
mode_nums = np.arange(1, n_modes+1)
bar_colors = ['#ff4444' if not r['pass'] else '#44ff88' for r in modal_results]
ax3.bar(mode_nums, freqs_hz, color=bar_colors, alpha=0.85)
ax3.axhline(f_axial_req, color='orange', lw=1.5, ls='--', label=f'F9 axial req {f_axial_req} Hz')
ax3.axhline(f_lat_req,   color='cyan',   lw=1.5, ls='--', label=f'F9 lateral req {f_lat_req} Hz')
ax3.set_title('Natural Frequencies (first 10 modes)', color='white', fontsize=10)
ax3.set_xlabel('Mode number', color='white', fontsize=8)
ax3.set_ylabel('Frequency [Hz]', color='white', fontsize=8)
ax3.tick_params(colors='white')
ax3.legend(fontsize=7, facecolor='#2a2a4e', labelcolor='white')
ax3.spines[:].set_color('#444')

# --- 8d. Displaced shape – axial load ---
ax4 = fig.add_subplot(2, 3, 4, projection='3d')
ax4.set_facecolor('#1a1a2e')
U_ax = static_results['Axial (6g)']['U']
scale = 50
disp_nodes = nodes + scale * U_ax.reshape(-1, 3)
for (i, j) in elements:
    ax4.plot([nodes[i,0], nodes[j,0]], [nodes[i,1], nodes[j,1]], [nodes[i,2], nodes[j,2]],
             'c-', lw=0.5, alpha=0.3)
    ax4.plot([disp_nodes[i,0], disp_nodes[j,0]],
             [disp_nodes[i,1], disp_nodes[j,1]],
             [disp_nodes[i,2], disp_nodes[j,2]], 'r-', lw=1.2, alpha=0.9)
ax4.set_title(f'Displaced Shape – Axial 6g (×{scale})', color='white', fontsize=10)
ax4.tick_params(colors='white', labelsize=7)
for pane in [ax4.xaxis.pane, ax4.yaxis.pane, ax4.zaxis.pane]:
    pane.fill = False
ax4.set_xlabel('X [m]', color='white', fontsize=7)
ax4.set_ylabel('Y [m]', color='white', fontsize=7)
ax4.set_zlabel('Z [m]', color='white', fontsize=7)

# --- 8e. Displaced shape – lateral load ---
ax5 = fig.add_subplot(2, 3, 5, projection='3d')
ax5.set_facecolor('#1a1a2e')
U_la = static_results['Lateral (2g)']['U']
disp_nodes_la = nodes + scale * U_la.reshape(-1, 3)
for (i, j) in elements:
    ax5.plot([nodes[i,0], nodes[j,0]], [nodes[i,1], nodes[j,1]], [nodes[i,2], nodes[j,2]],
             'c-', lw=0.5, alpha=0.3)
    ax5.plot([disp_nodes_la[i,0], disp_nodes_la[j,0]],
             [disp_nodes_la[i,1], disp_nodes_la[j,1]],
             [disp_nodes_la[i,2], disp_nodes_la[j,2]], '#ff88aa', lw=1.2, alpha=0.9)
ax5.set_title(f'Displaced Shape – Lateral 2g (×{scale})', color='white', fontsize=10)
ax5.tick_params(colors='white', labelsize=7)
for pane in [ax5.xaxis.pane, ax5.yaxis.pane, ax5.zaxis.pane]:
    pane.fill = False
ax5.set_xlabel('X [m]', color='white', fontsize=7)
ax5.set_ylabel('Y [m]', color='white', fontsize=7)
ax5.set_zlabel('Z [m]', color='white', fontsize=7)

# --- 8f. MoS summary bar chart ---
ax6 = fig.add_subplot(2, 3, 6)
ax6.set_facecolor('#1a1a2e')
names  = list(static_results.keys())
mos_vals = [static_results[n]['margin'] for n in names]
bar_c  = ['#44ff88' if m > 0 else '#ff4444' for m in mos_vals]
bars   = ax6.barh(names, mos_vals, color=bar_c, alpha=0.85)
ax6.axvline(0, color='white', lw=1)
ax6.set_title('Margin of Safety – Static Cases', color='white', fontsize=10)
ax6.set_xlabel('MoS  (>0 = PASS)', color='white', fontsize=8)
ax6.tick_params(colors='white')
ax6.spines[:].set_color('#444')
for bar, v in zip(bars, mos_vals):
    ax6.text(v + 0.02, bar.get_y() + bar.get_height()/2,
             f'{v:+.3f}', va='center', color='white', fontsize=8)

plt.suptitle('REAVER – Primary Structure Analysis  |  CFRP Q-I, OD=40mm, t=2mm',
             color='white', fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/reaver_structural_analysis.png',
            dpi=150, bbox_inches='tight', facecolor='#1a1a2e')
plt.close()
print("  Plot saved → reaver_structural_analysis.png")
print()
print("  To change tube sizing, edit OD_m and t_m at the top of the script.")
print("  To change load factors, edit ax_g and lat_g.")
print("=" * 65)
