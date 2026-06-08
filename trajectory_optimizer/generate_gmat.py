#!/usr/bin/env python3
"""
generate_gmat.py
----------------
Generates a GMAT script for a combined plane-change + altitude maneuver
between two circular orbits, using either:
  - Impulsive burns (ImpulsiveBurn)
  - Chemical finite burns (ChemicalThruster + FiniteBurn, ~minutes each)
  - Electric low-thrust spiral (ElectricThruster + single continuous FiniteBurn)

Usage:
    python generate_gmat.py

Edit the USER INPUTS section, then run. The script writes OUTPUT to disk.

Architecture per propulsion type:

  IMPULSIVE / CHEMICAL:
    Two-burn Hohmann transfer (or single burn if SMA unchanged).
    Burn 1 at u1 on the initial orbit: raises apogee + begins plane rotation.
    Burn 2 at the antipodal point on the transfer ellipse: circularises +
    completes plane rotation. DC varies burn components (impulsive) or
    durations (chemical finite) to meet SMA, INC, RAAN goals.

  ELECTRIC (low-thrust spiral):
    A single continuous burn throughout. The thrust direction is fixed in VNB
    as the analytically computed unit vector of (dVV, dVN) -- the same vector
    that a Hohmann burn would use, averaged across both burns and weighted by
    DeltaV. The spacecraft spirals slowly from the initial to the final orbit.
    DC varies the single burn duration to meet INC. SMA and RAAN are
    achieved by the analytically set thrust direction and the orbit evolution.
    No coast phase -- the Hohmann geometry is invalid for low-thrust.

  If SMA_init == SMA_final the transfer collapses to a pure plane change:
    Impulsive/Chemical: single burn (open-loop finite or DC-targeted impulsive).
    Electric: single continuous burn, duration from rocket equation.
"""

import math, os, sys

# ===========================================================
# USER INPUTS -- orbit
# ===========================================================
SMA_init   = 42164.2  # km   initial semi-major axis
INC_init   =     0.06960000000589468    # deg  initial inclination
RAAN_init  =    89.7    # deg  initial RAAN

SMA_final = 42875.0
# SMA_final  = 42164.169    # km   final semi-major axis
INC_final  =    7    # deg  final inclination
RAAN_final =    64    # deg  final RAAN

# ===========================================================
# USER INPUTS -- propulsion
# ===========================================================
# Propulsion type: 'impulsive', 'chemical', or 'electric'
PROP_TYPE  = 'chemical'

# --- Spacecraft ---
DRY_MASS   = 2440.0 + 200    # kg
FUEL_MASS  =  100.0    # kg   (initial propellant mass)

# --- Chemical thruster (used when PROP_TYPE = 'chemical') ---
CHEM_THRUST = 88    # N    constant thrust level
CHEM_ISP    = 300.0    # s    specific impulse

# --- Electric thruster (used when PROP_TYPE = 'electric') ---
ELEC_THRUST = 0.075      # N    constant thrust level
ELEC_ISP    = 1600.0   # s    specific impulse
ELEC_POWER  = 2      # kW   available power to thruster

# ===========================================================
# Other settings
# ===========================================================
EPOCH  = '01 Jan 2025 00:00:00.000'
OUTPUT = 'GEO_RAAN_Change.script'
# ===========================================================

MU   = 398600.4418   # km^3/s^2
G0   = 9.80665       # m/s^2
G0km = G0 * 1e-3     # km/s^2

# -----------------------------------------------------------
# Select active propulsion parameters
# -----------------------------------------------------------
def active_isp():
    if PROP_TYPE == 'electric':
        return ELEC_ISP
    elif PROP_TYPE == 'chemical':
        return CHEM_ISP
    else:
        return CHEM_ISP   # impulsive uses same formula for fuel estimate

def active_thrust_N():
    if PROP_TYPE == 'electric':
        return ELEC_THRUST
    return CHEM_THRUST

# -----------------------------------------------------------
# Vector helpers
# -----------------------------------------------------------
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def dot(a,b):   return sum(x*y for x,y in zip(a,b))
def norm(a):    return math.sqrt(dot(a,a))
def normalize(a): m=norm(a); return [x/m for x in a]
def scale(a,s): return [x*s for x in a]
def add(a,b):   return [x+y for x,y in zip(a,b)]

# -----------------------------------------------------------
# Propulsion calculations
# -----------------------------------------------------------
def fuel_kg(dv_kms, m0_kg, isp_s):
    """Propellant mass from Tsiolkovsky equation."""
    return m0_kg * (1.0 - math.exp(-abs(dv_kms) / (isp_s * G0km)))

def burn_duration_s(dv_kms, m0_kg, thrust_N, isp_s):
    """
    Burn duration estimate from rocket equation + constant thrust.
    mdot = F / (Isp * g0)
    dm   = m0 * (1 - exp(-dV / (Isp*g0)))
    t    = dm / mdot
    """
    dv_ms  = abs(dv_kms) * 1000.0           # convert to m/s
    mdot   = thrust_N / (isp_s * G0)        # kg/s
    dm     = m0_kg * (1.0 - math.exp(-dv_ms / (isp_s * G0)))
    return dm / mdot                          # seconds

def thrust_direction(dVV, dVN):
    """Unit vector of (dVV, dVN) burn in VNB for ThrustDirection fields."""
    mag = math.sqrt(dVV**2 + dVN**2)
    if mag < 1e-12:
        return 0.0, 1.0
    return dVV/mag, dVN/mag

def electric_params(SMA1, INC1_deg, RAAN1_deg, SMA2, INC2_deg, RAAN2_deg,
                    thrust_N, isp_s, m0_kg):
    """
    Compute all parameters for a two-phase electric low-thrust maneuver.

    Phase 1 -- altitude raise: continuous in-track (V) burn.
               ThrustDirection: V=+/-1, N=0.
               DC varies durV -> Achieve SC.Earth.SMA = SMA2.

    Phase 2 -- plane change: dual half-orbit out-of-plane burns, no coasting.
               +N over the half-orbit centred on the node line (ThrusterN),
               then -N over the antipodal half-orbit (ThrusterAntiN). The burn
               sign flips together with cos(u), so BOTH halves add the same
               inclination -- this halves the maneuver time versus a single
               burn-arc-per-orbit scheme and wastes no time coasting. The For
               loop repeats this N_ARCS times.

               Closed-loop targeting (set up in the mission sequence):
                 Vary N_ARCS   (orbit count)   -> Achieve SC.INC
                 Vary ta_start (arc-centre TA) -> Achieve SC.RAAN
               N_ARCS sets the plane-change MAGNITUDE; ta_start sets its
               DIRECTION (which inertial longitude the new node lands on).
               These are nearly orthogonal, so the 2x2 Jacobian is well
               conditioned.

    Key physics: a continuous fixed N-direction burn nets ZERO inclination
    change per orbit (cos(u) averages to zero); flipping the sign each half-
    orbit rectifies cos(u) so the inclination accumulates monotonically.
    """
    MU_loc = 398600.4418

    i1 = math.radians(INC1_deg); O1 = math.radians(RAAN1_deg)
    i2 = math.radians(INC2_deg); O2 = math.radians(RAAN2_deg)

    # --- Phase 1: altitude ---
    a_t  = (SMA1 + SMA2) / 2.0
    v1c  = math.sqrt(MU_loc / SMA1)
    v2c  = math.sqrt(MU_loc / SMA2)
    vtp  = math.sqrt(MU_loc * (2.0/SMA1 - 1.0/a_t))
    vta  = math.sqrt(MU_loc * (2.0/SMA2 - 1.0/a_t))
    dv_alt = abs(vtp - v1c) + abs(v2c - vta)          # km/s
    dm1    = fuel_kg(dv_alt, m0_kg, isp_s)
    t1_dur = burn_duration_s(dv_alt, m0_kg, thrust_N, isp_s)
    m1     = m0_kg - dm1
    # sign of V thrust: +1 to raise orbit, -1 to lower
    td1_v  = 1.0 if SMA2 >= SMA1 else -1.0

    # --- Optimal burn arc centre for plane change ---
    # Centre the arcs on the line of nodes BETWEEN the initial and final
    # planes (axis = h1 x h2): the axis about which h must rotate to go
    # h1 -> h2. Crossing it out-of-plane is the most efficient push AND
    # lands the resulting node at the correct RAAN. This is robust for
    # near-equatorial starts, where the old GVE linearisation
    #     tan(u) = (dRAAN/dINC)*sin(i1)
    # collapses to ~0 (sin(i1)->0) and wrongly parks the node at the
    # STARTING longitude instead of moving it to the target RAAN.
    #
    # NOTE: this only sets the INITIAL GUESS for the arc centre. The mission
    # sequence closes the loop on RAAN by letting the DC vary ta_start.
    dINC = INC2_deg - INC1_deg
    def _hvec(inc_deg, raan_deg):
        ii = math.radians(inc_deg); oo = math.radians(raan_deg)
        return [math.sin(ii)*math.sin(oo), -math.sin(ii)*math.cos(oo), math.cos(ii)]
    h1 = _hvec(INC1_deg, RAAN1_deg)
    h2 = _hvec(INC2_deg, RAAN2_deg)
    axis = cross(h1, h2)
    if norm(axis) < 1e-12:                       # planes identical / antiparallel
        u_opt = 90.0 if (RAAN2_deg - RAAN1_deg) >= 0 else 270.0
        td2_n = 1.0 if dINC >= 0 else -1.0
    else:
        axis = normalize(axis)
        # basis of the initial orbit plane: P toward asc. node, Q 90 deg ahead
        P = [math.cos(O1), math.sin(O1), 0.0]
        Q = [-math.cos(i1)*math.sin(O1), math.cos(i1)*math.cos(O1), math.sin(i1)]
        u_opt = math.degrees(math.atan2(dot(axis, Q), dot(axis, P))) % 360.0
        # choose the node crossing where a +N push rotates h1 toward h2
        target = [h2[k] - h1[k] for k in range(3)]
        r0, v0 = state_on_circle(SMA1, i1, O1, math.radians(u_opt))
        if dot(cross(r0, normalize(cross(r0, v0))), target) < 0:
            u_opt = (u_opt + 180.0) % 360.0
        td2_n = 1.0   # +N over the arc centred on u_opt rotates h1 -> h2

    # --- Phase 2: plane change ---
    cos_a  = (math.cos(i1)*math.cos(i2) +
               math.sin(i1)*math.sin(i2)*math.cos(O2 - O1))
    alpha  = math.acos(max(-1.0, min(1.0, cos_a)))
    dv_plane = 2.0 * v2c * math.sin(alpha / 2.0)      # km/s

    dm2    = fuel_kg(dv_plane, m1, isp_s)
    t2_burn = burn_duration_s(dv_plane, m1, thrust_N, isp_s)  # total on-thrust time

    T2          = 2.0 * math.pi * math.sqrt(SMA2**3 / MU_loc)
    ARC_DEG     = 180.0      # deg burned per half-orbit (must tile to 360 with the anti-arc)
    arc_frac    = ARC_DEG / 360.0
    t_pos_arc   = arc_frac * T2                     # +N burn time per orbit
    t_neg_arc   = (1.0 - arc_frac) * T2             # -N (anti) burn time per orbit
    # Initial-guess orbit count. Note: with the anti-burn BOTH halves add
    # inclination, so the real per-orbit plane change is ~2x the continuous-
    # equivalent baked into t2_burn; the DC corrects N_ARCS from this seed.
    n_arcs      = math.ceil(t2_burn / (t_pos_arc + t_neg_arc))
    t2_elapsed  = n_arcs * T2                        # total elapsed time

    # TA window for arc: [u_opt - ARC_DEG/2, u_opt + ARC_DEG/2]
    ta_start = (u_opt - ARC_DEG / 2.0) % 360.0
    ta_end   = (u_opt + ARC_DEG / 2.0) % 360.0

    return dict(
        dv_alt=dv_alt, dm1=dm1, t1_dur=t1_dur, m1=m1, td1_v=td1_v,
        dv_plane=dv_plane, dm2=dm2, t2_burn=t2_burn,
        u_opt=u_opt, td2_n=td2_n,
        T2=T2, ARC_DEG=ARC_DEG, t_pos_arc=t_pos_arc, t_neg_arc=t_neg_arc,
        n_arcs=n_arcs, t2_elapsed=t2_elapsed,
        ta_start=ta_start, ta_end=ta_end,
        fuel_total=dm1+dm2,
        same_alt=(abs(SMA2-SMA1)<0.1),
    )

# -----------------------------------------------------------
# Orbital mechanics helpers
# -----------------------------------------------------------
def state_on_circle(SMA, INC_r, RAAN_r, u_r):
    cO,sO = math.cos(RAAN_r), math.sin(RAAN_r)
    ci,si = math.cos(INC_r),  math.sin(INC_r)
    def rot(v):
        return [cO*v[0]-sO*ci*v[1]+sO*si*v[2],
                sO*v[0]+cO*ci*v[1]-cO*si*v[2],
                si*v[1]+ci*v[2]]
    vc = math.sqrt(MU/SMA)
    r  = scale(rot([math.cos(u_r),  math.sin(u_r), 0.0]), SMA)
    v  = scale(rot([-math.sin(u_r), math.cos(u_r), 0.0]), vc)
    return r, v

def vnb_burn(r, v, dVV, dVN):
    V = normalize(v)
    N = normalize(cross(r, v))
    dv = add(scale(V, dVV), scale(N, dVN))
    return add(v, dv)

def rv_to_inc_raan_sma(r, v):
    rmag = norm(r); vmag = norm(v)
    h    = cross(r, v); hmag = norm(h)
    eps  = vmag**2/2.0 - MU/rmag
    sma  = -MU/(2.0*eps)
    inc  = math.degrees(math.acos(max(-1.0, min(1.0, h[2]/hmag))))
    n    = cross([0,0,1], h)
    raan = math.degrees(math.atan2(n[1], n[0])) % 360.0
    return inc, raan, sma

def propagate_to_antipode(r0, v0):
    h_vec = cross(r0, v0); h_hat = normalize(h_vec)
    rmag  = norm(r0); vmag = norm(v0)
    eps   = vmag**2/2.0 - MU/rmag
    sma   = -MU/(2.0*eps)
    e_vec = scale([v0[1]*h_vec[2]-v0[2]*h_vec[1] - MU*r0[0]/rmag,
                   v0[2]*h_vec[0]-v0[0]*h_vec[2] - MU*r0[1]/rmag,
                   v0[0]*h_vec[1]-v0[1]*h_vec[0] - MU*r0[2]/rmag], 1.0/MU)
    ecc = norm(e_vec)
    if ecc < 1e-9:
        r_a = scale(normalize(r0), -sma)
        v_a = scale(normalize(cross(h_hat, normalize(r_a))), math.sqrt(MU/sma))
        return r_a, v_a
    e_hat  = normalize(e_vec)
    r_a    = scale(scale(e_hat, -1.0), sma*(1.0+ecc))
    v_a    = scale(normalize(cross(h_hat, normalize(r_a))), norm(h_vec)/norm(r_a))
    return r_a, v_a

# -----------------------------------------------------------
# Analytical burn solver
# -----------------------------------------------------------
def plane_angle(i1, O1, i2, O2):
    c = math.cos(i1)*math.cos(i2) + math.sin(i1)*math.sin(i2)*math.cos(O2-O1)
    return math.acos(max(-1.0, min(1.0, c)))

def optimal_split(SMA1, SMA2, alpha_total):
    a_t   = (SMA1+SMA2)/2.0
    v1    = math.sqrt(MU/SMA1)
    v2    = math.sqrt(MU/SMA2)
    vtp   = math.sqrt(MU*(2.0/SMA1 - 1.0/a_t))
    vta   = math.sqrt(MU*(2.0/SMA2 - 1.0/a_t))
    def tdv(a1):
        a2 = alpha_total-a1
        return (math.sqrt(vtp**2+v1**2-2*vtp*v1*math.cos(a1)) +
                math.sqrt(v2**2+vta**2-2*v2*vta*math.cos(a2)))
    lo,hi = 0.0, alpha_total; gr = (math.sqrt(5)+1)/2
    for _ in range(300):
        c=hi-(hi-lo)/gr; d=lo+(hi-lo)/gr
        if tdv(c)<tdv(d): hi=d
        else: lo=c
        if abs(hi-lo)<1e-14: break
    a1=(lo+hi)/2; a2=alpha_total-a1
    return a1, a2, vtp, vta

def simulate(SMA1, INC1r, RAAN1r, u1_deg, dVV1, dVN1, dVV2, dVN2):
    r1,v1 = state_on_circle(SMA1, INC1r, RAAN1r, math.radians(u1_deg))
    v1b   = vnb_burn(r1, v1, dVV1, dVN1)
    r2,v2 = propagate_to_antipode(r1, v1b)
    v2b   = vnb_burn(r2, v2, dVV2, dVN2)
    return rv_to_inc_raan_sma(r2, v2b)

def solve_two_burn(SMA1, INC1_deg, RAAN1_deg, SMA2, INC2_deg, RAAN2_deg,
                   dVV1_i, dVN1_i, dVV2_i, dVN2_i):
    i1 = math.radians(INC1_deg); O1 = math.radians(RAAN1_deg)
    best_u, best_e = 0.0, 1e9
    for u_try in range(0, 360, 10):
        ic,rc,sc = simulate(SMA1,i1,O1,float(u_try),dVV1_i,dVN1_i,dVV2_i,dVN2_i)
        e = abs((RAAN2_deg-rc+180)%360-180) + abs(ic-INC2_deg)
        if e < best_e: best_e=e; best_u=float(u_try)
    u1=best_u; vv1=dVV1_i; vn1=dVN1_i; vv2=dVV2_i; vn2=dVN2_i
    DU=1e-5; DN=1e-7; DV=1e-7
    for _ in range(500):
        i0,r0,s0 = simulate(SMA1,i1,O1,u1,vv1,vn1,vv2,vn2)
        ei=INC2_deg-i0; er=(RAAN2_deg-r0+180)%360-180; es=SMA2-s0
        if abs(ei)<1e-8 and abs(er)<1e-8 and abs(es)<1e-5: break
        iu,ru,su = simulate(SMA1,i1,O1,u1+DU,vv1,vn1,vv2,vn2)
        in_,rn,sn= simulate(SMA1,i1,O1,u1,vv1,vn1+DN,vv2,vn2)
        iv,rv,sv = simulate(SMA1,i1,O1,u1,vv1,vn1,vv2+DV,vn2)
        J = [[(iu-i0)/DU,  (in_-i0)/DN, (iv-i0)/DV],
             [((ru-r0+180)%360-180)/DU, ((rn-r0+180)%360-180)/DN, ((rv-r0+180)%360-180)/DV],
             [(su-s0)/DU,  (sn-s0)/DN,  (sv-s0)/DV]]
        b=[ei,er,es]
        M=[J[i][:]+[b[i]] for i in range(3)]
        for col in range(3):
            mx=max(range(col,3),key=lambda r:abs(M[r][col]))
            M[col],M[mx]=M[mx],M[col]
            if abs(M[col][col])<1e-20: break
            for row in range(3):
                if row!=col:
                    f=M[row][col]/M[col][col]
                    M[row]=[M[row][k]-f*M[col][k] for k in range(4)]
        dx=[M[i][3]/M[i][i] if abs(M[i][i])>1e-20 else 0 for i in range(3)]
        lim=[3.0,0.05,0.05]
        dx=[max(-lim[k],min(lim[k],dx[k])) for k in range(3)]
        u1=(u1+dx[0])%360; vn1+=dx[1]; vv2+=dx[2]
    i_f,r_f,s_f=simulate(SMA1,i1,O1,u1,vv1,vn1,vv2,vn2)
    ok=(abs(i_f-INC2_deg)<0.01 and abs((RAAN2_deg-r_f+180)%360-180)<0.01 and abs(s_f-SMA2)<1.0)
    return u1,vv1,vn1,vv2,vn2,ok

def solve_single_burn(SMA, INC1_deg, RAAN1_deg, INC2_deg, RAAN2_deg):
    i1=math.radians(INC1_deg); O1=math.radians(RAAN1_deg)
    sign=1.0 if INC2_deg>=INC1_deg else -1.0
    cos_a=(math.cos(i1)*math.cos(math.radians(INC2_deg))+
           math.sin(i1)*math.sin(math.radians(INC2_deg))*
           math.cos(math.radians(RAAN2_deg-RAAN1_deg)))
    dv0=2*math.sqrt(MU/SMA)*math.sin(math.acos(max(-1,min(1,cos_a)))/2)
    dRAAN=RAAN2_deg-RAAN1_deg; dINC=INC2_deg-INC1_deg
    if abs(dINC)>1e-6:
        u0=math.degrees(math.atan((dRAAN/dINC)*math.sin(i1)))%360
    else:
        u0=90.0 if dRAAN>=0 else 270.0
    def apply_N(u_r, dv_N):
        r,v=state_on_circle(SMA,i1,O1,u_r)
        N=normalize(cross(r,v))
        v2=add(v,scale(N,dv_N))
        return rv_to_inc_raan_sma(r,v2)
    u_deg=u0; dv_N=sign*dv0
    for _ in range(500):
        i0,r0,_=apply_N(math.radians(u_deg),dv_N)
        ei=INC2_deg-i0; er=(RAAN2_deg-r0+180)%360-180
        if abs(ei)<1e-8 and abs(er)<1e-8: break
        DU=1e-5; DN=1e-7
        iu,ru,_=apply_N(math.radians(u_deg)+DU,dv_N)
        in_,rn,_=apply_N(math.radians(u_deg),dv_N+DN)
        J=[[( iu-i0)/DU,(in_-i0)/DN],
           [((ru-r0+180)%360-180)/DU,((rn-r0+180)%360-180)/DN]]
        det=J[0][0]*J[1][1]-J[0][1]*J[1][0]
        if abs(det)<1e-20: break
        du=(J[1][1]*ei-J[0][1]*er)/det; dn=(-J[1][0]*ei+J[0][0]*er)/det
        if abs(du)>3: du*=3/abs(du)
        if abs(dn)>0.05: dn*=0.05/abs(dn)
        u_deg=(u_deg+math.degrees(du))%360; dv_N+=dn
    i_f,r_f,_=apply_N(math.radians(u_deg),dv_N)
    ok=abs(i_f-INC2_deg)<0.01 and abs((RAAN2_deg-r_f+180)%360-180)<0.01
    return u_deg, dv_N, ok

# -----------------------------------------------------------
# GMAT script writer
# -----------------------------------------------------------
def write_script(path, p):
    L=[]; A=L.append
    def S(s=''): A(s)

    isp    = active_isp()
    thrust = active_thrust_N()
    m0     = DRY_MASS + FUEL_MASS
    finite = PROP_TYPE in ('chemical','electric')

    # ---- header ----
    S('%----------------------------------------------------------')
    S('%  GMAT Orbital Maneuver Script')
    S('%  Auto-generated by generate_gmat.py')
    S(f'%  Propulsion type : {PROP_TYPE.upper()}')
    S('%')
    S('%  Initial orbit:')
    S(f'%    SMA  = {p["SMA1"]:.3f} km    INC = {p["INC1"]:.4f} deg    RAAN = {p["RAAN1"]:.4f} deg')
    S('%  Final orbit:')
    S(f'%    SMA  = {p["SMA2"]:.3f} km    INC = {p["INC2"]:.4f} deg    RAAN = {p["RAAN2"]:.4f} deg')
    S('%')
    if p['two_burn']:
        S(f'%  Burn 1 at u = {p["u1"]:.4f} deg  (coast {p["prop1"]:.1f} s)')
        S(f'%    dVV = {p["dVV1"]:+.6f} km/s   dVN = {p["dVN1"]:+.6f} km/s')
        S(f'%    |dV| = {math.sqrt(p["dVV1"]**2+p["dVN1"]**2):.6f} km/s')
        if finite:
            S(f'%    Thrust direction: V={p["td1_1"]:.6f}  N={p["td2_1"]:.6f}')
            S(f'%    Burn duration est: {p["dur1"]:.1f} s  ({p["dur1"]/60:.1f} min)')
        S(f'%  Transfer half-period: {p["T_trans"]:.1f} s')
        S(f'%  Burn 2 at u = {(p["u1"]+180)%360:.4f} deg')
        S(f'%    dVV = {p["dVV2"]:+.6f} km/s   dVN = {p["dVN2"]:+.6f} km/s')
        S(f'%    |dV| = {math.sqrt(p["dVV2"]**2+p["dVN2"]**2):.6f} km/s')
        if finite:
            S(f'%    Thrust direction: V={p["td1_2"]:.6f}  N={p["td2_2"]:.6f}')
            S(f'%    Burn duration est: {p["dur2"]:.1f} s  ({p["dur2"]/60:.1f} min)')
    else:
        S(f'%  Single burn at u = {p["u1"]:.4f} deg  (coast {p["prop1"]:.1f} s)')
        S(f'%    dVN = {p["dVN1"]:+.6f} km/s  (pure out-of-plane)')
        if finite:
            S(f'%    Burn duration est: {p["dur1"]:.1f} s  ({p["dur1"]/60:.1f} min)')
    S(f'%  Total DeltaV = {p["dv_total"]:.4f} km/s')
    S(f'%  Fuel estimate = {p["fuel"]:.1f} kg  (Isp={isp:.0f}s)')
    if PROP_TYPE == 'chemical':
        S(f'%  Thrust = {CHEM_THRUST:.1f} N')
    elif PROP_TYPE == 'electric':
        S(f'%  Thrust = {ELEC_THRUST:.3f} N    Power = {ELEC_POWER:.1f} kW')
        S(f'%  NOTE: Electric burn durations are very long (days to weeks).')
        S(f'%  GMAT integrates the trajectory during the burn, so gravity')
        S(f'%  losses and orbit evolution are fully captured.')
    S('%----------------------------------------------------------')
    S()

    # ---- spacecraft ----
    S('Create Spacecraft SC;')
    S('SC.DateFormat        = UTCGregorian;')
    S(f'SC.Epoch             = \'{EPOCH}\';')
    S('SC.CoordinateSystem  = EarthMJ2000Eq;')
    S('SC.DisplayStateType  = Keplerian;')
    S(f'SC.SMA  = {p["SMA1"]:.6f};')
    S( 'SC.ECC  = 0.0001;')
    S(f'SC.INC  = {p["INC1"]:.6f};')
    S(f'SC.RAAN = {p["RAAN1"]:.6f};')
    S( 'SC.AOP  = 0.0;')
    S( 'SC.TA   = 0.0;')
    S(f'SC.DryMass  = {DRY_MASS:.1f};')
    S( 'SC.Cd       = 2.2;')
    S( 'SC.Cr       = 1.8;')
    S( 'SC.DragArea = 15.0;')
    S( 'SC.SRPArea  = 15.0;')
    S()

    # ---- tank ----
    if PROP_TYPE == 'chemical':
        S('Create ChemicalTank MainTank;')
        S('MainTank.AllowNegativeFuelMass = true;')
        S(f'MainTank.FuelMass              = {FUEL_MASS:.1f};')
        S('MainTank.Pressure              = 1500.0;')
        S('MainTank.Temperature           = 20.0;')
        S('MainTank.RefTemperature        = 20.0;')
        S('MainTank.Volume                = 0.5;')
        S('MainTank.FuelDensity           = 1000.0;')
        S('MainTank.PressureModel         = PressureRegulated;')
        S('SC.Tanks = {MainTank};')
    elif PROP_TYPE == 'electric':
        S('Create ElectricTank MainTank;')
        S('MainTank.AllowNegativeFuelMass = true;')
        S(f'MainTank.FuelMass              = {FUEL_MASS:.1f};')
        S('SC.Tanks = {MainTank};')
    else:
        # impulsive -- still needs a tank for mass decrement
        S('Create ChemicalTank MainTank;')
        S('MainTank.AllowNegativeFuelMass = true;')
        S(f'MainTank.FuelMass              = {FUEL_MASS:.1f};')
        S('MainTank.Pressure              = 1500.0;')
        S('MainTank.Temperature           = 20.0;')
        S('MainTank.RefTemperature        = 20.0;')
        S('MainTank.Volume                = 0.5;')
        S('MainTank.FuelDensity           = 1000.0;')
        S('MainTank.PressureModel         = PressureRegulated;')
        S('SC.Tanks = {MainTank};')
    S()

    # ---- power system (electric only) ----
    if PROP_TYPE == 'electric':
        S('% Electric power system required for ElectricThruster')
        S('Create SolarPowerSystem SPS;')
        S(f'SPS.InitialMaxPower          = {ELEC_POWER:.2f};   % kW')
        S( 'SPS.AnnualDecayRate          = 0.0;')
        S( 'SPS.Margin                   = 0.05;')
        S( 'SPS.BusCoeff1                = 0.3;')
        S( 'SPS.ShadowModel              = None;')
        S( 'SPS.ShadowBodies             = {};')
        S( 'SC.PowerSystem               = SPS;')
        S()

    # ---- thrusters (finite burn types only) ----
    n_burns = 2 if p['two_burn'] else 1
    if PROP_TYPE == 'chemical':
        for k in range(1, n_burns+1):
            td1 = p[f'td1_{k}']; td2 = p[f'td2_{k}']
            S(f'Create ChemicalThruster Thruster{k};')
            S(f'Thruster{k}.CoordinateSystem  = Local;')
            S(f'Thruster{k}.Origin            = Earth;')
            S(f'Thruster{k}.Axes              = VNB;')
            S(f'Thruster{k}.ThrustDirection1  = {td1:.6f};   % V')
            S(f'Thruster{k}.ThrustDirection2  = {td2:.6f};   % N')
            S(f'Thruster{k}.ThrustDirection3  = 0.0;          % B')
            S(f'Thruster{k}.DutyCycle         = 1.0;')
            S(f'Thruster{k}.ThrustScaleFactor = 1.0;')
            S(f'Thruster{k}.GravitationalAccel = {G0:.5f};')
            S(f'Thruster{k}.C1                = {CHEM_THRUST:.2f};  % N  constant thrust')
            S(f'Thruster{k}.C2                = 0.0;')
            S(f'Thruster{k}.C3                = 0.0;')
            S(f'Thruster{k}.K1                = {CHEM_ISP:.2f};  % s  Isp')
            S(f'Thruster{k}.K2                = 0.0;')
            S(f'Thruster{k}.K3                = 0.0;')
            S(f'Thruster{k}.MixRatio          = [1.0];')
            S(f'Thruster{k}.Tank              = {{MainTank}};')
            S()
        thrusters = '{' + ', '.join(f'Thruster{k}' for k in range(1,n_burns+1)) + '}'
        S(f'SC.Thrusters = {thrusters};')
        S()
    elif PROP_TYPE == 'electric':
        # Two thrusters with different directions for the two phases:
        # ThrusterV: in-track (V=1,N=0) for altitude raise
        # ThrusterN: out-of-plane (V=0,N=+/-1) for plane change
        ep = p['elec']
        S('% Electric Phase 1 thruster: pure in-track for altitude change')
        S('Create ElectricThruster ThrusterV;')
        S('ThrusterV.CoordinateSystem        = Local;')
        S('ThrusterV.Origin                  = Earth;')
        S('ThrusterV.Axes                    = VNB;')
        S(f'ThrusterV.ThrustDirection1        = {ep["td1_v"]:.1f};   % V (in-track)')
        S( 'ThrusterV.ThrustDirection2        = 0.0;')
        S( 'ThrusterV.ThrustDirection3        = 0.0;')
        S( 'ThrusterV.DutyCycle               = 1.0;')
        S( 'ThrusterV.ThrustScaleFactor       = 1.0;')
        S( 'ThrusterV.ThrustModel             = ConstantThrustAndIsp;')
        S(f'ThrusterV.MaximumUsablePower      = {ELEC_POWER:.3f};')
        S( 'ThrusterV.MinimumUsablePower      = 0.001;')
        S(f'ThrusterV.ConstantThrust          = {ELEC_THRUST:.4f};')
        S(f'ThrusterV.Isp                     = {ELEC_ISP:.2f};')
        S( 'ThrusterV.Tank                    = {MainTank};')
        S()
        S('% Electric Phase 2 thruster: pure out-of-plane for plane change')
        S('Create ElectricThruster ThrusterN;')
        S('ThrusterN.CoordinateSystem        = Local;')
        S('ThrusterN.Origin                  = Earth;')
        S('ThrusterN.Axes                    = VNB;')
        S( 'ThrusterN.ThrustDirection1        = 0.0;')
        S(f'ThrusterN.ThrustDirection2        = {ep["td2_n"]:.1f};   % N (out-of-plane)')
        S( 'ThrusterN.ThrustDirection3        = 0.0;')
        S( 'ThrusterN.DutyCycle               = 1.0;')
        S( 'ThrusterN.ThrustScaleFactor       = 1.0;')
        S( 'ThrusterN.ThrustModel             = ConstantThrustAndIsp;')
        S(f'ThrusterN.MaximumUsablePower      = {ELEC_POWER:.3f};')
        S( 'ThrusterN.MinimumUsablePower      = 0.001;')
        S(f'ThrusterN.ConstantThrust          = {ELEC_THRUST:.4f};')
        S(f'ThrusterN.Isp                     = {ELEC_ISP:.2f};')
        S( 'ThrusterN.Tank                    = {MainTank};')

        S('Create ElectricThruster ThrusterAntiN;')
        S('ThrusterAntiN.CoordinateSystem        = Local;')
        S('ThrusterAntiN.Origin                  = Earth;')
        S('ThrusterAntiN.Axes                    = VNB;')
        S( 'ThrusterAntiN.ThrustDirection1        = 0.0;')
        S(f'ThrusterAntiN.ThrustDirection2        = {-ep["td2_n"]:.1f};   % N (out-of-plane)')
        S( 'ThrusterAntiN.ThrustDirection3        = 0.0;')
        S( 'ThrusterAntiN.DutyCycle               = 1.0;')
        S( 'ThrusterAntiN.ThrustScaleFactor       = 1.0;')
        S( 'ThrusterAntiN.ThrustModel             = ConstantThrustAndIsp;')
        S(f'ThrusterAntiN.MaximumUsablePower      = {ELEC_POWER:.3f};')
        S( 'ThrusterAntiN.MinimumUsablePower      = 0.001;')
        S(f'ThrusterAntiN.ConstantThrust          = {ELEC_THRUST:.4f};')
        S(f'ThrusterAntiN.Isp                     = {ELEC_ISP:.2f};')
        S( 'ThrusterAntiN.Tank                    = {MainTank};')
        S()
        S('SC.Thrusters = {ThrusterV, ThrusterN, ThrusterAntiN};')
        S()

    # ---- burns ----
    if not finite:
        # Impulsive burns
        for k, (dVV, dVN, label) in enumerate(
                [(p['dVV1'],p['dVN1'],'Burn1'),(p['dVV2'],p['dVN2'],'Burn2')]
                if p['two_burn'] else
                [(0.0,p['dVN1'],'Burn1')], start=1):
            S(f'Create ImpulsiveBurn {label};')
            S(f'{label}.CoordinateSystem = Local;')
            S(f'{label}.Origin           = Earth;')
            S(f'{label}.Axes             = VNB;')
            S(f'{label}.Element1         = {dVV:.6f};  % V km/s')
            S(f'{label}.Element2         = {dVN:.6f};  % N km/s')
            S(f'{label}.Element3         = 0.0;')
            S(f'{label}.DecrementMass    = true;')
            S(f'{label}.Tank             = {{MainTank}};')
            S(f'{label}.Isp              = {CHEM_ISP:.1f};')
            S(f'{label}.GravitationalAccel = {G0:.5f};')
            S()
    else:
        # Finite burns
        if PROP_TYPE == 'chemical':
            n_finite = n_burns
            for k in range(1, n_finite+1):
                S(f'Create FiniteBurn Burn{k};')
                S(f'Burn{k}.Thrusters = {{Thruster{k}}};')
                S()
        else:  # electric: two named burns, one per thruster
            S('Create FiniteBurn BurnV;')
            S('BurnV.Thrusters = {ThrusterV};')
            S()
            S('Create FiniteBurn BurnN;')
            S('BurnN.Thrusters = {ThrusterN};')
            S()
            S('Create FiniteBurn BurnAntiN;')
            S('BurnAntiN.Thrusters = {ThrusterAntiN};')
            S()

    # ---- force model ----
    S('Create ForceModel FM;')
    S('FM.CentralBody               = Earth;')
    S('FM.PrimaryBodies             = {Earth};')
    S('FM.GravityField.Earth.Degree = 2;')
    S('FM.GravityField.Earth.Order  = 0;   % J2 disabled')
    S('FM.Drag.AtmosphereModel      = None;')
    S('FM.SRP                       = Off;')
    S('FM.RelativisticCorrection    = Off;')
    S('FM.ErrorControl              = RSSStep;')
    S()
    S('Create Propagator Prop;')
    S('Prop.FM              = FM;')
    S('Prop.Type            = RungeKutta89;')
    S('Prop.InitialStepSize = 60;')
    S('Prop.Accuracy        = 1e-11;')
    S('Prop.MinStep         = 0.001;')
    S('Prop.MaxStep         = 900;')
    S('Prop.MaxStepAttempts = 50;')
    S()

    # ---- DC ----
    S('Create DifferentialCorrector DC;')
    S('DC.ShowProgress      = true;')
    S('DC.ReportStyle       = Verbose;')
    S('DC.ReportFile        = \'DCReport.txt\';')
    # Electric phase 2 solves 2 coupled goals (INC + RAAN) over a months-long
    # low-thrust arc, so allow more iterations than the impulsive/chemical case.
    S(f'DC.MaximumIterations = {100 if PROP_TYPE == "electric" else 50};')
    S('DC.DerivativeMethod  = ForwardDifference;')
    S('DC.Algorithm         = NewtonRaphson;')
    S()

    # ---- output ----
    S('Create ReportFile RF;')
    S('RF.SolverIterations = None;')
    S('RF.Filename         = \'PlaneChange_Report.txt\';')
    S('RF.Add              = {SC.UTCGregorian, SC.Earth.SMA, SC.ECC, SC.INC, SC.RAAN, SC.AOP, SC.TA, SC.TotalMass};')
    S('RF.WriteHeaders     = true;')
    S('RF.LeftJustify      = On;')
    S('RF.ZeroFill         = Off;')
    S('RF.ColumnWidth      = 23;')
    S('RF.Precision        = 10;')
    S()
    S('Create OrbitView OrbitView1;')
    S('OrbitView1.SolverIterations       = None;')
    S('OrbitView1.UpperLeft              = [ 0.0 0.0 ];')
    S('OrbitView1.Size                   = [ 0.8 0.8 ];')
    S('OrbitView1.RelativeZOrder         = 50;')
    S('OrbitView1.Maximized              = false;')
    S('OrbitView1.Add                    = {SC, Earth};')
    S('OrbitView1.CoordinateSystem       = EarthMJ2000Eq;')
    S('OrbitView1.DrawObject             = [ true true ];')
    S('OrbitView1.DataCollectFrequency   = 1;')
    S('OrbitView1.UpdatePlotFrequency    = 50;')
    S('OrbitView1.NumPointsToRedraw      = 0;')
    S('OrbitView1.ShowPlot               = true;')
    S('OrbitView1.MaxPlotPoints          = 20000;')
    S('OrbitView1.ShowLabels             = true;')
    S('OrbitView1.ViewPointReference     = Earth;')
    S('OrbitView1.ViewPointVector        = [ 0 0 90000 ];')
    S('OrbitView1.ViewDirection          = Earth;')
    S('OrbitView1.ViewScaleFactor        = 1.0;')
    S('OrbitView1.ViewUpCoordinateSystem = EarthMJ2000Eq;')
    S('OrbitView1.ViewUpAxis             = Z;')
    S('OrbitView1.EclipticPlane          = Off;')
    S('OrbitView1.XYPlane                = On;')
    S('OrbitView1.WireFrame              = Off;')
    S('OrbitView1.Axes                   = On;')
    S('OrbitView1.Grid                   = Off;')
    S('OrbitView1.SunLine                = Off;')
    S('OrbitView1.UseInitialView         = On;')
    S('OrbitView1.StarCount              = 7000;')
    S('OrbitView1.EnableStars            = On;')
    S('OrbitView1.EnableConstellations   = On;')
    S()

    # ---- variables ----
    if PROP_TYPE == 'chemical' and p['two_burn']:
        S(f'Create Variable dur1 dur2 t_trans;')
        S(f'dur1 = {p["dur1"]:.1f};   % s  Burn 1 duration initial guess')
        S(f'dur2 = {p["dur2"]:.1f};   % s  Burn 2 duration initial guess')
        S(f't_trans = {p["T_trans"]:.1f};   % s  Burn 2 duration initial guess')
        S()
    elif PROP_TYPE == 'electric':
        ep = p['elec']
        S('Create Variable durV;')
        S(f'durV = {ep["t1_dur"]:.1f};  % s  Phase 1 altitude burn duration')
        S('Create Variable N_ARCS;')
        S(f'N_ARCS = {ep["n_arcs"]};    % Phase 2 arc count   (DC-varied -> INC)')
        S('Create Variable ta_start;')
        S(f'ta_start = {ep["ta_start"]:.4f};  % deg +N arc start TA (DC-varied -> RAAN)')
        S('Create Variable I;')
        S()

    # ====================================================
    # Mission sequence
    # ====================================================
    S('BeginMissionSequence;')
    S()
    S('Report RF SC.UTCGregorian SC.Earth.SMA SC.ECC SC.INC SC.RAAN SC.AOP SC.TA SC.TotalMass;')
    S()

    if PROP_TYPE == 'electric':
        ep = p['elec']
        S('%--------------------------------------------------')
        S('% ELECTRIC PHASE 1: altitude raise (continuous V burn)')
        S('% DC varies duration to achieve target SMA.')
        S('%--------------------------------------------------')
        if not ep['same_alt']:
            S('Target DC {SolveMode = Solve, ExitMode = DiscardAndContinue};')
            S()
            S(f'   Vary DC(durV = {ep["t1_dur"]:.1f}, {{Perturbation = {max(3600.0,ep["t1_dur"]*0.01):.1f}, Lower = 1.0, Upper = 1e9, MaxStep = {max(86400.0,ep["t1_dur"]*0.3):.1f}}});')
            S()
            S('   BeginFiniteBurn BurnV(SC);')
            S('   Propagate Prop(SC) {SC.ElapsedSecs = durV};')
            S('   EndFiniteBurn BurnV(SC);')
            S()
            S('   Propagate Prop(SC) {SC.ElapsedSecs = 600.0};')
            S()
            S(f'   Achieve DC(SC.Earth.SMA = {p["SMA2"]:.3f}, {{Tolerance = 1.0}});')
            S()
            S('EndTarget;')
            S()
            S('Report RF SC.UTCGregorian SC.Earth.SMA SC.ECC SC.INC SC.RAAN SC.AOP SC.TA SC.TotalMass;')
            S()
        S('%--------------------------------------------------')
        S('% ELECTRIC PHASE 2: plane change via dual half-orbit burns')
        S(f'% +N over {ep["ARC_DEG"]:.0f} deg centred on the node line, then -N over the')
        S('% antipodal half. Both halves add inclination (the burn sign flips')
        S('% with cos(u)) AND steer the node, so no time is spent coasting.')
        S(f'% Seed: N_ARCS={ep["n_arcs"]} (~{ep["t2_elapsed"]/86400:.1f} days), arc start TA={ep["ta_start"]:.1f} deg.')
        S('%')
        S('% CLOSED-LOOP targeting (2 vars, 2 goals, nearly orthogonal):')
        S('%   Vary N_ARCS    (orbit count)    -> Achieve SC.INC   (magnitude)')
        S('%   Vary ta_start  (arc-centre TA)  -> Achieve SC.RAAN  (direction)')
        S('%--------------------------------------------------')
        S('Target DC {SolveMode = Solve, ExitMode = DiscardAndContinue};')
        S()
        S(f'   Vary DC(N_ARCS   = {ep["n_arcs"]}, {{Perturbation = 1, Lower = 1, Upper = 10000, MaxStep = 10}});')
        S(f'   Vary DC(ta_start = {ep["ta_start"]:.4f}, {{Perturbation = 0.5, Lower = 0.0, Upper = 360.0, MaxStep = 30.0}});')
        S()
        S('   % Phase the orbit so the +N arc is centred on the node line.')
        S('   Propagate Prop(SC) {SC.TA = ta_start};')
        S('   For I = 1 : 1 : N_ARCS;')
        S('      % +N half-orbit: build plane change across the node line')
        S('      BeginFiniteBurn BurnN(SC);')
        S(f'      Propagate Prop(SC) {{SC.ElapsedSecs = {ep["t_pos_arc"]:.1f}}};')
        S('      EndFiniteBurn BurnN(SC);')
        S('      % -N half-orbit: opposite sign across the antinode adds the')
        S('      % SAME inclination instead of coasting, ~halving total time.')
        S('      BeginFiniteBurn BurnAntiN(SC);')
        S(f'      Propagate Prop(SC) {{SC.ElapsedSecs = {ep["t_neg_arc"]:.1f}}};')
        S('      EndFiniteBurn BurnAntiN(SC);')
        S('   EndFor;')
        S()
        S('   Propagate Prop(SC) {SC.ElapsedSecs = 600.0};')
        S()
        S(f'   Achieve DC(SC.INC  = {p["INC2"]:.6f}, {{Tolerance = 0.05}});')
        S(f'   Achieve DC(SC.RAAN = {p["RAAN2"]:.6f}, {{Tolerance = 0.1}});')
        S()
        S('EndTarget;')
    else:
        S(f'% Coast to burn 1: u = {p["u1"]:.4f} deg')
        S(f'Propagate Prop(SC) {{SC.ElapsedSecs = {p["prop1"]:.1f}}};')
        S()
        S('Report RF SC.UTCGregorian SC.Earth.SMA SC.ECC SC.INC SC.RAAN SC.AOP SC.TA SC.TotalMass;')
        S()
        # ---- impulsive / chemical targeting ----
        if p['two_burn']:
            if not finite:
                S('Target DC {SolveMode = Solve, ExitMode = DiscardAndContinue};')
                S()
                S(f'   Vary DC(Burn1.Element1 = {p["dVV1"]:.6f}, {{Perturbation = 1e-5, Lower = -5.0, Upper = 5.0, MaxStep = 0.02}});')
                S(f'   Vary DC(Burn1.Element2 = {p["dVN1"]:.6f}, {{Perturbation = 1e-5, Lower = -5.0, Upper = 5.0, MaxStep = 0.02}});')
                S(f'   Vary DC(Burn2.Element1 = {p["dVV2"]:.6f}, {{Perturbation = 1e-5, Lower = -5.0, Upper = 5.0, MaxStep = 0.02}});')
                S(f'   Vary DC(Burn2.Element2 = {p["dVN2"]:.6f}, {{Perturbation = 1e-5, Lower = -5.0, Upper = 5.0, MaxStep = 0.02}});')
                S()
                S('   Maneuver Burn1(SC);')
                S()
                S(f'   % Coast transfer orbit to apogee: {p["T_trans"]:.1f} s (half transfer period)')
                S(f'   Propagate Prop(SC) {{SC.ElapsedSecs = {p["T_trans"]:.1f}}};')
                S()
                S('   Maneuver Burn2(SC);')
                S()
                S('   Propagate Prop(SC) {SC.ElapsedSecs = 10.0};')
                S()
                S(f'   Achieve DC(SC.Earth.SMA = {p["SMA2"]:.3f}, {{Tolerance = 0.1}});')
                S(f'   Achieve DC(SC.INC       = {p["INC2"]:.6f},  {{Tolerance = 0.01}});')
                S(f'   Achieve DC(SC.RAAN      = {p["RAAN2"]:.6f}, {{Tolerance = 0.01}});')
                S()
                S('EndTarget;')
            else:
                # chemical two-burn finite
                S('% Finite burn targeting: DC varies burn durations.')
                S('% Thrust directions are fixed analytically on Thruster1/2.')
                S('% Vary: dur1, dur2  --  Achieve: SMA, INC')
                S('Target DC {SolveMode = Solve, ExitMode = DiscardAndContinue};')
                S()
                S(f'   Vary DC(dur1 = {p["dur1"]:.1f}, {{Perturbation = 1.0, Lower = 1.0, Upper = 1e8, MaxStep = {max(60.0, p["dur1"]*0.1):.1f}}});')
                S(f'   Vary DC(dur2 = {p["dur2"]:.1f}, {{Perturbation = 1.0, Lower = 1.0, Upper = 1e8, MaxStep = {max(60.0, p["dur2"]*0.1):.1f}}});')
                S(f'   Vary DC(t_trans = {p["T_trans"]:.1f}, {{Perturbation = 1.0, Lower = 1.0, Upper = 1e8, MaxStep = {max(60.0, p["dur2"]*0.1):.1f}}});')
                S()
                S('   BeginFiniteBurn Burn1(SC);')
                S('   Propagate Prop(SC) {SC.ElapsedSecs = dur1};')
                S('   EndFiniteBurn Burn1(SC);')
                S()
                S(f'   % Coast transfer orbit to apogee: {p["T_trans"]:.1f} s (half transfer period)')
                S(f'   Propagate Prop(SC) {{SC.ElapsedSecs = t_trans}};')
                S()
                S('   BeginFiniteBurn Burn2(SC);')
                S('   Propagate Prop(SC) {SC.ElapsedSecs = dur2};')
                S('   EndFiniteBurn Burn2(SC);')
                S()
                S('   Propagate Prop(SC) {SC.ElapsedSecs = 60.0};')
                S()
                S(f'   Achieve DC(SC.Earth.SMA = {p["SMA2"]:.3f}, {{Tolerance = 0.5}});')
                S(f'   Achieve DC(SC.INC       = {p["INC2"]:.6f},  {{Tolerance = 0.01}});')
                S(f'   Achieve DC(SC.RAAN      = {p["RAAN2"]:.6f}, {{Tolerance = 0.01}});')
                S()
                S('EndTarget;')
        else:
            # single burn (impulsive or chemical finite, same altitude)
            if not finite:
                S('Target DC {SolveMode = Solve, ExitMode = DiscardAndContinue};')
                S()
                S(f'   Vary DC(Burn1.Element2 = {p["dVN1"]:.6f}, {{Perturbation = 1e-5, Lower = -2.0, Upper = 2.0, MaxStep = 0.05}});')
                S()
                S('   Maneuver Burn1(SC);')
                S()
                S('   Propagate Prop(SC) {SC.ElapsedSecs = 10.0};')
                S()
                S(f'   Achieve DC(SC.RAAN = {p["RAAN2"]:.6f}, {{Tolerance = 0.01}});')
                S(f'   Achieve DC(SC.INC  = {p["INC2"]:.6f},  {{Tolerance = 0.01}});')
                S()
                S('EndTarget;')
            else:
                # chemical single finite burn: open-loop
                S('% Single chemical finite burn: open-loop (no targeter).')
                S('% Direction and duration are set analytically.')
                S('BeginFiniteBurn Burn1(SC);')
                S(f'Propagate Prop(SC) {{SC.ElapsedSecs = {p["dur1"]:.1f}}};')
                S('EndFiniteBurn Burn1(SC);')
                S()
                S('Propagate Prop(SC) {SC.ElapsedSecs = 60.0};')

    S()
    S('Report RF SC.UTCGregorian SC.Earth.SMA SC.ECC SC.INC SC.RAAN SC.AOP SC.TA SC.TotalMass;')
    S()
    S('Propagate Prop(SC) {SC.ElapsedDays = 1.0};')
    S()
    S('Report RF SC.UTCGregorian SC.Earth.SMA SC.ECC SC.INC SC.RAAN SC.AOP SC.TA SC.TotalMass;')

    with open(path,'w',encoding='ascii') as f:
        f.write('\n'.join(L)+'\n')

# -----------------------------------------------------------
# Main
# -----------------------------------------------------------
def main():
    if PROP_TYPE not in ('impulsive','chemical','electric'):
        print(f'ERROR: PROP_TYPE must be impulsive, chemical, or electric. Got: {PROP_TYPE}')
        sys.exit(1)

    print('='*62)
    print('  GMAT Orbital Maneuver Script Generator')
    print('='*62)
    print(f'  Propulsion : {PROP_TYPE.upper()}')
    if PROP_TYPE == 'chemical':
        print(f'  Thrust     : {CHEM_THRUST:.1f} N     Isp : {CHEM_ISP:.1f} s')
    elif PROP_TYPE == 'electric':
        print(f'  Thrust     : {ELEC_THRUST:.3f} N    Isp : {ELEC_ISP:.1f} s    Power : {ELEC_POWER:.1f} kW')
    print(f'  Initial    : SMA={SMA_init} km  INC={INC_init} deg  RAAN={RAAN_init} deg')
    print(f'  Final      : SMA={SMA_final} km  INC={INC_final} deg  RAAN={RAAN_final} deg')
    print()

    if (abs(INC_final-INC_init)<1e-6 and abs(RAAN_final-RAAN_init)<1e-6
            and abs(SMA_final-SMA_init)<0.1):
        print('ERROR: Orbits are identical.'); sys.exit(1)

    isp    = active_isp()
    thrust = active_thrust_N()
    finite = PROP_TYPE in ('chemical','electric')
    i1=math.radians(INC_init); O1=math.radians(RAAN_init)
    T1=2*math.pi*math.sqrt(SMA_init**3/MU)
    same_alt=abs(SMA_final-SMA_init)<0.1
    m0=DRY_MASS+FUEL_MASS

    if same_alt:
        print('  Same-altitude: single-burn pure plane change.')
        u1,dVN1,conv=solve_single_burn(SMA_init,INC_init,RAAN_init,INC_final,RAAN_final)
        dVV1=0.0; dVV2=0.0; dVN2=0.0
        prop1=(u1%360)/360*T1; T_trans=0.0
        dv_total=abs(dVN1)
        fuel=fuel_kg(dv_total,m0,isp)
        two_burn=False
        td1_1,td2_1=thrust_direction(dVV1,dVN1)
        td1_2,td2_2=0.0,1.0
        dur1=burn_duration_s(dv_total,m0,thrust,isp) if finite else 0.0
        dur2=0.0
        # electric: two-phase arc-burn parameters (same-alt case, no phase 1)
        elec = electric_params(SMA_init,INC_init,RAAN_init,SMA_final,INC_final,RAAN_final,
                               ELEC_THRUST,ELEC_ISP,m0) if PROP_TYPE=='electric' else {}
    else:
        print('  Altitude change: two-burn Hohmann + plane change.')
        i2=math.radians(INC_final); O2=math.radians(RAAN_final)
        alpha=plane_angle(i1,O1,i2,O2)
        a1_opt,a2_opt,vtp,vta=optimal_split(SMA_init,SMA_final,alpha)
        a_t=(SMA_init+SMA_final)/2; v1c=math.sqrt(MU/SMA_init); v2c=math.sqrt(MU/SMA_final)
        T_trans=math.pi*math.sqrt(a_t**3/MU)
        dVV1=vtp*math.cos(a1_opt)-v1c; dVN1=vtp*math.sin(a1_opt)
        dVV2=v2c*math.cos(a2_opt)-vta; dVN2=v2c*math.sin(a2_opt)
        print('  Solving for burn location...')
        u1,dVV1,dVN1,dVV2,dVN2,conv=solve_two_burn(
            SMA_init,INC_init,RAAN_init,SMA_final,INC_final,RAAN_final,
            dVV1,dVN1,dVV2,dVN2)
        prop1=(u1%360)/360*T1
        dv1=math.sqrt(dVV1**2+dVN1**2); dv2=math.sqrt(dVV2**2+dVN2**2)
        dv_total=dv1+dv2
        fuel=fuel_kg(dv_total,m0,isp)
        two_burn=True
        td1_1,td2_1=thrust_direction(dVV1,dVN1)
        td1_2,td2_2=thrust_direction(dVV2,dVN2)
        if finite:
            dur1=burn_duration_s(dv1,m0,thrust,isp)
            m_after1=m0-fuel_kg(dv1,m0,isp)
            dur2=burn_duration_s(dv2,m_after1,thrust,isp)
        else:
            dur1=dur2=0.0
        # electric: two-phase arc-burn parameters
        elec = electric_params(SMA_init,INC_init,RAAN_init,SMA_final,INC_final,RAAN_final,
                               ELEC_THRUST,ELEC_ISP,m0) if PROP_TYPE=='electric' else {}

    if not conv:
        print('  WARNING: solver did not fully converge. Check GMAT output.')
    else:
        print('  Converged.')

    print()
    print(f'  {"--- Burn 1 ---"}')
    print(f'  u1         = {u1:.4f} deg   coast = {prop1:.1f} s ({prop1/3600:.3f} h)')
    print(f'  dVV (V)    = {dVV1:+.6f} km/s')
    print(f'  dVN (N)    = {dVN1:+.6f} km/s')
    print(f'  |dV1|      = {math.sqrt(dVV1**2+dVN1**2):.6f} km/s')
    if finite:
        print(f'  Direction  : V={td1_1:.4f}  N={td2_1:.4f}')
        print(f'  Duration   = {dur1:.1f} s ({dur1/60:.1f} min)')
    if two_burn:
        print()
        print(f'  {"--- Transfer ---"}')
        print(f'  Half-period = {T_trans:.1f} s ({T_trans/3600:.3f} h)')
        print()
        print(f'  {"--- Burn 2 ---"}')
        print(f'  u2         = {(u1+180)%360:.4f} deg')
        print(f'  dVV (V)    = {dVV2:+.6f} km/s')
        print(f'  dVN (N)    = {dVN2:+.6f} km/s')
        print(f'  |dV2|      = {math.sqrt(dVV2**2+dVN2**2):.6f} km/s')
        if finite:
            print(f'  Direction  : V={td1_2:.4f}  N={td2_2:.4f}')
            print(f'  Duration   = {dur2:.1f} s ({dur2/60:.1f} min)')
    print()
    print(f'  Total |DeltaV| = {dv_total:.4f} km/s')
    print(f'  Fuel estimate  = {fuel:.1f} kg  (Isp={isp:.0f}s)')
    if finite and PROP_TYPE=='electric':
        print(f'  NOTE: Electric burn durations are very long -- see script comments.')

    if PROP_TYPE == 'electric':
        print()
        print('  --- Electric two-phase plan ---')
        if not elec.get('same_alt',True):
            print(f'  Phase 1 (altitude): continuous V burn')
            print(f'    Duration   = {elec["t1_dur"]:.1f} s ({elec["t1_dur"]/86400:.2f} days)')
            print(f'    Fuel       = {elec["dm1"]:.1f} kg')
        print(f'  Phase 2 (plane chg): {elec["n_arcs"]} arcs of {elec["ARC_DEG"]:.0f} deg (seed)')
        print(f'    Arc centre : TA = {elec["u_opt"]:.1f} deg  (start TA = {elec["ta_start"]:.1f} deg)')
        print(f'    N sign     : {elec["td2_n"]:+.0f}  ({"raise" if elec["td2_n"]>0 else "lower"} inclination)')
        print(f'    Elapsed    ~ {elec["t2_elapsed"]/86400:.1f} days')
        print(f'    Fuel       = {elec["dm2"]:.1f} kg')
        print(f'  Total fuel   = {elec["fuel_total"]:.1f} kg  (Isp={isp:.0f}s)')
        print(f'  Closed-loop  : Vary N_ARCS->INC, Vary ta_start->RAAN')

    params=dict(
        SMA1=SMA_init,INC1=INC_init,RAAN1=RAAN_init,
        SMA2=SMA_final,INC2=INC_final,RAAN2=RAAN_final,
        u1=u1,prop1=prop1,T_trans=T_trans,
        dVV1=dVV1,dVN1=dVN1,dVV2=dVV2,dVN2=dVN2,
        td1_1=td1_1,td2_1=td2_1,td1_2=td1_2,td2_2=td2_2,
        dur1=dur1,dur2=dur2,
        elec=elec,
        dv_total=dv_total,fuel=fuel,two_burn=two_burn)

    out=os.path.join(os.path.dirname(os.path.abspath(__file__)),OUTPUT)
    write_script(out,params)
    print()
    print(f'  Script written to: {out}')
    print('='*62)

if __name__=='__main__':
    main()