"""
REAVER RPO - Simulation Entry Point
===================================
Runs Part 1 (debris RPO) and Part 2 (RH docking cycles), prints a phase-by-
phase breakdown and the tumble-rate sweep, then emits the five RPO dV constants
ready to paste into config.py (Part 3 of the build checklist).

    cd c:\\Projects\\DSE\\REAVER
    python trajectory_optimizer/rpo_run.py
"""
import numpy as np

from config import MS_DRY, MS_ISP, MS_VEX
from rpo_config import OMEGA_REQ_DPS, TUMBLE_SWEEP_DPS
from rpo_debris_sim import simulate_debris_rpo, sweep_tumble_rates
from rpo_rh_sim import simulate_all_rh_cycles

# EUTE 12 WEST A (worst-case target) is catalogue index 6.  Placing it first in
# the visit order gives the heaviest chaser mass coincident with the worst-case
# target — the conservative posture for the debris-RPO propellant headline.
WORST_CASE_SEQ = [6, 0, 1, 2, 3]


def debris_chaser_masses(seq=WORST_CASE_SEQ):
    """
    Real MS mass at each debris arrival, taken from the optimizer's forward-pass
    mass timeline (mass_vec[i] = mass on arrival at debris i, before tug release).
    mass_vec[0] is 'wet mass minus first-transfer propellant'.
    """
    import nominal_mission as nm
    r = nm.evaluate_sequence(list(seq))
    return np.asarray(r['mass_vec'][:5], float)


def _rule(char="-", n=72):
    print(char * n)


def report_debris():
    print("\n" + "=" * 72)
    print("PART 1 - DEBRIS RPO  (EUTE 12 WEST A, 2700 kg, worst case)")
    print("=" * 72)

    masses = debris_chaser_masses()
    m_worst = float(masses[0])      # heaviest chaser (first debris) -> most propellant
    res = simulate_debris_rpo(omega_dps=OMEGA_REQ_DPS, chaser_mass=m_worst)

    print(f"\nTumble rate: {OMEGA_REQ_DPS:.1f} deg/s (1 rpm)   "
          f"chaser mass: {m_worst:.0f} kg (worst-case, first debris)")
    _rule()
    print(f"{'Phase':<34s} {'dV':>8s}        {'Prop':>7s}      {'Time':>7s}")
    _rule()
    for p in res['phases']:
        print(p)
    _rule()
    print(f"{'DV_RPO_DEBRIS (incl. +50% abort)':<34s} "
          f"{res['dv_debris']:8.3f} m/s   {res['prop_debris']:7.4f} kg")
    print(f"{'DV_RPO_DETUMBLE (momentum dump)':<34s} "
          f"{res['dv_detumble']:8.3f} m/s   {res['prop_detumble']:7.4f} kg")
    print(f"{'Total phase duration':<34s} {res['duration']/3600:8.2f} h")
    _rule()

    print("\nPer-leg propellant (real arrival masses from mass timeline):")
    tot_prop = 0.0
    for i, m in enumerate(masses, 1):
        leg = simulate_debris_rpo(omega_dps=OMEGA_REQ_DPS, chaser_mass=float(m))
        prop = leg['prop_debris'] + leg['prop_detumble']
        tot_prop += prop
        print(f"  Debris {i}: chaser {m:7.0f} kg  ->  "
              f"dV {leg['dv_debris'] + leg['dv_detumble']:5.2f} m/s   prop {prop:5.3f} kg")
    print(f"  Total debris-RPO propellant over 5 legs: {tot_prop:.2f} kg")
    _rule()

    print("\nTumble-rate sweep (Conings & Mooij Table 5 benchmark):")
    for r in sweep_tumble_rates(TUMBLE_SWEEP_DPS):
        print(f"  w = {r['omega_dps']:4.1f} deg/s  ->  "
              f"DV_DEBRIS = {r['dv_debris']:6.3f} m/s   "
              f"DV_DETUMBLE = {r['dv_detumble']:6.3f} m/s   "
              f"prop = {r['prop_debris'] + r['prop_detumble']:6.3f} kg")

    return res


def report_rh():
    print("\n" + "=" * 72)
    print("PART 2 - RECYCLING-HUB RPO  (MS x5 docking cycles)")
    print("=" * 72)

    cycles = simulate_all_rh_cycles(ms_mass_floor=MS_DRY, n_cycles=5)

    _rule()
    print(f"{'Cycle':>5s} {'MS mass':>9s} {'MEET dV':>9s} {'DOCK dV':>9s} "
          f"{'RH dV':>8s} {'Prop':>8s}")
    _rule()
    tot_meet = tot_dock = tot_prop = 0.0
    for i, c in enumerate(cycles, 1):
        prop = c['prop_meet'] + c['prop_dock']
        tot_meet += c['dv_meet']; tot_dock += c['dv_dock']; tot_prop += prop
        print(f"{i:>5d} {c['ms_mass']:8.0f}k {c['dv_meet']:8.3f} "
              f"{c['dv_dock']:8.3f} {c['dv_rh']:8.3f} {prop:7.4f}k")
    _rule()
    n = len(cycles)
    meet_avg = tot_meet / n
    dock_avg = tot_dock / n
    print(f"{'avg':>5s} {'':9s} {meet_avg:8.3f} {dock_avg:8.3f} "
          f"{meet_avg + dock_avg:8.3f} {tot_prop:7.4f}k")
    _rule()
    print(f"\nCycle 1 (heaviest) RH dV : {cycles[0]['dv_rh']:.3f} m/s")
    print(f"Cycle {n} (lightest) RH dV : {cycles[-1]['dv_rh']:.3f} m/s")

    return {'meet_avg': meet_avg, 'dock_avg': dock_avg,
            'rh_avg': meet_avg + dock_avg, 'total_prop': tot_prop}


def emit_constants(debris, rh):
    dv_debris = debris['dv_debris']
    dv_detumble = debris['dv_detumble']
    dv_meet = rh['meet_avg']
    dv_dock = rh['dock_avg']
    dv_rh = dv_meet + dv_dock

    print("\n" + "=" * 72)
    print("PART 3 - RPO dV CONSTANTS  (paste into config.py)")
    print("=" * 72)
    print(f"DV_RPO_DEBRIS   = {dv_debris:6.2f}   # m/s  MS alone, tumbling 2700 kg @ 1 rpm (+50% abort)")
    print(f"DV_RPO_DETUMBLE = {dv_detumble:6.2f}   # m/s  combined body momentum dump (REQ-ACS-M4)")
    print(f"DV_RPO_RH_MEET  = {dv_meet:6.2f}   # m/s  MS alone meets tug+debris within 500 m")
    print(f"DV_RPO_RH_DOCK  = {dv_dock:6.2f}   # m/s  MS+tug+debris to RH port (CoM offset dominant)")
    print(f"DV_RPO_RH       = {dv_rh:6.2f}   # m/s  total per RH cycle")

    print("\nConsistency checks:")
    debris_total = dv_debris + dv_detumble
    ok_order = dv_rh > debris_total
    print(f"  DV_RPO_RH ({dv_rh:.2f}) > DV_RPO_DEBRIS+DETUMBLE "
          f"({debris_total:.2f})  ->  {'PASS' if ok_order else 'FAIL'}")

    total_rpo_dv = 5 * debris_total + 5 * dv_rh
    m_prop = MS_DRY * (np.exp(total_rpo_dv / MS_VEX) - 1.0)
    print(f"  Total mission RPO dV : 5x({debris_total:.1f}) + 5x({dv_rh:.1f}) "
          f"= {total_rpo_dv:.1f} m/s")
    print(f"  Upper-bound propellant @ MS_DRY={MS_DRY:.0f} kg, "
          f"Isp={MS_ISP:.0f}s : {m_prop:.1f} kg")
    _rule("=", 72)


def main():
    debris = report_debris()
    rh = report_rh()
    emit_constants(debris, rh)


if __name__ == "__main__":
    main()
