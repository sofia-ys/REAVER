'''here is just all the copy/paste pieces to edit the other files to get what we want'''

'''VARYING DELTA V'''
# trajectory_optimizer\reaver_optimizer.py 
# if __name__ == '__main__':

#     t0=time.time()
#     res = evaluate_all_sequences()

#     if res:
#         nf  = res['n_feasible']
#         wp  = int(np.argmax(res['prop_used']))
#         wtp = int(np.argmax(res['total_prop']))

#         import pandas as pd
#         dv_table = pd.DataFrame({
#             "mission_days": res["mission_day"],
#             "total_delta_v_m_s": res["total_dv"]
#         })
#         dv_table.to_csv(
#             "mission_days_and_delta_v_table.csv",
#             sep=",",
#             index=False,
#             decimal="."
#         )
#         print("Saved mission days and delta-v table to mission_days_and_delta_v_table.csv")

#     print("\n  Done.")

'''VARYING MS DRY MASS'''
# trajectory_optimizer\nominal_mission.py
# LINE 154
# # t_legs  = np.array([T_LEG_FINITE[nf[i], nt[i]] for i in range(6)])
#     '''MASS VARYING VERSION'''
#     t_legs = np.zeros(6)

#     mass = ms_wet0
#     mass_vec = []
#     t = 0.0
#     tug_starts = []

#     for i in range(6):
#         fi = nf[i]
#         ti = nt[i]

#         d1 = DV1[fi, ti]
#         d2 = DV2[fi, ti]
#         dph = DV_PH[fi, ti]

#         leg_time = 0.0

#         # Burn 1
#         burn_time = finite_burn_time(d1, mass, SMA_ALL[fi])
#         t += burn_time
#         leg_time += burn_time
#         mass = mass * np.exp(-d1 / MS_VEX)

#         # Transfer coast
#         t += T_TR[fi, ti]
#         leg_time += T_TR[fi, ti]

#         # Burn 2
#         burn_time = finite_burn_time(d2, mass, SMA_ALL[ti])
#         t += burn_time
#         leg_time += burn_time
#         mass = mass * np.exp(-d2 / MS_VEX)

#         # Burn 3
#         burn_time = finite_burn_time(dph, mass, SMA_ALL[ti])
#         t += burn_time
#         leg_time += burn_time
#         mass = mass * np.exp(-dph / MS_VEX)

#         # Phasing coast
#         t += T_PH[fi, ti]
#         leg_time += T_PH[fi, ti]

#         t_legs[i] = leg_time
#         mass_vec.append(mass)

#         if i < 5:
#             mass = mass * np.exp(-(DV_RPO_DEBRIS + DV_RPO_DETUMBLE) / MS_VEX)
#             t += T_OPS
#             tug_starts.append(t)
#             mass -= tug_mwets[i]


# if __name__ == '__main__':

#     import pandas as pd
#     ms_dry_mass_range = range(698, 957, 10)
#     all_tables = []   # collect results from every MS_DRY value here
#     for MS_DRY in ms_dry_mass_range:
#         print(f"===========DONT GIVE UP: {MS_DRY}=============")

#         t0=time.time()
#         res = evaluate_all_sequences()

#         if res:
#             nf  = res['n_feasible']
#             wp  = int(np.argmax(res['prop_used']))
#             wtp = int(np.argmax(res['total_prop']))
            
#             ms_dry_mass_table = pd.DataFrame({
#                 "ms_dry_mass_kg": MS_DRY,
#                 "mission_days": res["mission_day"]
#             })
#             all_tables.append(ms_dry_mass_table)

#     final_table = pd.concat(all_tables, ignore_index=True)
#     final_table.to_csv(
#         "mission_days_for_ms_dry_mass_range.csv",
#         sep=",",
#         index=False,
#         decimal="."
#     )
#     print("Saved mission days table to mission_days_for_ms_dry_mass_range.csv")

#     print("\n  Done.")


'''VARYING TUG DRY MASS'''
# trajectory_optimizer\reaver_core.py
# def compute_tug_spirals(rh_raan_deg=RH_RAAN, tug_dry=TUG_DRY):
#     tug_dv    = np.zeros(N_DEB)
#     tug_time  = np.zeros(N_DEB)
#     tug_mprop = np.zeros(N_DEB)
#     tug_mwet  = np.zeros(N_DEB)

#     v2 = np.sqrt(MU / RH_SMA)
#     i2, o2 = RH_INC, rh_raan_deg

#     for k in range(N_DEB):
#         m_pl  = tug_dry + MASS[k]
#         m_wet = m_pl * 1.35

#         i1, o1 = INC[k], RAAN[k]
#         v1 = np.sqrt(MU / SMA[k])

#         cos_d = (
#             np.cos(i1*D2R)*np.cos(i2*D2R)
#             + np.sin(i1*D2R)*np.sin(i2*D2R)*np.cos((o2-o1)*D2R)
#         )

#         dth  = np.arccos(np.clip(cos_d, -1.0, 1.0))
#         dv_e = np.sqrt(v1**2 + v2**2 - 2*v1*v2*np.cos(np.pi/2*dth))

#         for _ in range(60):
#             mf   = m_wet * np.exp(-dv_e / TUG_VEX)
#             mnew = m_pl + (m_wet - mf)
#             if abs(mnew - m_wet) < 0.05:
#                 break
#             m_wet = 0.6*m_wet + 0.4*mnew

#         t_s = (m_wet * TUG_VEX / TUG_THR) * (1 - np.exp(-dv_e / TUG_VEX))

#         tug_dv[k]    = dv_e
#         tug_time[k]  = t_s / DAY
#         tug_mprop[k] = m_wet - m_pl
#         tug_mwet[k]  = m_wet

#     return tug_dv, tug_time, tug_mprop, tug_mwet

# trajectory_optimizer\nominal_mission.py
# def evaluate_sequence(seq, tug_dry=TUG_DRY, tug_mprop=TUG_MPROP, tug_time=TUG_TIME):
#     nf = [RH_IDX] + list(seq)
#     nt = list(seq) + [RH_IDX]

#     tug_prop_uniform = float(tug_mprop[list(seq)].max())
#     tug_mwets = np.full(5, tug_dry + tug_prop_uniform)

# LINE 229
# tug_arrive = np.array([tug_starts[i] + tug_time[seq[i]] for i in range(5)])

# if __name__ == '__main__':
#     import pandas as pd
#     from reaver_core import compute_tug_spirals

#     FORCED_IDS = [804, 884, 697, 488, 443]

#     id_map = {row[0]: i for i, row in enumerate(_RAW)}
#     debris = [id_map[x] for x in FORCED_IDS]
#     seq = optimise_ordering(debris)

#     tug_dry_mass_range = range(130, 270, 1)
#     rows = []

#     for tug_dry in tug_dry_mass_range:
#         print(f"=========== TUG_DRY: {tug_dry} kg ===========")

#         tug_dv, tug_time, tug_mprop, tug_mwet = compute_tug_spirals(
#             rh_raan_deg=RH_RAAN,
#             tug_dry=tug_dry
#         )

#         r = evaluate_sequence(
#             seq,
#             tug_dry=tug_dry,
#             tug_mprop=tug_mprop,
#             tug_time=tug_time
#         )

#         rows.append({
#         "tug_dry_mass_kg": tug_dry,
#         "mission_days": r["mission_day"]
#         })

#     final_table = pd.DataFrame(rows)

#     final_table.to_csv(
#         "mission_days_for_tug_dry_mass_range.csv",
#         sep=",",
#         index=False,
#         decimal="."
#     )

#     print("Saved mission days table to mission_days_for_tug_dry_mass_range.csv")
#     print("\n  Done.")