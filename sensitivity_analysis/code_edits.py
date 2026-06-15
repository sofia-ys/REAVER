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


'''VARYING MS ISP'''
# trajectory_optimizer\nominal_mission.py
# def evaluate_sequence(seq, ms_isp=MS_ISP):
#     ms_vex = ms_isp * G0

#     nf = [RH_IDX] + list(seq)
#     nt = list(seq) + [RH_IDX]

#     tug_prop_uniform = float(TUG_MPROP[list(seq)].max())
#     tug_mwets = np.full(5, TUG_DRY + tug_prop_uniform)

#     # Backward pass
#     m = MS_DRY
#     m = m * np.exp(DV_LEG[seq[-1], RH_IDX] / ms_vex)

#     for i in range(3, -1, -1):
#         m += tug_mwets[i+1]
#         m = m * np.exp(DV_LEG[seq[i], seq[i+1]] / ms_vex)

#     m += tug_mwets[0]
#     m = m * np.exp(DV_LEG[RH_IDX, seq[0]] / ms_vex)

#     ms_prop_orbital = m - MS_DRY - tug_mwets.sum()

#     rcs_allocation = MS_RCS_MARGIN * ms_prop_orbital
#     ms_prop = ms_prop_orbital + rcs_allocation
#     ms_wet0 = MS_DRY + tug_mwets.sum() + ms_prop

#     # Forward pass: mass + timing (orbital burns + RPO burns interleaved)
#     dv_legs = np.array([DV_LEG[nf[i], nt[i]] for i in range(6)])
#     t_legs  = np.array([T_LEG_FINITE[nf[i], nt[i]] for i in range(6)])

#     mass = ms_wet0
#     mass_vec = []
#     t = 0.0
#     tug_starts = []
#     for i in range(6):
#         mass = mass * np.exp(-dv_legs[i] / ms_vex)
#         mass_vec.append(mass)
#         t += t_legs[i]
#         if i < 5:
#             # RPO at debris: V-bar inspect + tumble-axis capture + detumble
#             mass = mass * np.exp(-(DV_RPO_DEBRIS + DV_RPO_DETUMBLE) / ms_vex)
#             t += T_OPS
#             tug_starts.append(t)
#             mass -= tug_mwets[i]
#     # 5 × RH proximity docking (tug-meet + dock) as each tug+debris arrives
#     for _ in range(5):
#         mass = mass * np.exp(-DV_RPO_RH / ms_vex)
#     ms_return = t

#     rcs_margin  = mass - MS_DRY   # remaining RCS propellant at end of mission
#     tug_arrive  = np.array([tug_starts[i] + TUG_TIME[seq[i]] for i in range(5)])

#     # Sequential RH handovers: each must wait for (a) the previous handover to finish
#     # and (b) the tug to arrive — no two handovers overlap.
#     arrive_order = np.argsort(tug_arrive)
#     handover = np.zeros(5)
#     t_rh = ms_return
#     for rank in range(5):
#         i = int(arrive_order[rank])
#         t_start = max(t_rh, float(tug_arrive[i]))
#         handover[i] = t_start + T_OPS
#         t_rh = handover[i]
#     mission_day = handover.max()
#     tot_dv      = dv_legs.sum()

#     return dict(
#         seq             = seq,
#         ms_isp = ms_isp,
#         ms_vex = ms_vex,
#         ms_prop_orbital = ms_prop_orbital,
#         rcs_allocation  = rcs_allocation,
#         ms_prop         = ms_prop,         # total initial prop (orbital + RCS)
#         ms_wet0         = ms_wet0,
#         tug_mwets       = tug_mwets,
#         dv_legs         = dv_legs,
#         t_legs          = t_legs,
#         mass_vec        = np.array(mass_vec),
#         tug_starts      = np.array(tug_starts),
#         tug_arrive      = tug_arrive,
#         handover        = handover,
#         ms_return       = ms_return,
#         mission_day     = mission_day,
#         tot_dv          = tot_dv,
#         feasible        = mission_day <= MAX_DAYS,
#         rcs_margin      = rcs_margin,
#     )

# if __name__ == '__main__':
#     import pandas as pd

#     print("="*70)
#     print("  REAVER — MS ISP SENSITIVITY")
#     print("="*70)

#     ms_isp_range = range(240, 260, 1)

#     FORCED_IDS = [804, 884, 697, 488, 443]

#     id_map = {row[0]: i for i, row in enumerate(_RAW)}
#     debris = [id_map[x] for x in FORCED_IDS]

#     # Optimise order once using nominal MS_ISP
#     seq = optimise_ordering(debris)

#     rows = []

#     for ms_isp in ms_isp_range:
#         print(f"=========== MS_ISP: {ms_isp} s ===========")

#         r = evaluate_sequence(seq, ms_isp=ms_isp)

#         rows.append({
#             "ms_isp_s": ms_isp,
#             "mission_days": r["mission_day"]
#         })

#     final_table = pd.DataFrame(rows)

#     final_table.to_csv(
#         "mission_days_for_ms_isp_range.csv",
#         sep=",",
#         index=False,
#         decimal="."
#     )

#     print("Saved MS Isp sensitivity table to mission_days_for_ms_isp_range.csv")
#     print("\n  Done.")


'''VARYING TUG ISP'''
# trajectory_optimizer\reaver_core.py
# def compute_tug_spirals(rh_raan_deg=RH_RAAN, tug_isp=TUG_ISP):
#     tug_vex = tug_isp * G0

#     tug_dv    = np.zeros(N_DEB)
#     tug_time  = np.zeros(N_DEB)
#     tug_mprop = np.zeros(N_DEB)
#     tug_mwet  = np.zeros(N_DEB)

#     v2 = np.sqrt(MU / RH_SMA)
#     i2, o2 = RH_INC, rh_raan_deg

#     for k in range(N_DEB):
#         m_pl  = TUG_DRY + MASS[k]
#         m_wet = m_pl * 1.35

#         i1, o1 = INC[k], RAAN[k]
#         v1 = np.sqrt(MU / SMA[k])

#         cos_d = (
#             np.cos(i1*D2R)*np.cos(i2*D2R) +
#             np.sin(i1*D2R)*np.sin(i2*D2R)*np.cos((o2-o1)*D2R)
#         )

#         dth  = np.arccos(np.clip(cos_d, -1.0, 1.0))
#         dv_e = np.sqrt(v1**2 + v2**2 - 2*v1*v2*np.cos(np.pi/2*dth))

#         for _ in range(60):
#             mf   = m_wet * np.exp(-dv_e / tug_vex)
#             mnew = m_pl + (m_wet - mf)
#             if abs(mnew - m_wet) < 0.05:
#                 break
#             m_wet = 0.6*m_wet + 0.4*mnew

#         t_s = (m_wet * tug_vex / TUG_THR) * (1 - np.exp(-dv_e / tug_vex))

#         tug_dv[k]    = dv_e
#         tug_time[k]  = t_s / DAY
#         tug_mprop[k] = m_wet - m_pl
#         tug_mwet[k]  = m_wet

#     return tug_dv, tug_time, tug_mprop, tug_mwet

# trajectory_optimizer\nominal_mission.py
# def evaluate_sequence(seq, tug_mprop=TUG_MPROP, tug_time=TUG_TIME):
#     nf = [RH_IDX] + list(seq)
#     nt = list(seq) + [RH_IDX]
#     tug_prop_uniform = float(tug_mprop[list(seq)].max())
#     tug_mwets = np.full(5, TUG_DRY + tug_prop_uniform)

#     # Backward pass: orbital prop sized so final mass = MS_DRY (transfers only)
#     m = MS_DRY
#     m = m * np.exp(DV_LEG[seq[-1], RH_IDX] / MS_VEX)
#     for i in range(3, -1, -1):
#         m += tug_mwets[i+1]
#         m  = m * np.exp(DV_LEG[seq[i], seq[i+1]] / MS_VEX)
#     m += tug_mwets[0]
#     m  = m * np.exp(DV_LEG[RH_IDX, seq[0]] / MS_VEX)
#     ms_prop_orbital = m - MS_DRY - tug_mwets.sum()

#     # 10 % RCS margin added at launch — covers all proximity operations (RPO)
#     rcs_allocation = MS_RCS_MARGIN * ms_prop_orbital
#     ms_prop        = ms_prop_orbital + rcs_allocation   # total initial propellant
#     ms_wet0        = MS_DRY + tug_mwets.sum() + ms_prop

#     # Forward pass: mass + timing, recalculated for current tug mass
#     dv_legs = np.array([DV_LEG[nf[i], nt[i]] for i in range(6)])
#     t_legs = np.zeros(6)

#     mass = ms_wet0
#     mass_vec = []
#     t = 0.0
#     tug_starts = []

#     for i in range(6):
#         fi = nf[i]
#         ti = nt[i]

#         d1  = DV1[fi, ti]
#         d2  = DV2[fi, ti]
#         dph = DV_PH[fi, ti]

#         leg_time = 0.0

#         # Burn 1
#         burn_time = finite_burn_time(d1, mass, SMA_ALL[fi])
#         t += burn_time
#         leg_time += burn_time
#         mass = mass * np.exp(-d1 / MS_VEX)

#         # Hohmann transfer coast
#         t += T_TR[fi, ti]
#         leg_time += T_TR[fi, ti]

#         # Burn 2
#         burn_time = finite_burn_time(d2, mass, SMA_ALL[ti])
#         t += burn_time
#         leg_time += burn_time
#         mass = mass * np.exp(-d2 / MS_VEX)

#         # Burn 3: phasing burn
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
#             # RPO at debris
#             mass = mass * np.exp(-(DV_RPO_DEBRIS + DV_RPO_DETUMBLE) / MS_VEX)

#             # Fixed operations time
#             t += T_OPS

#             # Tug starts after release
#             tug_starts.append(t)

#             # Tug released from mothership
#             mass -= tug_mwets[i]

#     # 5 × RH proximity docking (tug-meet + dock) as each tug+debris arrives
#     for _ in range(5):
#         mass = mass * np.exp(-DV_RPO_RH / MS_VEX)
#     ms_return = t

#     rcs_margin  = mass - MS_DRY   # remaining RCS propellant at end of mission
#     tug_arrive = np.array([tug_starts[i] + tug_time[seq[i]] for i in range(5)])

#     # Sequential RH handovers: each must wait for (a) the previous handover to finish
#     # and (b) the tug to arrive — no two handovers overlap.
#     arrive_order = np.argsort(tug_arrive)
#     handover = np.zeros(5)
#     t_rh = ms_return
#     for rank in range(5):
#         i = int(arrive_order[rank])
#         t_start = max(t_rh, float(tug_arrive[i]))
#         handover[i] = t_start + T_OPS
#         t_rh = handover[i]
#     mission_day = handover.max()
#     tot_dv      = dv_legs.sum()

#     return dict(
#         seq             = seq,
#         ms_prop_orbital = ms_prop_orbital,
#         rcs_allocation  = rcs_allocation,
#         ms_prop         = ms_prop,         # total initial prop (orbital + RCS)
#         ms_wet0         = ms_wet0,
#         tug_mwets       = tug_mwets,
#         dv_legs         = dv_legs,
#         t_legs          = t_legs,
#         mass_vec        = np.array(mass_vec),
#         tug_starts      = np.array(tug_starts),
#         tug_arrive      = tug_arrive,
#         handover        = handover,
#         ms_return       = ms_return,
#         mission_day     = mission_day,
#         tot_dv          = tot_dv,
#         feasible        = mission_day <= MAX_DAYS,
#         rcs_margin      = rcs_margin,
#     )

# if __name__ == '__main__':
#     import pandas as pd
#     from reaver_core import compute_tug_spirals

#     print("="*70)
#     print("  REAVER — TUG ISP SENSITIVITY")
#     print("="*70)

#     tug_isp_range = range(1400, 1651, 1)

#     FORCED_IDS = [804, 884, 697, 488, 443]

#     id_map = {row[0]: i for i, row in enumerate(_RAW)}
#     debris = [id_map[x] for x in FORCED_IDS]

#     # Optimise order once using nominal tug model
#     seq = optimise_ordering(debris)

#     rows = []

#     for tug_isp in tug_isp_range:
#         print(f"=========== TUG_ISP: {tug_isp} s ===========")

#         tug_dv, tug_time, tug_mprop, tug_mwet = compute_tug_spirals(
#             rh_raan_deg=RH_RAAN,
#             tug_isp=tug_isp
#         )

#         r = evaluate_sequence(
#             seq,
#             tug_mprop=tug_mprop,
#             tug_time=tug_time
#         )

#         rows.append({
#             "tug_isp_s": tug_isp,
#             "mission_days": r["mission_day"]
#         })

#     final_table = pd.DataFrame(rows)

#     final_table.to_csv(
#         "mission_days_for_tug_isp_range.csv",
#         sep=",",
#         index=False,
#         decimal="."
#     )

#     print("Saved tug Isp sensitivity table to mission_days_for_tug_isp_range.csv")
#     print("\n  Done.")


'''VARYING DEBRIS MASS'''
# trajectory_optimizer\reaver_core.py
# def compute_tug_spirals(rh_raan_deg=RH_RAAN, debris_mass=MASS):
#     tug_dv    = np.zeros(N_DEB)
#     tug_time  = np.zeros(N_DEB)
#     tug_mprop = np.zeros(N_DEB)
#     tug_mwet  = np.zeros(N_DEB)

#     v2 = np.sqrt(MU / RH_SMA)
#     i2, o2 = RH_INC, rh_raan_deg

#     for k in range(N_DEB):
#         m_pl  = TUG_DRY + debris_mass[k]
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

# # trajectory_optimizer\nominal_mission.py
# def evaluate_sequence(seq, tug_mprop=TUG_MPROP, tug_time=TUG_TIME):
#     nf = [RH_IDX] + list(seq)
#     nt = list(seq) + [RH_IDX]
#     tug_prop_uniform = float(tug_mprop[list(seq)].max())
#     tug_mwets = np.full(5, TUG_DRY + tug_prop_uniform)

#     # Backward pass: orbital prop sized so final mass = MS_DRY (transfers only)
#     m = MS_DRY
#     m = m * np.exp(DV_LEG[seq[-1], RH_IDX] / MS_VEX)
#     for i in range(3, -1, -1):
#         m += tug_mwets[i+1]
#         m  = m * np.exp(DV_LEG[seq[i], seq[i+1]] / MS_VEX)
#     m += tug_mwets[0]
#     m  = m * np.exp(DV_LEG[RH_IDX, seq[0]] / MS_VEX)
#     ms_prop_orbital = m - MS_DRY - tug_mwets.sum()

#     # 10 % RCS margin added at launch — covers all proximity operations (RPO)
#     rcs_allocation = MS_RCS_MARGIN * ms_prop_orbital
#     ms_prop        = ms_prop_orbital + rcs_allocation   # total initial propellant
#     ms_wet0        = MS_DRY + tug_mwets.sum() + ms_prop

#     # Forward pass: mass + timing (orbital burns + RPO burns interleaved)
#     dv_legs = np.array([DV_LEG[nf[i], nt[i]] for i in range(6)])
#     t_legs  = np.array([T_LEG_FINITE[nf[i], nt[i]] for i in range(6)])

#     mass = ms_wet0
#     mass_vec = []
#     t = 0.0
#     tug_starts = []
#     for i in range(6):
#         mass = mass * np.exp(-dv_legs[i] / MS_VEX)
#         mass_vec.append(mass)
#         t += t_legs[i]
#         if i < 5:
#             # RPO at debris: V-bar inspect + tumble-axis capture + detumble
#             mass = mass * np.exp(-(DV_RPO_DEBRIS + DV_RPO_DETUMBLE) / MS_VEX)
#             t += T_OPS
#             tug_starts.append(t)
#             mass -= tug_mwets[i]
#     # 5 × RH proximity docking (tug-meet + dock) as each tug+debris arrives
#     for _ in range(5):
#         mass = mass * np.exp(-DV_RPO_RH / MS_VEX)
#     ms_return = t

#     rcs_margin  = mass - MS_DRY   # remaining RCS propellant at end of mission
#     tug_arrive = np.array([tug_starts[i] + tug_time[seq[i]] for i in range(5)])

#     # Sequential RH handovers: each must wait for (a) the previous handover to finish
#     # and (b) the tug to arrive — no two handovers overlap.
#     arrive_order = np.argsort(tug_arrive)
#     handover = np.zeros(5)
#     t_rh = ms_return
#     for rank in range(5):
#         i = int(arrive_order[rank])
#         t_start = max(t_rh, float(tug_arrive[i]))
#         handover[i] = t_start + T_OPS
#         t_rh = handover[i]
#     mission_day = handover.max()
#     tot_dv      = dv_legs.sum()

#     return dict(
#         seq             = seq,
#         ms_prop_orbital = ms_prop_orbital,
#         rcs_allocation  = rcs_allocation,
#         ms_prop         = ms_prop,         # total initial prop (orbital + RCS)
#         ms_wet0         = ms_wet0,
#         tug_mwets       = tug_mwets,
#         dv_legs         = dv_legs,
#         t_legs          = t_legs,
#         mass_vec        = np.array(mass_vec),
#         tug_starts      = np.array(tug_starts),
#         tug_arrive      = tug_arrive,
#         handover        = handover,
#         ms_return       = ms_return,
#         mission_day     = mission_day,
#         tot_dv          = tot_dv,
#         feasible        = mission_day <= MAX_DAYS,
#         rcs_margin      = rcs_margin,
#     )

# if __name__ == '__main__':
#     import pandas as pd
#     from reaver_core import compute_tug_spirals
#     from itertools import product

#     print("="*70)
#     print("  REAVER — DEBRIS MASS SENSITIVITY")
#     print("="*70)


#     FORCED_IDS = [804, 884, 697, 488, 443]

#     id_map = {row[0]: i for i, row in enumerate(_RAW)}
#     debris = [id_map[x] for x in FORCED_IDS]

#     # Optimise order once using nominal catalogue masses
#     seq = optimise_ordering(debris)

#     m_debris = list(range(1000, 2501, 200))
#     m_debris += [2440, 2035, 2087, 1585.73, 1902.87]
#     m_debris_range = product(m_debris, repeat=5)

#     rows = []

#     for mass_combo in product(m_debris, repeat=5):
#         print(f"=========== debris masses: {mass_combo} kg ===========")

#         debris_mass_custom = MASS.copy()

#         # Assign custom masses to the selected debris, in visit order
#         for idx, new_mass in zip(seq, mass_combo):
#             debris_mass_custom[idx] = new_mass

#         # Recompute only tug spirals with changed debris masses
#         tug_dv, tug_time, tug_mprop, tug_mwet = compute_tug_spirals(
#             rh_raan_deg=RH_RAAN,
#             debris_mass=debris_mass_custom
#         )

#         r = evaluate_sequence(
#             seq,
#             tug_mprop=tug_mprop,
#             tug_time=tug_time
#         )

#         rows.append({
#             "m_debris_1_kg": mass_combo[0],
#             "m_debris_2_kg": mass_combo[1],
#             "m_debris_3_kg": mass_combo[2],
#             "m_debris_4_kg": mass_combo[3],
#             "m_debris_5_kg": mass_combo[4],
#             "total_debris_mass_kg": sum(mass_combo),
#             "mission_days": r["mission_day"]
#         })

#     final_table = pd.DataFrame(rows)

#     final_table.to_csv(
#         "mission_days_for_debris_mass_range.csv",
#         sep=",",
#         index=False,
#         decimal="."
#     )

#     print("Saved debris mass sensitivity table to mission_days_for_debris_mass_range.csv")
#     print("\n  Done.")