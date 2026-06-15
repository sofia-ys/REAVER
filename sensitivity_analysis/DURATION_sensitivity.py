import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

parameters = {}  # initialising our dictionary to store all the varions

'''VARYING DELTA V'''
# this is already done in trajectory_optimizer\\reaver_optimizer.py so we just rip the values from there
dv_table = pd.read_csv("sensitivity_analysis\\mission_days_and_delta_v_table.csv")
dv_table = dv_table.to_numpy()
dv_mission_days = dv_table[:, 0]

parameters["$\\Delta V$"] = dv_mission_days


'''VARYING MS DRY MASS'''
# CBE is 698.472, system-level contingency is 955.6597
ms_m_dry_table = pd.read_csv("sensitivity_analysis\\mission_days_for_ms_dry_mass_range.csv")
ms_m_dry_table = ms_m_dry_table.to_numpy()
ms_m_dry_mission_days = ms_m_dry_table[:, 1]

parameters["$m_{dry,ms}$"] = ms_m_dry_mission_days

'''VARYING TUG DRY MASS'''
tug_m_dry_table = pd.read_csv("sensitivity_analysis\\mission_days_for_tug_dry_mass_range.csv")
tug_m_dry_table = tug_m_dry_table.to_numpy()
tug_m_dry_mission_days = tug_m_dry_table[:, 1]

parameters["$m_{dry,tug}$"] = tug_m_dry_mission_days


'''maybe some margin on tug/ms wet mass'''

'''VARYING ISP'''
# CHEMICAL PROPELLANT MANUFACTURER RANGE: 243 - 256
ms_isp_table = pd.read_csv("sensitivity_analysis\\mission_days_for_ms_isp_range.csv")
ms_isp_table = ms_isp_table.to_numpy()
ms_isp_mission_days = ms_isp_table[:, 1]

parameters["$I_{sp,ms}$"] = ms_isp_mission_days

# ELECTRIC PROPELLANT MANUFACTURER RANGE: 1400 - 1650
tug_isp_table = pd.read_csv("sensitivity_analysis\\mission_days_for_tug_isp_range.csv")
tug_isp_table = tug_isp_table.to_numpy()
tug_isp_mission_days = tug_isp_table[:, 1]

parameters["$I_{sp,tug}$"] = tug_isp_mission_days


'''VARYING DEBRIS MASS'''
# every combination for debris 1000-2500 kg in increments of 200kg
m_debris_table = pd.read_csv("sensitivity_analysis\\mission_days_for_debris_mass_range_full.csv")
m_debris_table = m_debris_table.to_numpy()
m_debris_mission_days = m_debris_table[:, -1]

parameters["$m_{debris}$"] = m_debris_mission_days

'''VARYING POWER CONSUMPTION'''
power = 300


'''PLOTTING'''
plt.axvline(x=297.3, linestyle="--", label="Design value")
plt.axvline(x=365.0, linestyle="--", color="red", label="Constraint value")
plt.boxplot(parameters.values(), tick_labels=parameters.keys(), orientation='horizontal')
plt.xlabel(f"Mission duration [days]")
plt.ylabel(f"Technical parameter")
plt.legend(loc='upper right')
plt.show()