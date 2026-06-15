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


'''VARYING ISP'''
chemical_isp_range = range(243, 256, 1)  
electric_isp_range = range(1400, 1650, 10)


'''VARYING POWER CONSUMPTION'''
power = 300


'''VARYING DEBRIS MASS'''
m_debris = range(1000, 2001, 100)
m_debris_range = product(m_debris, repeat=5)


'''PLOTTING'''
plt.axvline(x=297, linestyle="--", label="Design value")
plt.axvline(x=365, linestyle="--", color="red", label="Constraint value")
plt.boxplot(parameters.values(), tick_labels=parameters.keys(), orientation='horizontal')
plt.xlabel(f"Mission duration [days]")
plt.ylabel(f"Technical parameter")
plt.legend()
plt.show()