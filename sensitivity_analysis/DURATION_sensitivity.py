import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

parameters = {}  # initialising our dictionary to store all the varions

# VARYING DELTA V -- this is already done in trajectory_optimizer\\reaver_optimizer.py so we just rip the values from there
# the table to get the values we wanted is:
# import pandas as pd
# dv_table = pd.DataFrame({
#     "mission_days": res["mission_day"],
#     "total_delta_v_m_s": res["total_dv"]
# })
# dv_table.to_csv(
#     "mission_days_and_delta_v_table.csv",
#     sep=",",
#     index=False,
#     decimal="."
# )
# print("Saved mission days and delta-v table to mission_days_and_delta_v_table.csv")

dv_table = pd.read_csv("sensitivity_analysis\mission_days_and_delta_v_table.csv")
dv_table.to_csv(
    "mission_days_and_delta_v_table.csv",
    sep=",",
    index=False,
    decimal="."
)
dv_table = dv_table.to_numpy()
mission_days = dv_table[:, 0]
total_delta_v = dv_table[:, 1]

parameters["$\Delta V$"] = mission_days

# VARYING DRY MASS
m_dry = 600

# VARYING PROPELLANT MASS
m_prop = 2300

# VARYING WET MASS
m_wet = 2900

# VARYING ISP
chemical_isp_range = range(243, 256, 1)  
electric_isp_range = range(1400, 1650, 10)


# VARYING POWER CONSUMPTION
power = 300


# VARYING DEBRIS MASS
m_debris = range(1000, 2001, 100)
m_debris_range = product(m_debris, repeat=5)


# PLOTTING
plt.boxplot(parameters.values(), labels=parameters.keys(), orientation='horizontal')
plt.axvline(x=365, linestyle="--", color="red")
plt.xlabel(f"Mission duration [days]")
plt.ylabel(f"Technical parameter")
plt.show()