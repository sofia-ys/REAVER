import numpy as np
from itertools import product

# VARYING DELTA V
delta_v = 1200

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
