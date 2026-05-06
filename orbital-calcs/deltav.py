import numpy as np

'''HOW TO USE:
inputs: initial and final radius and inclination
1)  create an instance of the class:
    deltaV = DeltaV(r_initial, r_final, i_initial, i_final)
2)  run the combined deltav manouevre function (does altitude and inclinatino change)
    delta_v = deltaV.combined()
3)  yay you now have a maybe wrong value for delta v 
'''

class DeltaV:
    def __init__(self, r_initial, r_final, i_initial, i_final):
        # Earth parameters
        self.mu_Earth = 398600.4415  # [km^2 / s^2]
        self.r_Earth = 6378.137  # [km]
        self.v_initial = np.sqrt(self.mu_Earth/r_initial)
        self.v_final = np.sqrt(self.mu_Earth/r_final)
        self.delta_i = i_final - i_initial
    
    # inclination and altitude change
    def combined(self):
        deltaV = self.v_initial**2 + self.v_final**2 - 2 * self.v_initial * self.v_final * np.cos(self.delta_i)
        return deltaV
    
    def phase(self):
        return "hello"
    
# recycling hub parameters
r_rh = 42164  # GEO radius [km]
i_rh = 7 * np.pi/180  # inclination [rad]