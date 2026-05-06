import numpy as np

class DeltaV:
    def _innit_(self, r_initial, r_final, i_initial, i_final):
        # Earth parameters
        mu_Earth = 398600.4415  # [km^2 / s^2]
        r_Earth = 6378.137  # [km]
        self.v_initial = np.sqrt(mu_Earth/r_initial)
        self.v_final = np.sqrt(mu_Earth/r_final)
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