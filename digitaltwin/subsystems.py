import numpy as np
from utils import Component

class Subsystem():
    def __init__(self) -> None:
        return

    def mass(self):
        """Calculate the mass of the subsystem"""
        return
    
# example mass budgets: page 82 table 20, page 83 table 22

class Propulsion(Subsystem):
    def __init__(self, m_dry, n_targets, dv_list, m_debris_list, Isp):
        self.m_dry = m_dry
        self.n_targets = n_targets
        self.dv_list = dv_list
        self.m_debris_list = m_debris_list
        self.Isp = Isp
        self.m_prop, self.m_rcs = self.prop_mass_multiTarget()  # actually get the values
        return
    
    # rocket equation to find prop mass
    def propellant_m(self, dv, Isp, m_final):
        return m_final * (np.e**(1000*dv/(Isp * 9.80665)) - 1)
    
    # Dry mass of a chemical RCS system [p. 166 ADSEE-1222] -- based on prop mass
    def rcs_mass(self, m_prop):
        return 0.178 * m_prop + 7.69
    
    # prop mass for multi-target missions where the s/c has to carry all its prop mass
    def sequence_prop_mass(self):    
        m_prop = 0
        for i in range(self.n_targets + 1):
                m_final = self.m_dry + sum(self.m_debris_list[:(self.n_targets - i)]) + m_prop 
                m_prop += self.propellant_m(self.dv_list[-(i+1)], self.Isp, m_final)
        return m_prop
    
    # calculating propellant mass and dry mass of propulsion system 
    def prop_mass_multiTarget(self):  # TODO: change the name maybe? naming convention?
        m_prop, _ = self.sequence_prop_mass(self.m_dry, self.n_targets, self.dv_list, self.m_debris_list)  # initialising our propellant mass estimates
        m_prop_prev = 0 
        while abs(m_prop - m_prop_prev) > 1:  # convergence condition
                m_prop_prev = m_prop  # setting the previous estimate so we can compare
                m_prop = self.sequence_prop_mass(self.m_dry + self.rcs_mass(m_prop), self.n_targets, self.dv_list, self.m_debris_list)
        m_rcs = self.rcs_mass(m_prop)
        return m_prop, m_rcs
    
    def _base_mass_items(self):
        return [
            Component("propellant_mass", self.m_prop),
            Component("feed_system_mass", self.m_rcs,
                      note="dry mass of propulsion system, all inclusive (not component-wise)")
        ]



class CaptureSystem(Subsystem):
    """
    Capture system including Capture Mechanism and Rendezvous, Proximity Operations sensors (RPO)
    """
    def __init__(
            self
    ):
        return

class AOCS(Subsystem):
    """
    Attitude and Orbit Control Subsystem (AOCS)

    """
    def __init__(
            self,
            # Actuators
            n_reaction_wheels: int,
            n_mangetorquers: int,
            n_control_moment_gyros: int,
            n_rcs_thrusters: int, #TODO: is rcs included in propulsion or AOCS?

            # Sensors
            n_star_trackers: int,
            n_sun_sensors: int,
            n_horizon_sensors: int,
            n_gps_receivers: int,
            n_imus: int,
            n_magnetometers: int,

    ) -> None:
        return

class EPS(Subsystem):
    # EPS is 20-50% of s/c dry mass (page 131)
    # power source mass esimtation (page 166)
    def __init__(
        self,
        battery_capacity_kwh: float,
        solar_power_w: float,
        solar_cell_efficiency: float,
        #TODO: EXPAND
        ) -> None:
        return
    
    def power_output_mass(self):
        specific_power = 1/()
        return self.power_output * specific_power

class Structures(Subsystem):
    def __init__(self):
        return

class TCS(Subsystem):
    def __init__(self, coating_area, mli_area, osr_area, heatPipe_length, radiator_area, heater_area, m_dry):
        # initial estimate of dry mass for preliminary mass fraction calculations
        self.m_dry = m_dry

        # all parameters that influence TCS sizing
        self.coating_area = coating_area
        self.mli_area = mli_area
        self.osr_area = osr_area
        self.heatPipe_length = heatPipe_length
        self.radiator_area = radiator_area
        self.heater_area = heater_area

        # ADSEE p. 124 values for thermal materials
        self.coating_mass = 0.24  # [kg/m2] Paints/coatings that modify the emissivity and/or absorptance of a surface
        self.mli_mass = 0.3  # [kg/m2] Multi-layer insulation (MLI); Multiple layers of thin foils
        self.osr_mass = 1  # [kg/m2] Second Surface Mirrors or Optical Surface Reflectors (OSR)
        self.heatPipe_mass = 0.33  # [kg/m] Heat pipes: To transport heat from surfaces of high temperature to surfaces with lower temperature
        # Typical mass for active, deployable radiators is in range [5-15 kg/m2]
        self.radiator_mass = 15  # [kg/m2] Radiators: active radiators uses heat pipes to transport heat from a hot spot to a cold spot where the heat can be radiated out into space 
        self.heater_mass = 2  # [kg/m2] Heaters: To provide local heating for instance to prevent propellants from freezing or a drastic reduction in battery capacity
        return
    
    def _base_mass_items(self):
        return [
            Component("coating_mass", self.coating_mass * self.coating_area,
                      note="check density of chosen paint"),
            Component("mli_mass", self.mli_mass * self.mli_area,
                      note="check density of chosen insulator"),
            Component("osr_mass", self.osr_mass * self.osr_area,
                      note="check density of chosen insulator"),
            Component("heatPipe_mass", self.heatPipe_mass * self.heatPipe_length,
                      note="check length of chosen heat pipes"),
            Component("radiator_mass", self.radiator_mass * self.radiator_area,
                      note="check"),
            Component("heater_mass", self.heater_mass * self.heater_area,
                      note="check")
        ]
    
    def _preliminary_mass(self):
        return 0.05 * dry_mass
    
    # TCS is 2-5% of s/c dry mass (page 115)
    

class TTC(Subsystem):
    def __init__(self,
                frequency_band: str,
                n_antennas: int,):
        return

# harness (wiring, cables, etc) is 3-10% of dry mass (page 143)
class CDH(Subsystem):
    def __init__(
                self,
                harness_mass: float):
        return
