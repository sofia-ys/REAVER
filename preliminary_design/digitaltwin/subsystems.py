from typing import Optional

import numpy as np
from utils import Component

def _apply_contingency(base_mass: float, contingency: float) -> float:
    """Return mass including contingency margin."""
    return base_mass * (1.0 + contingency)

class Subsystem():

    def __init__(self, contingency: Optional[float] = None) -> None:
        self.DEFAULT_CONTINGENCY = 0.2 
        self.contingency = contingency if contingency is not None else self.DEFAULT_CONTINGENCY
        return

    def mass(self) -> float:
        """Total subsystem mass -- not including yet contingency [kg]."""
        base = sum(item.mass_kg for item in self._base_mass_items())
        return base # _apply_contingency(base, self.contingency)


#==================================================================
# MASS DEPENDENT SUBSYSTEMS
#==================================================================

# example mass budgets: page 82 table 20, page 83 table 22
class Propulsion(Subsystem):
    def __init__(self, m_dry:float, n_targets:int, dv_list:list, m_debris_list:list, Isp:float) -> None:
        super().__init__(contingency=0.1)

        self.m_dry = m_dry  # initial guess
        self.n_targets = n_targets
        self.dv_list = dv_list
        self.m_debris_list = m_debris_list
        self.Isp = Isp

        self.prop_mass_multiTarget()  # actually get the values
    
    # rocket equation to find prop mass
    def propellant_m(self, dv, Isp, m_final):
        return m_final * (np.e**(1000*dv/(Isp * 9.80665)) - 1)
    
    # Dry mass of a chemical RCS system [p. 166 ADSEE-1222] -- based on prop mass
    def rcs_mass(self, m_prop):
        return 0.178 * m_prop + 7.69
    
    # prop mass for multi-target missions where the s/c has to carry all its prop mass
    def sequence_prop_mass(self, m_dry, n_targets, dv_list, m_debris_list):    
        m_prop = 0
        for i in range(n_targets + 1):
                m_final = m_dry + sum(m_debris_list[:(n_targets - i)]) + m_prop 
                m_prop += self.propellant_m(dv_list[-(i+1)], self.Isp, m_final)
        return m_prop
    
    # calculating propellant mass and dry mass of propulsion system 
    def prop_mass_multiTarget(self):  # TODO: change the name maybe? naming convention?
        m_prop = self.sequence_prop_mass(self.m_dry, self.n_targets, self.dv_list, self.m_debris_list)  # initialising our propellant mass estimates
        m_prop_prev = 0 
        while abs(m_prop - m_prop_prev) > 1:  # convergence condition
                m_prop_prev = m_prop  # setting the previous estimate so we can compare
                m_prop = self.sequence_prop_mass(self.m_dry + self.rcs_mass(m_prop), self.n_targets, self.dv_list, self.m_debris_list)
        m_rcs = self.rcs_mass(m_prop)

        # saving the values we want into our class
        self.m_prop = m_prop 
        self.m_rcs = m_rcs
    
    def _base_mass_items(self):
        return [
            # Component("propellant_mass", self.m_prop),
            Component("feed_system_mass", self.m_rcs,
                      note="dry mass of propulsion system, all inclusive (not component-wise)")
        ]


class Structures(Subsystem):
    # no formula offered in ADSEE -> can just use percentage of dry mass as estimation 
    # dry mass range 50-1750 kg, mass of the structures (and mechanisms) system is roughly about 20 25% of the dry mass with a relative standard error of 21.7% (p. 104)
    # BASELINE REPORT USES 25% hence that value is taken
    def __init__(self, m_dry):
        super().__init__()
        self.m_structures = 0.25 * m_dry

    def _base_mass_items(self):
        return [
            Component("structures_mass", self.m_structures,
                      note="historical data based estimation")
        ]


class TCS(Subsystem):
    '''Thermal Control System (TCS) -- preliminary estimate'''
    # taking 3% of dry mass being TCS (from baseline report)
    def __init__(self, m_dry):
        super().__init__(contingency=0.1)
        self.m_tcs = 0.03 * m_dry

    def _base_mass_items(self):
        return [
            Component("tcs_mass", self.m_tcs,
                      note="historical data based estimation")
        ]
    
    '''IF WE WANT TO MAKE IT MORE DETAILED (later)'''
    # def __init__(self, coating_area, mli_area, osr_area, heatPipe_length, radiator_area, heater_area, m_dry):
    #     # initial estimate of dry mass for preliminary mass fraction calculations
    #     self.m_dry = m_dry

    #     # all parameters that influence TCS sizing
    #     self.coating_area = coating_area
    #     self.mli_area = mli_area
    #     self.osr_area = osr_area
    #     self.heatPipe_length = heatPipe_length
    #     self.radiator_area = radiator_area
    #     self.heater_area = heater_area

    #     # ADSEE p. 124 values for thermal materials
    #     self.coating_mass = 0.24  # [kg/m2] Paints/coatings that modify the emissivity and/or absorptance of a surface
    #     self.mli_mass = 0.3  # [kg/m2] Multi-layer insulation (MLI); Multiple layers of thin foils
    #     self.osr_mass = 1  # [kg/m2] Second Surface Mirrors or Optical Surface Reflectors (OSR)
    #     self.heatPipe_mass = 0.33  # [kg/m] Heat pipes: To transport heat from surfaces of high temperature to surfaces with lower temperature
    #     # Typical mass for active, deployable radiators is in range [5-15 kg/m2]
    #     self.radiator_mass = 15  # [kg/m2] Radiators: active radiators uses heat pipes to transport heat from a hot spot to a cold spot where the heat can be radiated out into space
    #     self.heater_mass = 2  # [kg/m2] Heaters: To provide local heating for instance to prevent propellants from freezing or a drastic reduction in battery capacity
    #     return

    # def _base_mass_items(self):
    #     return [
    #         Component("coating_mass", self.coating_mass * self.coating_area,
    #                   note="check density of chosen paint"),
    #         Component("mli_mass", self.mli_mass * self.mli_area,
    #                   note="check density of chosen insulator"),
    #         Component("osr_mass", self.osr_mass * self.osr_area,
    #                   note="check density of chosen insulator"),
    #         Component("heatPipe_mass", self.heatPipe_mass * self.heatPipe_length,
    #                   note="check length of chosen heat pipes"),
    #         Component("radiator_mass", self.radiator_mass * self.radiator_area,
    #                   note="check"),
    #         Component("heater_mass", self.heater_mass * self.heater_area,
    #                   note="check")
    #     ]

    # def _preliminary_mass(self):
    #     return 0.05 * self.dry_mass

    # TCS is 2-5% of s/c dry mass (page 115)

#==================================================================
# MASS INDEPENDENT SUBSYSTEMS
#==================================================================

class CaptureSystem(Subsystem):
    """
    Capture system including Capture Mechanism and Rendezvous, Proximity Operations sensors (RPO)
    """
    # i just wanna make this work, we can fix this more later
    # Capture system should be independent of mass, i think?
    def __init__(self, capture_type: str) -> None:
        super().__init__()
        # based on the baseline report 30% for capture
        # using random AI numbers <3 
        capture_masses = {
            "robotic_manipulator": 200,
            "gecko_adhesive": 25,
            "clamp": 50,
            "payload_adapter": 60,
            "mev_mechanism": 80
        }
        self.m_captureSys = capture_masses[capture_type]

    def _base_mass_items(self):
        return [
            Component("capture_system_mass", self.m_captureSys, note="should probs do for dif capture types")
        ]

    
    
class AOCS(Subsystem):
    """
    Attitude and Orbit Control Subsystem (AOCS)

    """
    def __init__(self, sc_type:str="complex") -> None:
        super().__init__(contingency=0.1)        
        self.sc_type = sc_type

        # list of aocs sensors as per the ones chosen in the midterm report, manually adjust these
        if self.sc_type == "complex":
            self.sensors = {
                # component_name: (qty., approx total mass [kg])
                "star_tracker": (2, 4),
                "gyro_unit": (2, 9),
                "coarse_sun_sensors": (8, 1.7),
                "fine_sun_sensors": (4, 0.2),
                "lidar": (2, 30.6),
                "cameras": (4, 2.1),
                "gps_receiver": (1, 1.2),  # TODO: shouldn't we have more for redundancy?
            }
            self.actuators = {
                "reaction_wheels": (4, 19.4),
                "thrusters": (12, 4),
                "solar_array_devices": (2, 3.3)
                }
        else:
            self.sensors = {
                "star_tracker": (1, 0.27),  # TODO: wouldn't we always want at least 2 for redundancy?
                "gyro_unit": (1, 0.06),
                "sun_sensors": (6, 0.22),
                "cameras": (2, 0.65),
            }
            self.actuators = {
                "reaction_wheels": (4, 4.4),
                "thrusters": (8, 3.1),
                "solar_array_devices": (1, 1.5)
            }

        self.mass_sum()
    
    def mass_sum(self):
        self.m_aocs_sensors = 0
        self.m_aocs_actuators = 0

        for sensor in self.sensors.values():
            self.m_aocs_sensors += sensor[1]
        for actuator in self.actuators.values():
            self.m_aocs_actuators += actuator[1]

        return
    
    def _base_mass_items(self):
        return [
            Component("aocs_sensors", self.m_aocs_sensors, note=f"for {self.sc_type} spacecraft"),
            Component("aocs_actuators", self.m_aocs_actuators, note=f"for {self.sc_type} spacecraft")
        ]
    
class EPS(Subsystem):
    """
    Electrical Power System (EPS), includes power generation, storage and handling.
    Assumes PV system including batteries for default value of specific power
    """

    # using random AI numbers <3
    def __init__(self, mission_type):
        super().__init__(contingency=0.2)  # if full solar array sizing -- contingency = 10
        eps_masses = {
            "STR": 150,
            "MTR": 200,
            "PSW": 300,
            "M&T": 750,
            "ISC": 1000
        }

        self.m_eps = eps_masses[mission_type]


    def _base_mass_items(self):
        return [
            Component("eps_mass", self.m_eps, note="can definitely be improved, im lazy"),
        ]

    '''so im being very lazy but (p. 136) in adsee actually goes over this'''
    # EPS is 20-50% of s/c dry mass (page 131)
    # power source mass esimtation (page 166)

    # EPS is not necessarily mass dependant, but probably mass related
    # Could potentially use a statistical relation of power-drymass instead of required power (how will we get this)
    # def __init__(
    #     self,
    #     power_output: float, # [W] continuous required power output of the power system
    #     specific_power: float = 12, #[W/kg] Photo-voltaic system (incl. batteries) range 7-12 We/kg
    #     # battery_capacity_kwh: float,
    #     # solar_power_w: float,
    #     # solar_cell_efficiency: float,
    #     #TODO: EXPAND
    #     ) -> None:
    #     self.power_output = power_output
    #     self.specific_power = specific_power
    #     return
    
    # def power_output_mass(self):
    #     specific_mass = 1 / self.specific_power
    #     return self.power_output * specific_mass


class TTC(Subsystem):
    def __init__(self, m_dry):
        super().__init__(contingency=0.1)
        # 2% from BASELINE REPORT
        self.m_ttc = m_dry * 0.02

    def _base_mass_items(self):
        return [
            Component("ttc_mass", self.m_ttc, note="can definitely be improved, im lazy"),
        ]
    
    ''' so the complicated explanation is on (p. 220)
        but i can't be bothered.'''
    # def __init__(self, tx_RF_power, tx_density, RF_power_req):
    #     self.tx_RF_power = tx_RF_power
    #     self.tx_density = tx_density
    #     self.RF_power_req = RF_power_req
    #     return
    
    # def mass_transponder(self):
    #     # low RF powers (up to about 10 W RF power)
    #     if self.RF_power_req == "low":
    #         return self.tx_RF_power / self.tx_density
    #     else:
    #         return 

# harness (wiring, cables, etc) is 3-10% of dry mass (page 143)
class CDH(Subsystem):
    """
    Command & Data Handling (CDH), consisting of
    - on-board computer (obc): standard obc
    - payload processor (pp): processing of vbn and robotics
    - remote terminal unit (rtu),
    - solid-state recorder,
    - data harness
    """

    # TODO: Consider COTS components, probably additional computing capabilities required for machine learning models used by visual-based navigation, possibly GPUs
    # TODO: payload processing might be better handled in the payload subsystem
    # TODO: maybe use different redundancy values for specific components
    # For redundancy, I think core spacecraft systems should be flight-proven radiation hardened components (so not single-board new-space systems)
    # but they provide insufficient computing capabilites for complex ML/AI vision-based navigation [1]
    # thus, separate (and probably redundant) high-performance gpus should be included, but I think these should
    # be included in the VBN system, thus in the payload subsystem.

    # [1] https://arc.aiaa.org/doi/pdf/10.2514/1.I010555

    # newspace (all in one): https://www.beyondgravity.com/sites/default/files/media_document/2026-02/cOBC_fact_sheet_2026-01-27.pdf

    def __init__(self, n_redundancy: int = 2):
        super().__init__(contingency=0.1)

        self.data_harness_mass = 0  # TODO: According to [Brown] the harness mass is in range 3-10% of on-orbit dry mass of the spacecraft. (p. 143)
        self.obc_mass   = 5.4   # [kg] https://www.beyondgravity.com/sites/default/files/media_document/2026-02/cOBC_fact_sheet_2026-01-27.pdf
        self.pp_mass    = 3.6   # [kg] https://www.beyondgravity.com/sites/default/files/media_document/2026-02/Satellites_FoX_Payload_Processor_Datasheet.pdf
        self.ssr_mass   = 0.75  # [kg] https://www.satnow.com/search/solid-state-recorders/filters?page=1&country=global&sorbit=;GEO;
        self.rtu_mass   = 17    # [kg] https://www.beyondgravity.com/sites/default/files/media_document/2023-11/Remote-Terminal-Unit.PDF
        self.n_redundancy = n_redundancy
        return

    def _base_mass_items(self):
        return [
            Component("data_harness_mass", self.data_harness_mass), #TODO should data harness scale with redundancy?
            Component("obc_mass", self.n_redundancy * self.obc_mass),
            Component("pp_mass", self.n_redundancy * self.pp_mass),
            Component("ssr_mass", self.n_redundancy * self.ssr_mass),
            Component("rtu_mass", self.n_redundancy * self.rtu_mass),
        ]