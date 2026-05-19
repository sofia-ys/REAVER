from typing import Optional
from dataclasses import dataclass
import numpy as np

@dataclass
class Component:
    '''A single named spacecraft component'''
    name: str
    mass_kg: float
    note: str = ""

def _apply_contingency(base_mass: float, contingency: float) -> float:
    '''Return mass including a subsystem contingency margin'''
    return base_mass * (1.0 + contingency)

class Subsystem():
    '''Generic functions for each subsystem'''
    def __init__(self, contingency: Optional[float] = None) -> None:
        self.DEFAULT_CONTINGENCY = 0.2 
        self.contingency = contingency if contingency is not None else self.DEFAULT_CONTINGENCY
        return

    def mass(self) -> float:
        """Total subsystem mass -- not including yet contingency [kg]."""
        base = sum(item.mass_kg for item in self._base_mass_items())
        return _apply_contingency(base, self.contingency)  # m_subsystem 

class Propulsion(Subsystem):
    '''Propulsion, includes tanks'''
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
    '''Structures'''
    def __init__(self, m_dry):
        super().__init__()
        self.m_structures = 0.25 * m_dry

    def _base_mass_items(self):
        return [
            Component("structures_mass", self.m_structures,
                      note="historical data based estimation")
        ]


class TCS(Subsystem):
    '''Thermal Control System (TCS)'''
    def __init__(self, m_dry):
        super().__init__(contingency=0.1)
        self.m_tcs = 0.03 * m_dry

    def _base_mass_items(self):
        return [
            Component("tcs_mass", self.m_tcs,
                      note="historical data based estimation")
        ]
    

class CaptureSystem(Subsystem):
    '''Capture system including Capture Mechanism, and Rendezvous & Proximity Operations (RPO) sensors'''
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
    '''Attitude and Orbit Control System (AOCS)'''
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
    '''Electrical Power System (EPS), includes power generation, storage and handling'''
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


class TTC(Subsystem):
    def __init__(self, m_dry):
        super().__init__(contingency=0.1)
        self.m_ttc = m_dry * 0.02

    def _base_mass_items(self):
        return [
            Component("ttc_mass", self.m_ttc, note="can definitely be improved, im lazy"),
        ]

class CDH(Subsystem):
    '''Command & Data Handling (CDH), includes
        - on-board computer (obc): standard obc
        - payload processor (pp): processing of vbn and robotics
        - remote terminal unit (rtu)
        - solid-state recorder
        - data harness'''

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