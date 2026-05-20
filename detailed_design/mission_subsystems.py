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
        return base # _apply_contingency(base, self.contingency)  # 

class Propulsion(Subsystem):
    '''Propulsion, includes tanks'''
    def __init__(self, m_dry, sc_type, m_wet_t=None) -> None:
        super().__init__(contingency=0.1)
        self.sc_type = sc_type
        self.m_dry = m_dry  # [kg] initial guess for dry mass

        self.m_debris = 1800  # [kg] mass of the collected debris, assuming all are same (for now)
        self.power = None

        if sc_type == "ms":
            self.Isp_ms = 253  # [s] chemical propulsion for the mothership
            self.dv_ms = [0.306, 0.306, 0.306, 0.306, 0.306, 0.306]  # [km/s] average case dv grand tour
            self.m_wet_t = m_wet_t  # we take a fixed value for the tug wet mass so we can use it, but we don't iterate it
            self.m_prop_ms, self.m_rcs_ms = self.find_m_prop(sc_type="ms")  # [kg] prop mass and prop sys dry mass for the mothership
        else:  # values we use for the tug
            self.Isp_t = 1660  # [s] electric propulsion for the tug
            self.dv_t = 0.52  # [km/s] average case dv debris -> RH for each tug
            self.m_prop_t, self.m_rcs_t = self.find_m_prop(sc_type="tug")  # [kg] prop mass and prop sys dry mass for the tug
            thrust = self.calculate_thrust(burn_duration_days=230) #230 days for tug
            self.power = self.calculate_power(thrust)

    def calculate_thrust(self, burn_duration_days):
        """
        Calculate the thrust required based off of propellant mass and burn duration
        m_dot = m_prop / t_b
        v_e = Isp * g0
        dp = v_e*dm
        F = dp/dt = m_dot * v_e
        """
        t_b = burn_duration_days * 24 * 60 * 60  # [s]
        if not self.m_prop_t:
            raise ValueError('Propellant mass is not defined')

        return (self.m_prop_t / t_b) * self.Isp_t * 9.80665

    def calculate_power(self, thrust, prop_efficiency=0.7):
        """
        Calculate the electric power required by propulsion system, using propulsion system efficiency.
        Default value of efficiency is 0.7 based off of NEXT-C thruster
        """
        self.power = (thrust * self.Isp_t * 9.80665) / (2 * prop_efficiency)
        return self.power

    def propellant_m(self, dv, Isp, m_final):
        return m_final * (np.e**(1000*dv/(Isp * 9.80665)) - 1)
    
    def rcs_mass(self, m_prop):
        return 0.178 * m_prop + 7.69

    def tug_m_prop(self, m_dry_t):
        m_final = m_dry_t + self.m_debris
        tug_m_prop = self.propellant_m(self.dv_t, self.Isp_t, m_final)
        return tug_m_prop
    
    '''FOR FOUR TUGS (MS GOES BACK WITH DEBRIS)'''
    def ms_m_prop(self, m_dry_ms):
        m_prop = 0
        for i in range(len(self.dv_ms)):  # number of manouevres
            m_final = m_dry_ms + (i - 1) * self.m_wet_t + m_prop
            if i == 0:
                m_final += self.m_debris + self.m_wet_t  # for the last manouevre, the mothership also carries the last debris
            m_prop += self.propellant_m(self.dv_ms[-(i+1)], self.Isp_ms, m_final)
        return m_prop
    
    '''FOR FIVE TUGS (MS GOES BACK ALONE)'''
    # def ms_m_prop(self, m_dry_ms):
    #     m_prop = 0
    #     for i in range(len(self.dv_ms)):  # number of manouevres
    #         m_final = m_dry_ms + i * self.m_wet_t + m_prop
    #         m_prop += self.propellant_m(self.dv_ms[-(i+1)], self.Isp_ms, m_final)
    #     return m_prop
    
    def find_m_prop(self, sc_type):
        m_prop = self.ms_m_prop(self.m_dry) if sc_type=="ms" else self.tug_m_prop(self.m_dry)  # initialising our propellant mass estimates
        m_prop_prev = 0 
        while abs(m_prop - m_prop_prev) > 1:  # convergence condition
            m_prop_prev = m_prop  # setting the previous estimate so we can compare
            m_prop = self.ms_m_prop(self.m_dry + self.rcs_mass(m_prop)) if sc_type=="ms" else self.tug_m_prop(self.m_dry + self.rcs_mass(m_prop))
        m_rcs = self.rcs_mass(m_prop)
        return m_prop, m_rcs

    def _base_mass_items(self):
        if self.sc_type == "ms":
            return [Component("mothership_rcs_mass", self.m_rcs_ms)]
        else:  # for the tug
            return [Component("tug_rcs_mass", self.m_rcs_t)]
    
    def propellant_mass(self):
        if self.sc_type == "ms":
            return [Component("mothership_propellant_mass", self.m_prop_ms)]
        else:
            return [Component("tug_propellant_mass", self.m_prop_t)]

class AOCS(Subsystem):
    '''Attitude and Orbit Control System (AOCS)'''
    def __init__(self, sc_type:str) -> None:
        super().__init__(contingency=0.1)        
        self.sc_type = sc_type  # either mothership or tug

        # list of aocs sensors as per the ones chosen in the midterm report, manually adjust these
        if self.sc_type == "ms":
            self.sensors = {
                # component_name: (qty., approx total mass [kg], approx power [W])
                "star_tracker": (2, 4, 12),
                "gyro_unit": (2, 9, 27),
                "coarse_sun_sensors": (8, 1.7, 0),
                "fine_sun_sensors": (4, 0.2, 0),
                "lidar": (2, 30.6, 160),
                "cameras": (4, 2.1, 8),
                "gps_receiver": (1, 1.2, 8),  # TODO: shouldn't we have more for redundancy?
            }
            self.actuators = {
                "reaction_wheels": (4, 19.4, 80), #TODO: Check, says 80 ss in report, steady-state?
                "thrusters": (12, 4, 0), #TODO: POWER TBD
                "solar_array_devices": (2, 3.3, 0) #TODO: POWER TBD
            }
        else:  # sc_type == "tug"
            self.sensors = {
                "star_tracker": (1, 0.27, 1.5),  # TODO: wouldn't we always want at least 2 for redundancy?
                "gyro_unit": (1, 0.06, 1.5),
                "sun_sensors": (6, 0.22, 0.9),
                "cameras": (2, 0.65, 3.6),
            }
            self.actuators = {
                "reaction_wheels": (4, 4.4, 56),  # range 16-56, chose upper limit
                "thrusters": (8, 3.1, 0),  # TODO: POWER TBD
                "solar_array_devices": (1, 1.5, 0)  # TODO: POWER TBD
            }

        self.mass_sum()
        self.power_sum()
    
    def mass_sum(self):
        self.m_aocs_sensors = 0
        self.m_aocs_actuators = 0

        for sensor in self.sensors.values():
            self.m_aocs_sensors += sensor[1]
        for actuator in self.actuators.values():
            self.m_aocs_actuators += actuator[1]

        return
    def power_sum(self):
        self.power_aocs_sensors = 0
        self.power_aocs_actuators = 0

        for sensor in self.sensors.values():
            self.power_aocs_sensors += sensor[2]
        for actuator in self.actuators.values():
            self.power_aocs_actuators += actuator[2]

        return

    
    def _base_mass_items(self):
        return [
            Component("aocs_sensors", self.m_aocs_sensors, note=f"for {self.sc_type} spacecraft"),
            Component("aocs_actuators", self.m_aocs_actuators, note=f"for {self.sc_type} spacecraft")
        ]

class CaptureSystem(Subsystem):
    '''Capture system including Capture Mechanism, and Rendezvous & Proximity Operations (RPO) sensors'''
    def __init__(self, sc_type) -> None:
        super().__init__()

        # using random AI numbers <3 
        capture_masses = {
            "robotic_manipulator": 200,
            "mev_mechanism": 80
        }

        if sc_type == "ms":
            self.m_captureSys = capture_masses["robotic_manipulator"]
        else:  # for the tug 
            self.m_captureSys = capture_masses["mev_mechanism"]

    def _base_mass_items(self):
        return [
            Component("capture_system_mass", self.m_captureSys, note="get a better mass estimate pls")
        ]
    
class EPS(Subsystem):
    '''Electrical Power System (EPS), includes power generation, storage and handling'''
    def __init__(self, sc_type, prop_power: Optional[float] = None):
        super().__init__(contingency=0.2)  # if full solar array sizing -- contingency = 10
        
        if sc_type == "ms":
            self.P_eps = 2000 # W, from baseline report (probably needs to be revised)
            # self.m_eps = 120  # [kg]
        else:  # for the tug, since it uses electric propulsion
            # self.m_eps = 40  # [kg]
            self.P_eps = prop_power + 300 # Just prop power for now with 300W for everything else

        self.m_eps = self.solar_array_sizing() + self.battery_sizing() + self.pmd_sizing()
        self.m_eps = self.m_eps + self.harness_sizing()

    def solar_array_sizing(self):
        """ Solar array sizing (page 132 from ADSEE AE1222 reader pdf)"""
        lifetime = 1 # Lifetime is one year
        eff_d = 0.8
        degradation_factor = 0.01
        p_d = self.P_eps

        p_sa_bol = p_d / eff_d  # Eclipse (or nighttime efficiency is ignored for now) P_e / eff_e * t_e/t_d
        p_sa_eol = p_sa_bol * (1 - degradation_factor)**lifetime
        p_sp_sa = 120 # W/kg
        p_delta_sa = 5 # m^2/kg
        m_sa = p_sa_eol / p_sp_sa
        a_sa = p_sa_eol / p_delta_sa
        return m_sa

    def battery_sizing(self):
        """ Assume no battery for now, because no ecplise. Although it is probably relevant for peak power"""
        m_bat = 0
        return m_bat

    def pmd_sizing(self):
        m_pmd = 0.071 * self.P_eps * 0.15 # kg
        return m_pmd

    def harness_sizing(self):
        # m_har = 0.03 * self.m_eps # 3-10% of on-orbit dry mass (ADSEE reader)
        m_har = 0.25 * self.m_eps # 3-10% of on-orbit dry mass (ADSEE reader)
        return m_har

    def _base_mass_items(self):
        return [
            Component("eps_mass", self.m_eps, note="can definitely be improved, im lazy"),
        ]
    
class CDH(Subsystem):
    '''Command & Data Handling (CDH), includes
        - on-board computer (obc): standard obc
        - payload processor (pp): processing of vbn and robotics
        - remote terminal unit (rtu)
        - solid-state recorder
        - harness is included in EPS
        '''

    def __init__(self, sc_type):
        super().__init__(contingency=0.1)

        if sc_type == "ms":
            self.obc_mass   = 5.4   # [kg] dual-redundant https://www.beyondgravity.com/sites/default/files/media_document/2026-02/cOBC_fact_sheet_2026-01-27.pdf
            self.pp_mass    = 7.6   # [kg] dual-redundant https://www.beyondgravity.com/sites/default/files/media_document/2026-02/Satellites_FoX_Payload_Processor_Datasheet.pdf
            self.ssr_mass   = 0.75 * 2 # [kg] two units https://www.satnow.com/search/solid-state-recorders/filters?page=1&country=global&sorbit=;GEO;
            self.rtu_mass   = 17.0    # [kg] redundant https://www.beyondgravity.com/sites/default/files/media_document/2023-11/Remote-Terminal-Unit.PDF
            self.n_redundancy = 1  # most components have redundancy built-in
            self.m_cdh = self.n_redundancy * (self.obc_mass + self.pp_mass + self.ssr_mass) + self.rtu_mass
        else:
            self.m_cdh = 5.4 # kg already dual-redundant, https://www.beyondgravity.com/sites/default/files/media_document/2026-02/cOBC_fact_sheet_2026-01-27.pdf
        return

    def _base_mass_items(self):
        return [
            Component("m_cdh", self.m_cdh)
        ]

###########################################
# NOT YET REVISED FOR MOTHERSHIP AND TUGS #
###########################################

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

class TTC(Subsystem):
    def __init__(self, m_dry):
        super().__init__(contingency=0.1)
        self.m_ttc = m_dry * 0.02

    def _base_mass_items(self):
        return [
            Component("ttc_mass", self.m_ttc, note="can definitely be improved, im lazy"),
        ]