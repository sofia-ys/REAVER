class Subsystem():
    def __init__(self) -> None:
        return

    def mass(self):
        """Calculate the mass of the subsystem"""
        return
    
# example mass budgets: page 82 table 20, page 83 table 22

class Propulsion(Subsystem):
    def __init__(
            self,
            n_main_thrusters: int,
            feed_system_mass: float, #TODO: maybe this is propellant dependant?
            #TODO: EXPAND
    ) -> None:
        return


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

class Structures(Subsystem):
    def __init__(self):
        return

class TCS(Subsystem):
    def __init__(self, paint_area, mli_area, osr_area, heatPipe_length, radiator_area, heater_area):
        # ADSEE p. 124
        self.paint_mass = 0.24 * paint_area  # Paints/coatings that modify the emissivity and/or absorptance of a surface
        self.mli_mass = 0.3 * mli_area  # Multi-layer insulation (MLI); Multiple layers of thin foils
        self.osr_mass = 1 * osr_area  # Second Surface Mirrors or Optical Surface Reflectors (OSR)
        self.heatPipe_mass = 0.33 * heatPipe_length  # Heat pipes: To transport heat from surfaces of high temperature to surfaces with lower temperature
        # Typical mass for active, deployable radiators is in range [5-15 kg/m2]
        self.radiator_mass = 15 * radiator_area  # Radiators: active radiators uses heat pipes to transport heat from a hot spot to a cold spot where the heat can be radiated out into space 
        self.heater_mass = 2 * heater_area  # Heaters: To provide local heating for instance to prevent propellants from freezing or a drastic reduction in battery capacity
        return
    
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
