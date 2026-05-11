class Subsystem():
    def __init__(self) -> None:
        return

    def mass(self):
        """Calculate the mass of the subsystem"""
        return

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
    def __init__(
            self,
            battery_capacity_kwh: float,
            solar_power_w: float,
            solar_cell_efficiency: float,
            #TODO: EXPAND
    ) -> None:
        return

class Structures(Subsystem):

class TCS(Subsystem):

class TTC(Subsystem):

class CDH(Subsystem):