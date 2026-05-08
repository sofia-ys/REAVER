"""
Spacecraft Digital Twin — Dry Mass Sizing
==========================================
Scope: dry mass subsystems only (excludes propellant, tanks, and
       propellant-dependent structures, which scale with wet mass).

Mass roll-up hierarchy:
    Component → Subsystem → Spacecraft

Each subsystem exposes:
    - Configuration parameters (geometry, power, data rate, …)
    - .mass()           → float  [kg], computed from parameters
    - .mass_breakdown() → dict   {item: kg}
    - .contingency      → float  (fractional margin, e.g. 0.20 = 20 %)

Spacecraft.mass() sums subsystem masses and applies a system-level margin.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

@dataclass
class MassItem:
    """A single named mass contribution."""
    name: str
    mass_kg: float
    note: str = ""


def _apply_contingency(base_mass: float, contingency: float) -> float:
    """Return mass including contingency margin."""
    return base_mass * (1.0 + contingency)


# ---------------------------------------------------------------------------
# Subsystem base class
# ---------------------------------------------------------------------------

class Subsystem:
    """
    Abstract base for all subsystems.

    Attributes
    ----------
    contingency : float
        Fractional mass contingency (default 0.20 = 20 %).
        Applied on top of the computed base mass.
    heritage : str
        Design heritage note (e.g. 'TRL 5 - new design', 'flight proven').
        Informational only; may influence contingency choice.
    """

    DEFAULT_CONTINGENCY = 0.20

    def __init__(self, contingency: Optional[float] = None, heritage: str = ""):
        self.contingency = contingency if contingency is not None else self.DEFAULT_CONTINGENCY
        self.heritage = heritage

    # --- override in subclasses -------------------------------------------

    def _base_mass_items(self) -> list:
        """Return list of MassItems before contingency."""
        raise NotImplementedError

    # --- public API (do not override) -------------------------------------

    def mass(self) -> float:
        """Total subsystem mass including contingency [kg]."""
        base = sum(item.mass_kg for item in self._base_mass_items())
        return _apply_contingency(base, self.contingency)

    def mass_breakdown(self) -> dict:
        """
        Itemised mass breakdown.

        Returns dict with keys:
            'items'       - {name: base_mass_kg}
            'base_mass'   - float
            'contingency' - float (fraction)
            'total_mass'  - float (with contingency)
        """
        items = self._base_mass_items()
        base = sum(i.mass_kg for i in items)
        return {
            "items": {i.name: round(i.mass_kg, 4) for i in items},
            "base_mass": round(base, 4),
            "contingency": self.contingency,
            "total_mass": round(_apply_contingency(base, self.contingency), 4),
        }


# ---------------------------------------------------------------------------
# 1. Propulsion  (dry components only — no propellant, no tanks)
# ---------------------------------------------------------------------------

class Propulsion(Subsystem):
    """
    Dry propulsion hardware.

    Scope: thrusters, valves, feed-line components, and attitude-control
    thrusters.  Tanks and propellant are excluded (wet-mass domain).

    Parameters
    ----------
    n_main_thrusters       : number of main / delta-V thrusters
    mass_per_main_thruster : kg per unit (default: 0.5 kg, e.g. 22 N biprop)
    n_rcs_thrusters        : number of reaction-control / attitude thrusters
    mass_per_rcs_thruster  : kg per unit (default: 0.15 kg, monoprop 1 N)
    feed_system_mass       : kg — valves, regulators, lines (excl. tanks)
    """

    def __init__(
        self,
        n_main_thrusters: int = 1,
        mass_per_main_thruster: float = 0.5,
        n_rcs_thrusters: int = 4,
        mass_per_rcs_thruster: float = 0.15,
        feed_system_mass: float = 1.0,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.n_main_thrusters = n_main_thrusters
        self.mass_per_main_thruster = mass_per_main_thruster
        self.n_rcs_thrusters = n_rcs_thrusters
        self.mass_per_rcs_thruster = mass_per_rcs_thruster
        self.feed_system_mass = feed_system_mass

    def _base_mass_items(self) -> list:
        return [
            MassItem("main_thrusters",
                     self.n_main_thrusters * self.mass_per_main_thruster),
            MassItem("rcs_thrusters",
                     self.n_rcs_thrusters * self.mass_per_rcs_thruster),
            MassItem("feed_system", self.feed_system_mass,
                     note="valves, regulators, lines (no tanks)"),
        ]


# ---------------------------------------------------------------------------
# 2. Capture / Docking / Berthing mechanism
# ---------------------------------------------------------------------------

class Capture(Subsystem):
    """
    Capture, docking, or berthing mechanism.

    Parameters
    ----------
    mechanism_type   : 'passive' | 'active' | 'robotic_arm'
    mechanism_mass   : kg — structural/mechanical part of the interface
    actuator_mass    : kg — motors, latches (0 for passive)
    electronics_mass : kg — proximity sensors, cameras, controllers
    """

    MECHANISM_DEFAULTS = {
        "passive":     dict(mechanism_mass=2.0, actuator_mass=0.0,  electronics_mass=0.3),
        "active":      dict(mechanism_mass=3.5, actuator_mass=1.5,  electronics_mass=0.8),
        "robotic_arm": dict(mechanism_mass=8.0, actuator_mass=4.0,  electronics_mass=2.0),
    }

    def __init__(
        self,
        mechanism_type: str = "passive",
        mechanism_mass: Optional[float] = None,
        actuator_mass: Optional[float] = None,
        electronics_mass: Optional[float] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        defaults = self.MECHANISM_DEFAULTS.get(mechanism_type, self.MECHANISM_DEFAULTS["passive"])
        self.mechanism_type   = mechanism_type
        self.mechanism_mass   = mechanism_mass   if mechanism_mass   is not None else defaults["mechanism_mass"]
        self.actuator_mass    = actuator_mass    if actuator_mass    is not None else defaults["actuator_mass"]
        self.electronics_mass = electronics_mass if electronics_mass is not None else defaults["electronics_mass"]

    def _base_mass_items(self) -> list:
        return [
            MassItem("mechanism",   self.mechanism_mass,   note=self.mechanism_type),
            MassItem("actuators",   self.actuator_mass),
            MassItem("electronics", self.electronics_mass),
        ]


# ---------------------------------------------------------------------------
# 3. AOCS — Attitude and Orbit Control System
# ---------------------------------------------------------------------------

class AOCS(Subsystem):
    """
    Attitude and Orbit Control System hardware (sensors + actuators).

    Actuators
    ---------
    n_reaction_wheels : count  (3-axis or pyramid config)
    mass_per_rw       : kg
    n_magnetorquers   : count  (optional, LEO only)
    mass_per_mtq      : kg

    Sensors
    -------
    n_star_trackers   : count
    mass_per_st       : kg
    n_sun_sensors     : count
    mass_per_ss       : kg
    n_imu             : count  (IMU / gyro units)
    mass_per_imu      : kg
    gps_receiver      : bool
    gps_mass          : kg
    aocs_computer_mass: kg (0 if integrated into CDH)
    """

    def __init__(
        self,
        n_reaction_wheels: int = 3,
        mass_per_rw: float = 0.9,
        n_magnetorquers: int = 0,
        mass_per_mtq: float = 0.15,
        n_star_trackers: int = 1,
        mass_per_st: float = 0.35,
        n_sun_sensors: int = 4,
        mass_per_ss: float = 0.02,
        n_imu: int = 1,
        mass_per_imu: float = 0.25,
        gps_receiver: bool = True,
        gps_mass: float = 0.3,
        aocs_computer_mass: float = 0.0,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.n_reaction_wheels  = n_reaction_wheels
        self.mass_per_rw        = mass_per_rw
        self.n_magnetorquers    = n_magnetorquers
        self.mass_per_mtq       = mass_per_mtq
        self.n_star_trackers    = n_star_trackers
        self.mass_per_st        = mass_per_st
        self.n_sun_sensors      = n_sun_sensors
        self.mass_per_ss        = mass_per_ss
        self.n_imu              = n_imu
        self.mass_per_imu       = mass_per_imu
        self.gps_receiver       = gps_receiver
        self.gps_mass           = gps_mass
        self.aocs_computer_mass = aocs_computer_mass

    def _base_mass_items(self) -> list:
        items = [
            MassItem("reaction_wheels", self.n_reaction_wheels * self.mass_per_rw),
            MassItem("magnetorquers",   self.n_magnetorquers   * self.mass_per_mtq),
            MassItem("star_trackers",   self.n_star_trackers   * self.mass_per_st),
            MassItem("sun_sensors",     self.n_sun_sensors     * self.mass_per_ss),
            MassItem("imu",             self.n_imu             * self.mass_per_imu),
            MassItem("aocs_computer",   self.aocs_computer_mass),
        ]
        if self.gps_receiver:
            items.append(MassItem("gps_receiver", self.gps_mass))
        return items


# ---------------------------------------------------------------------------
# 4. EPS — Electrical Power System
# ---------------------------------------------------------------------------

class EPS(Subsystem):
    """
    Electrical Power System.

    Solar arrays sized from required power and specific power.
    Batteries sized from energy capacity and depth-of-discharge.

    Parameters
    ----------
    solar_array_power_w       : required BOL array output [W]
    specific_power_w_kg       : array specific power [W/kg]
                                (rigid crystalline Si ~80, GaAs rigid ~100,
                                 flexible ~150)
    array_structure_frac      : extra structure/harness as fraction of panel mass
    battery_capacity_wh       : total battery energy capacity [Wh]
    battery_specific_e_whkg   : cell specific energy [Wh/kg] (Li-ion ~150)
    battery_dod               : depth of discharge (sizing factor)
    n_battery_modules         : number of battery packs
    battery_module_overhead_kg: structural overhead per module [kg]
    pdcu_mass                 : kg — Power Distribution & Control Unit
    harness_mass              : kg — main power harness (EPS portion)
    """

    def __init__(
        self,
        solar_array_power_w: float = 500.0,
        specific_power_w_kg: float = 100.0,
        array_structure_frac: float = 0.20,
        battery_capacity_wh: float = 200.0,
        battery_specific_e_whkg: float = 150.0,
        battery_dod: float = 0.80,
        n_battery_modules: int = 1,
        battery_module_overhead_kg: float = 0.3,
        pdcu_mass: float = 2.5,
        harness_mass: float = 3.0,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.solar_array_power_w        = solar_array_power_w
        self.specific_power_w_kg        = specific_power_w_kg
        self.array_structure_frac       = array_structure_frac
        self.battery_capacity_wh        = battery_capacity_wh
        self.battery_specific_e_whkg    = battery_specific_e_whkg
        self.battery_dod                = battery_dod
        self.n_battery_modules          = n_battery_modules
        self.battery_module_overhead_kg = battery_module_overhead_kg
        self.pdcu_mass                  = pdcu_mass
        self.harness_mass               = harness_mass

    @property
    def solar_array_mass(self) -> float:
        base = self.solar_array_power_w / self.specific_power_w_kg
        return base * (1.0 + self.array_structure_frac)

    @property
    def battery_mass(self) -> float:
        cell_mass = (self.battery_capacity_wh / self.battery_dod) / self.battery_specific_e_whkg
        return cell_mass + self.n_battery_modules * self.battery_module_overhead_kg

    def _base_mass_items(self) -> list:
        return [
            MassItem("solar_arrays", self.solar_array_mass,
                     note=f"{self.solar_array_power_w} W @ {self.specific_power_w_kg} W/kg"),
            MassItem("batteries",    self.battery_mass,
                     note=f"{self.battery_capacity_wh} Wh, DoD={self.battery_dod}"),
            MassItem("pdcu",         self.pdcu_mass),
            MassItem("harness",      self.harness_mass),
        ]


# ---------------------------------------------------------------------------
# 5. Structures
# ---------------------------------------------------------------------------

class Structures(Subsystem):
    """
    Primary and secondary structure (dry, non-propellant-dependent).

    Sizing approach: structural mass is estimated as a fraction of the
    expected total dry mass (structural mass fraction method, common in
    early-phase sizing).  Alternatively, individual panel/frame masses
    may be specified directly.

    Parameters
    ----------
    primary_structure_mass   : kg — main load-bearing structure
    secondary_structure_mass : kg — brackets, panels, closeouts
    separation_system_mass   : kg — launch adapter ring / separation nuts
    mechanisms_mass          : kg — deployable hinges, hold-down releases
    """

    def __init__(
        self,
        primary_structure_mass: float = 10.0,
        secondary_structure_mass: float = 3.0,
        separation_system_mass: float = 1.5,
        mechanisms_mass: float = 1.0,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.primary_structure_mass   = primary_structure_mass
        self.secondary_structure_mass = secondary_structure_mass
        self.separation_system_mass   = separation_system_mass
        self.mechanisms_mass          = mechanisms_mass

    def _base_mass_items(self) -> list:
        return [
            MassItem("primary_structure",   self.primary_structure_mass),
            MassItem("secondary_structure", self.secondary_structure_mass),
            MassItem("separation_system",   self.separation_system_mass),
            MassItem("mechanisms",          self.mechanisms_mass),
        ]


# ---------------------------------------------------------------------------
# 6. TCS — Thermal Control System
# ---------------------------------------------------------------------------

class TCS(Subsystem):
    """
    Thermal Control System.

    Parameters
    ----------
    radiator_area_m2          : m2 — total radiator panel area
    radiator_areal_mass_kg_m2 : kg/m2 (Al radiator ~2-5 kg/m2)
    mli_blanket_mass          : kg — multi-layer insulation blankets
    heater_mass               : kg — electrical heaters + thermostats
    heat_pipe_mass            : kg — heat pipes / loop heat pipes
    phase_change_mass         : kg — phase-change material units (if used)
    """

    def __init__(
        self,
        radiator_area_m2: float = 0.5,
        radiator_areal_mass_kg_m2: float = 3.0,
        mli_blanket_mass: float = 1.0,
        heater_mass: float = 0.5,
        heat_pipe_mass: float = 0.8,
        phase_change_mass: float = 0.0,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.radiator_area_m2          = radiator_area_m2
        self.radiator_areal_mass_kg_m2 = radiator_areal_mass_kg_m2
        self.mli_blanket_mass          = mli_blanket_mass
        self.heater_mass               = heater_mass
        self.heat_pipe_mass            = heat_pipe_mass
        self.phase_change_mass         = phase_change_mass

    @property
    def radiator_mass(self) -> float:
        return self.radiator_area_m2 * self.radiator_areal_mass_kg_m2

    def _base_mass_items(self) -> list:
        return [
            MassItem("radiators",    self.radiator_mass,
                     note=f"{self.radiator_area_m2} m2"),
            MassItem("mli_blankets", self.mli_blanket_mass),
            MassItem("heaters",      self.heater_mass),
            MassItem("heat_pipes",   self.heat_pipe_mass),
            MassItem("pcm_units",    self.phase_change_mass),
        ]


# ---------------------------------------------------------------------------
# 7. TTC — Telemetry, Tracking & Command
# ---------------------------------------------------------------------------

class TTC(Subsystem):
    """
    Telemetry, Tracking & Command subsystem.

    Parameters
    ----------
    frequency_band    : 'UHF' | 'S' | 'X' | 'Ka'
    n_antennas        : number of antennas
    antenna_mass_each : kg per antenna (defaults by band if None)
    transceiver_mass  : kg — RF transceiver / transponder
    amplifier_mass    : kg — SSPA or TWTA
    rf_harness_mass   : kg — coax, waveguide, connectors
    """

    BAND_ANTENNA_DEFAULTS = {
        "UHF": 0.20,
        "S":   0.35,
        "X":   0.80,
        "Ka":  1.50,
    }

    def __init__(
        self,
        frequency_band: str = "S",
        n_antennas: int = 2,
        antenna_mass_each: Optional[float] = None,
        transceiver_mass: float = 1.2,
        amplifier_mass: float = 0.5,
        rf_harness_mass: float = 0.4,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.frequency_band    = frequency_band
        self.n_antennas        = n_antennas
        self.antenna_mass_each = (
            antenna_mass_each
            if antenna_mass_each is not None
            else self.BAND_ANTENNA_DEFAULTS.get(frequency_band, 0.35)
        )
        self.transceiver_mass  = transceiver_mass
        self.amplifier_mass    = amplifier_mass
        self.rf_harness_mass   = rf_harness_mass

    def _base_mass_items(self) -> list:
        return [
            MassItem("antennas",    self.n_antennas * self.antenna_mass_each,
                     note=f"{self.n_antennas}x {self.frequency_band}-band"),
            MassItem("transceiver", self.transceiver_mass),
            MassItem("amplifier",   self.amplifier_mass),
            MassItem("rf_harness",  self.rf_harness_mass),
        ]


# ---------------------------------------------------------------------------
# 8. CDH — Command & Data Handling
# ---------------------------------------------------------------------------

class CDH(Subsystem):
    """
    Command & Data Handling / On-Board Computer subsystem.

    Parameters
    ----------
    obc_mass           : kg — main on-board computer
    n_remote_terminals : number of RTUs / smart I/O cards
    mass_per_rtu       : kg per RTU
    mass_storage_mass  : kg — solid-state mass memory
    data_harness_mass  : kg — data bus cables (1553, SpaceWire, CAN, ...)
    """

    def __init__(
        self,
        obc_mass: float = 0.8,
        n_remote_terminals: int = 2,
        mass_per_rtu: float = 0.25,
        mass_storage_mass: float = 0.15,
        data_harness_mass: float = 0.6,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.obc_mass           = obc_mass
        self.n_remote_terminals = n_remote_terminals
        self.mass_per_rtu       = mass_per_rtu
        self.mass_storage_mass  = mass_storage_mass
        self.data_harness_mass  = data_harness_mass

    def _base_mass_items(self) -> list:
        return [
            MassItem("obc",              self.obc_mass),
            MassItem("remote_terminals", self.n_remote_terminals * self.mass_per_rtu),
            MassItem("mass_storage",     self.mass_storage_mass),
            MassItem("data_harness",     self.data_harness_mass),
        ]


# ---------------------------------------------------------------------------
# Spacecraft — top-level assembly
# ---------------------------------------------------------------------------

class SpaceCraft:
    """
    Spacecraft digital twin — dry mass sizing.

    Subsystems
    ----------
    propulsion : Propulsion  — dry propulsion hardware (no tanks/propellant)
    capture    : Capture     — docking / berthing mechanism
    aocs       : AOCS        — attitude & orbit control sensors/actuators
    eps        : EPS         — electrical power system
    structures : Structures  — primary & secondary structure
    tcs        : TCS         — thermal control
    ttc        : TTC         — telemetry, tracking & command
    cdh        : CDH         — command & data handling

    Parameters
    ----------
    system_margin : fractional system-level margin applied on top of all
                    subsystem masses (default 0.05 = 5 %).
                    This is separate from per-subsystem contingency.
    """

    SUBSYSTEM_NAMES = [
        "propulsion", "capture", "aocs", "eps",
        "structures", "tcs", "ttc", "cdh",
    ]

    def __init__(self, system_margin: float = 0.05):
        self.system_margin = system_margin

        self.propulsion: Optional[Propulsion] = None
        self.capture:    Optional[Capture]    = None
        self.aocs:       Optional[AOCS]       = None
        self.eps:        Optional[EPS]        = None
        self.structures: Optional[Structures] = None
        self.tcs:        Optional[TCS]        = None
        self.ttc:        Optional[TTC]        = None
        self.cdh:        Optional[CDH]        = None

    # ------------------------------------------------------------------
    # Mass roll-up
    # ------------------------------------------------------------------

    def _active_subsystems(self) -> dict:
        """Return {name: subsystem} for all non-None subsystems."""
        return {
            name: getattr(self, name)
            for name in self.SUBSYSTEM_NAMES
            if getattr(self, name) is not None
        }

    def mass(self) -> float:
        """
        Total spacecraft dry mass including system margin [kg].

        system_margin is applied once on top of the sum of all subsystem
        masses (each of which already includes their own contingency).
        """
        subsystem_total = sum(s.mass() for s in self._active_subsystems().values())
        return _apply_contingency(subsystem_total, self.system_margin)

    def mass_breakdown(self, verbose: bool = False) -> dict:
        """
        Full mass breakdown at subsystem (and optionally component) level.

        Parameters
        ----------
        verbose : if True, include per-item breakdown within each subsystem.

        Returns dict with structure:
            {
              'subsystems': {
                  'propulsion': {
                      'total_mass': float,
                      'base_mass':  float,
                      'contingency': float,
                      ['items': {...}]   # only if verbose=True
                  },
                  ...
              },
              'subsystem_total':  float,
              'system_margin':    float,
              'system_margin_kg': float,
              'total_dry_mass':   float,
            }
        """
        subsystem_data = {}
        subsystem_total = 0.0

        for name, subsys in self._active_subsystems().items():
            bd = subsys.mass_breakdown()
            entry = {
                "total_mass":  bd["total_mass"],
                "base_mass":   bd["base_mass"],
                "contingency": bd["contingency"],
            }
            if verbose:
                entry["items"] = bd["items"]
            subsystem_data[name] = entry
            subsystem_total += bd["total_mass"]

        total_dry = _apply_contingency(subsystem_total, self.system_margin)

        return {
            "subsystems":       subsystem_data,
            "subsystem_total":  round(subsystem_total, 3),
            "system_margin":    self.system_margin,
            "system_margin_kg": round(total_dry - subsystem_total, 3),
            "total_dry_mass":   round(total_dry, 3),
        }

    def summary(self, verbose: bool = False) -> str:
        """Human-readable mass budget table."""
        bd = self.mass_breakdown(verbose=verbose)
        lines = [
            "=" * 52,
            "  SPACECRAFT DRY MASS BUDGET",
            "=" * 52,
            f"  {'Subsystem':<22} {'Base':>7}  {'Cont.':>6}  {'Total':>7}",
            f"  {'-'*22} {'-'*7}  {'-'*6}  {'-'*7}",
        ]
        for name, data in bd["subsystems"].items():
            lines.append(
                f"  {name:<22} {data['base_mass']:>7.2f}  "
                f"{data['contingency']*100:>5.0f}%  {data['total_mass']:>7.2f}"
            )
            if verbose and "items" in data:
                for item_name, item_mass in data["items"].items():
                    if item_mass > 0:
                        lines.append(f"    {'> ' + item_name:<22} {item_mass:>7.2f}")

        lines += [
            f"  {'-'*22} {'-'*7}  {'-'*6}  {'-'*7}",
            f"  {'Subsystem total':<22} {'':>7}  {'':>6}  {bd['subsystem_total']:>7.2f}",
            f"  {'System margin':<22} {'':>7}  {bd['system_margin']*100:>5.0f}%  "
            f"{bd['system_margin_kg']:>7.2f}",
            "=" * 52,
            f"  {'TOTAL DRY MASS':<22} {'':>7}  {'':>6}  {bd['total_dry_mass']:>7.2f}  kg",
            "=" * 52,
        ]
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Quick-start example
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    sc = SpaceCraft(system_margin=0.05)

    sc.propulsion = Propulsion(
        n_main_thrusters=1,
        n_rcs_thrusters=8,
        feed_system_mass=1.2,
        contingency=0.15,
    )
    sc.capture = Capture(mechanism_type="active", contingency=0.20)
    sc.aocs = AOCS(
        n_reaction_wheels=4,   # pyramid config
        n_star_trackers=2,
        gps_receiver=True,
        contingency=0.15,
    )
    sc.eps = EPS(
        solar_array_power_w=800,
        specific_power_w_kg=120,
        battery_capacity_wh=300,
        contingency=0.10,
    )
    sc.structures = Structures(
        primary_structure_mass=12.0,
        secondary_structure_mass=4.0,
        contingency=0.10,
    )
    sc.tcs = TCS(
        radiator_area_m2=0.8,
        mli_blanket_mass=1.5,
        contingency=0.15,
    )
    sc.ttc = TTC(frequency_band="X", n_antennas=2, contingency=0.15)
    sc.cdh = CDH(obc_mass=1.2, n_remote_terminals=3, contingency=0.10)

    print(sc.summary(verbose=True))

