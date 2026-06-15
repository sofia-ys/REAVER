"""
REAVER RPO Budget Helpers
=========================
Tsiolkovsky conversions shared by the guidance and control layers, and a small
``PhaseResult`` container so every RPO phase reports a uniform
(ΔV, propellant, duration) triple.
"""
from dataclasses import dataclass
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
import numpy as np

from config import MS_VEX


def prop_from_dv(dv, m0, vex=MS_VEX):
    """Propellant burned for an impulsive ΔV at wet mass ``m0`` (Tsiolkovsky)."""
    return float(m0 * (1.0 - np.exp(-dv / vex)))


def dv_from_prop(m_prop, m0, vex=MS_VEX):
    """Equivalent ΔV for burning ``m_prop`` of propellant at wet mass ``m0``."""
    if m_prop >= m0:
        raise ValueError("propellant exceeds vehicle mass")
    return float(vex * np.log(m0 / (m0 - m_prop)))


@dataclass
class PhaseResult:
    """Uniform output of one RPO phase."""
    name: str
    dv: float          # m/s  (true ΔV for translation; equivalent ΔV for attitude)
    prop: float        # kg
    duration: float    # s
    note: str = ""

    def __str__(self):
        return (f"{self.name:<34s} {self.dv:8.3f} m/s   "
                f"{self.prop:7.4f} kg   {self.duration/60:7.2f} min")
