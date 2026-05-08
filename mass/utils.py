from dataclasses import dataclass


@dataclass
class MassItem:
    """A single named mass contribution."""
    name: str
    mass_kg: float
    note: str = ""