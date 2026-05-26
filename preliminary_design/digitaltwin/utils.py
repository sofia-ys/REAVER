from dataclasses import dataclass


@dataclass
class Component:
    """A single named spacecraft component"""
    name: str
    mass_kg: float
    note: str = ""