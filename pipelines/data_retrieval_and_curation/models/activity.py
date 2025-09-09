from dataclasses import dataclass

@dataclass
class Activity:
    compound_id: str
    target: str
    value: float
    units: str
    source: str
