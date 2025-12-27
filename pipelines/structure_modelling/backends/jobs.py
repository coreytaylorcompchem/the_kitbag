from dataclasses import dataclass
from pathlib import Path

@dataclass
class InferenceJob:
    tool_name: str
    tool_backend: object
    run_id: int
    seed: int
    device: int
    output_dir: Path
    refinement: dict | None
