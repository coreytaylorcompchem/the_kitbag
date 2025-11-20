import os
from pathlib import Path
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# OPI imports
from opi.core import Calculator
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import Task, Dft

# Manual mappings
DISPERSION_MAP = {
    "d3": "D3ZERO",
    "d3bj": "D3BJ",
    "d4": "D4",
}

RI_MAP = {
    "ri": "RI",
    "rijcosx": "RIJCOSX",
}

GRID_MAP = {
    "grid2": "Grid2",
    "grid3": "Grid3",
    "grid4": "Grid4",
    "grid5": "Grid5",
    "gridx5": "GridX5",
    "nofinalgrid": "NoFinalGrid",
}

REFERENCE_MAP = {
    "rhf": "RHF",
    "uhf": "UHF",
    "rohf": "ROHF",
}

# Map YAML scf_type to BlockScf scfmode
SCF_TYPES = {
    "pk": "direct",
    "rik": "rik",
    "gto": "gto",
    "gdm": "gdm",
}

class OrcaBackend:

    def __init__(self, global_config=None):
        self.global_config = global_config or {}
        self.working_dir = Path(
            self.global_config.get("output", {}).get("directory", ".")
        )

    def _make_calculator(self, xyz_file, basename, ncpu, ram):
        calc = Calculator(
            basename=basename,
            working_dir=str(self.working_dir)
        )
        calc.structure = Structure.from_xyz(xyz_file)
        calc.input.ncores = int(ncpu)
        calc.input.memory = int(ram)
        return calc

    def _apply_keywords(self, calc, step):
        # -------------------
        # DFT functional
        method = step.get("method", "").upper()
        if hasattr(Dft, method):
            func_keyword = getattr(Dft, method)
        else:
            raise ValueError(f"Unknown DFT functional '{method}'")

        # -------------------
        # Basis, Dispersion, Grid, RI
        basis = step.get("basis")
        disp_kw = DISPERSION_MAP.get(step.get("dispersion", "").lower())
        grid_kw = GRID_MAP.get(step.get("grid", "").lower())
        ri_kw = RI_MAP.get(step.get("ri", "").lower())

        # -------------------
        # Task: optimise / single point
        task = Task.OPT if step.get("task", "optimise") == "optimise" else Task.SP

        # -------------------
        # Build a single ORCA ! line
        # Only add each keyword once
        keywords_list = [func_keyword, basis, disp_kw, grid_kw, ri_kw, task]
        keywords_list = [kw for kw in keywords_list if kw]  # remove None
        def _keyword_to_str(kw):
            # If it's a SimpleKeyword, use its .keyword attribute
            if hasattr(kw, "keyword"):
                return kw.keyword
            return str(kw)
        
        calc.input.add_simple_keywords(" ".join(_keyword_to_str(kw) for kw in keywords_list))

        # -------------------
        # Maxcore / parallelization
        if "ram" in step:
            calc.input.maxcore = int(step["ram"])
        if "ncpu" in step:
            calc.input.nprocs = int(step["ncpu"])

        # -------------------
        # SCF block
        from opi.input.blocks.block_scf import BlockScf
        scf_kwargs = {}
        ref = step.get("reference")
        if ref:
            scf_kwargs["reference"] = REFERENCE_MAP[ref.lower()]
        if "maxiter" in step:
            scf_kwargs["maxiter"] = int(step["maxiter"])
        scf_type = step.get("scf_type")
        if scf_type:
            scf_kwargs["scfmode"] = SCF_TYPES.get(scf_type.lower())
            if scf_kwargs["scfmode"] is None:
                raise ValueError(f"Unknown SCF type '{scf_type}'")
        if scf_kwargs:
            calc.input.add_blocks(BlockScf(**scf_kwargs))

        # -------------------
        # Solvent block
        solv = step.get("solvent")
        if isinstance(solv, dict):
            model = solv.get("model", "").lower()
            eps = solv.get("epsilon", 80.0)
            if model == "cpcm":
                from opi.input.blocks.block_cpcm import BlockCpcm
                calc.input.add_blocks(BlockCpcm(epsilon=float(eps)))
            elif model == "cosmo":
                from opi.input.blocks.block_cosmors import BlockCosmors
                calc.input.add_blocks(BlockCosmors(epsilon=float(eps)))
            elif model == "smd":
                calc.input.add_simple_keywords("smd")

    # ---------------------------------------------------------
    def optimise(self, xyz_file, step, log_path=None):
        calc = self._make_calculator(
            xyz_file,
            basename="optimise",
            ncpu=step.get("ncpu", 4),
            ram=step.get("ram", 4000),
        )

        # DO NOT add Task.OPT here, handled inside _apply_keywords
        self._apply_keywords(calc, step)

        calc.write_input()
        return calc.run()

    # ---------------------------------------------------------
    def single_point(self, xyz_file, step, log_path=None):
        calc = self._make_calculator(
            xyz_file,
            basename="single_point",
            ncpu=step.get("ncpu", 4),
            ram=step.get("ram", 4000),
        )

        # DO NOT add Task.SP here, handled inside _apply_keywords
        self._apply_keywords(calc, step)

        calc.write_input()
        return calc.run()

