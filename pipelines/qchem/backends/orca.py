import os
from pathlib import Path
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# OPI imports
from opi.core import Calculator
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import Task, Dft

# SCF block imports
from opi.input.blocks.block_scf import BlockScf, DIIS, Shift


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

    # ---------------------------------------------------------
    # Normalization helper: ensures strings are clean and safe
    def _norm(self, kw):
        if kw is None:
            return None
        return str(kw).strip()

    # ---------------------------------------------------------
    def _apply_keywords(self, calc, step):
        # -------------------
        # DFT functional (only part still handled by OPI)
        method = step.get("method", "").upper()
        if hasattr(Dft, method):
            func_keyword = getattr(Dft, method)
        else:
            raise ValueError(f"Unknown DFT functional '{method}'")

        # -------------------
        # ORCA keywords come directly from YAML as valid ORCA keywords
        basis   = self._norm(step.get("basis"))
        disp_kw = self._norm(step.get("dispersion"))
        grid_kw = self._norm(step.get("grid"))      # use DEFGRID1/2/3
        ri_kw   = self._norm(step.get("ri"))

        # -------------------
        # Task keyword
        task = Task.OPT if step.get("task", "optimise") == "optimise" else Task.SP

        # -------------------
        # Build ORCA "! ..." line
        keyword_pieces = [func_keyword.keyword]
        for raw in (basis, disp_kw, grid_kw, ri_kw):
            if raw:
                keyword_pieces.append(raw)
        keyword_pieces.append(task.keyword)

        calc.input.add_simple_keywords(" ".join(keyword_pieces))

        # -------------------
        # Maxcore / parallelization
        if "ram" in step:
            calc.input.maxcore = int(step["ram"])
        if "ncpu" in step:
            calc.input.nprocs = int(step["ncpu"])

        # -------------------
        # SCF block
        scf_kwargs = {}

        # reference
        ref = step.get("reference")
        if ref:
            scf_kwargs["reference"] = self._norm(ref)

        # SCF options
        scf_opts = step.get("scf", {})
        for k, v in scf_opts.items():
            key = k.lower()
            if key == "diis":
                if isinstance(v, bool):
                    scf_kwargs["diis"] = DIIS(enabled=v)
                elif isinstance(v, dict):
                    scf_kwargs["diis"] = DIIS(**v)
                else:
                    raise ValueError(f"Invalid diis value: {v}")
            elif key == "shift":
                if isinstance(v, (float, int)):
                    scf_kwargs["shift"] = Shift(scf=v)
                elif isinstance(v, dict):
                    scf_kwargs["shift"] = Shift(**v)
                else:
                    raise ValueError(f"Invalid shift value: {v}")
            elif key == "scfmode":
                scf_kwargs["scfmode"] = v
            elif key == "maxiter":
                scf_kwargs["maxiter"] = int(v)
            else:
                # Pass other SCF keys as-is
                scf_kwargs[key] = v

        if scf_kwargs:
            calc.input.add_blocks(BlockScf(**scf_kwargs))

        # -------------------
        # Solvent block
        solv = step.get("solvent")
        if isinstance(solv, dict):
            model = solv.get("model", "").lower()
            eps = float(solv.get("epsilon", 80.0))
            solv_name = solv.get("name", "").upper() if solv.get("name") else None
            use_draco = bool(solv.get("draco", False))

            if model == "cpcm":
                from opi.input.blocks.block_cpcm import BlockCpcm

                kwargs = {"epsilon": eps}
                if solv_name:
                    kwargs["solvent"] = solv_name  # specify solvent by name
                if use_draco:
                    kwargs["draco"] = True  # enable DRACO

                calc.input.add_blocks(BlockCpcm(**kwargs))

            elif model == "smd":
                calc.input.add_simple_keywords("SMD")

    # ---------------------------------------------------------
    def optimise(self, xyz_file, step, log_path=None):
        calc = self._make_calculator(
            xyz_file,
            basename="optimise",
            ncpu=step.get("ncpu", 4),
            ram=step.get("ram", 4000),
        )
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
        self._apply_keywords(calc, step)
        calc.write_input()
        return calc.run()
