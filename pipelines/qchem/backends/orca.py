import os

from pathlib import Path
import json
import numpy as np
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# Core OPI imports
from opi.core import Calculator
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import Task, Dft

# Block imports used to build the import file
from opi.input.blocks.block_scf import BlockScf, DIIS, Shift
from opi.input.blocks.block_cpcm import BlockCpcm


class OrcaBackend:
    def __init__(self, global_config=None):
        self.global_config = global_config or {}
        self.working_dir = Path(
            self.global_config.get("output", {}).get("directory", ".")
        )

    def _wrap_result(self, calc, output):

        prop_path = Path(calc.working_dir) / "optimise.property.json"
        with open(prop_path, "r") as f:
            data = json.load(f)

        geoms = data.get("Geometries", [])
        last = geoms[-1]

        # Electronic energy fallback
        energy = (
            last.get("DFT_Energy", {}).get("finalEn")
            or last["Energy"][0]["totalEnergy"][0][0]
        )

        # Gibbs free energy from frequency calculation
        gibbs_energy = last.get("ThermalData", {}).get("GibbsFreeEnergy", energy)

        # coords
        cart = last["Geometry"]["Coordinates"]["Cartesians"]
        atoms = [row[0] for row in cart]
        coords_bohr = np.array([row[1:] for row in cart], dtype=float)
        BOHR_TO_ANG = 0.52917721092
        coords = coords_bohr * BOHR_TO_ANG

        return {
            "energy": energy,
            "gibbs_energy": gibbs_energy,
            "atoms": atoms,
            "coords": coords.tolist(),
            "charge": data.get("Calculation_Info", {}).get("Charge"),
            "multiplicity": data.get("Calculation_Info", {}).get("Mult"),
        }

    def _make_calculator(self, xyz_file, basename, ncpu, ram, charge=0, mult=1):
        calc = Calculator(
            basename=basename,
            working_dir=str(self.working_dir)
        )

        structure = Structure.from_xyz(xyz_file)

        # Explicitly set charge/multiplicity
        structure.charge = int(charge)
        structure.mult = int(mult)

        calc.structure = structure

        calc.input.ncores = int(ncpu)
        calc.input.memory = int(ram)

        return calc


    # Normalisation helper: ensures strings are clean and safe
    def _norm(self, kw):
        if kw is None:
            return None
        return str(kw).strip()

    def _apply_keywords(self, calc, step):
        # DFT functional (the only part still handled by OPI directly)
        method = step.get("method", "").upper()
        if hasattr(Dft, method):
            func_keyword = getattr(Dft, method)
        else:
            raise ValueError(f"Unknown DFT functional '{method}'")

        # ORCA keywords that come directly from YAML as valid ORCA keywords
        basis   = self._norm(step.get("basis"))
        disp_kw = self._norm(step.get("dispersion"))
        grid_kw = self._norm(step.get("grid"))      # use DEFGRID1/2/3, not the old Grid4, GridX5, etc.
        ri_kw   = self._norm(step.get("ri"))

        # Task keyword
        task_name = step.get("task", "optimise")
        task = Task.OPT if task_name=="optimise" else Task.SP

        # Build ORCA "! ..." input line
        keyword_pieces = [func_keyword.keyword]
        for raw in (basis, disp_kw, grid_kw, ri_kw):
            if raw:
                keyword_pieces.append(raw)
        keyword_pieces.append(task.keyword)

        # Add frequency calculation if requested
        if step.get("freq", False) and task_name=="optimise":
            keyword_pieces.append("FREQ")

        calc.input.add_simple_keywords(" ".join(keyword_pieces))

        if "ram" in step:
            calc.input.maxcore = int(step["ram"])
        if "ncpu" in step:
            calc.input.nprocs = int(step["ncpu"])

        # SCF block
        scf_kwargs = {}

        # reference wave function
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

        # Solvent block (still haven't got COSMO working but PCM covers more solvents anyway)
        solv = step.get("solvent")
        if isinstance(solv, dict):
            model = solv.get("model", "").lower()
            eps = float(solv.get("epsilon", 80.0))
            solv_name = solv.get("name", "").upper() if solv.get("name") else None
            use_draco = bool(solv.get("draco", False))

            if model == "cpcm":

                kwargs = {"epsilon": eps}
                if solv_name:
                    kwargs["solvent"] = solv_name  # specify solvent by name
                if use_draco:
                    kwargs["draco"] = True  # enable DRACO - adds a lot of cost and not sure if it's worth it

                calc.input.add_blocks(BlockCpcm(**kwargs))

            elif model == "smd":
                calc.input.add_simple_keywords("SMD")

        logger.warning(f"BACKEND RECEIVED SOLVENT BLOCK: {step.get('solvent')}")

    def optimise(self, xyz_file, step, log_path=None):
        calc = self._make_calculator(
            xyz_file,
            basename="optimise",
            ncpu=step.get("ncpu", 4),
            ram=step.get("ram", 4000),
            charge=step.get("charge", 0),
            mult=step.get("multiplicity", 1),
        )
        self._apply_keywords(calc, step)
        calc.write_input()
        calc.run()  # runs ORCA but returns nothing useful

        calc.create_jsons()

        # retrieve parsed results
        output = calc.get_output()
        output.parse()  # important!

        res = self._wrap_result(calc, output)
        # logger.warning(res) # debug for now

        return res

    def single_point(self, xyz_file, step, log_path=None):
        calc = self._make_calculator(
            xyz_file,
            basename="single_point",
            ncpu=step.get("ncpu", 4),
            ram=step.get("ram", 4000),
            charge=step.get("charge", 0),
            mult=step.get("multiplicity", 1),
        )
        self._apply_keywords(calc, step)
        calc.write_input()
        calc.run()
        calc.create_jsons()

        # retrieve parsed results
        output = calc.get_output()
        output.parse()  # important!

        res = self._wrap_result(calc, output)
        # logger.warning(res) # debug for now

        return res