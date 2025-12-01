import os
import threading
import time
import re
import sys
from pathlib import Path
import json
import numpy as np
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# Core OPI imports
from opi.core import Calculator
from opi.input.structures.structure import Structure
from opi.input.simple_keywords import Task, Dft
from opi.input.blocks.block_scf import BlockScf, DIIS, Shift
from opi.input.blocks.block_cpcm import BlockCpcm


class OrcaBackend:
    def __init__(self, global_config=None):
        self.global_config = global_config or {}
        self.working_dir = Path(
            self.global_config.get("output", {}).get("directory", ".")
        )

    def _wrap_result(self, calc, output):
        prop_path = Path(calc.working_dir) / f"{calc.basename}.property.json"
        if not prop_path.exists():
            raise FileNotFoundError(f"Property JSON not found: {prop_path}")

        with open(prop_path, "r") as f:
            data = json.load(f)

        geoms = data.get("Geometries", [])
        last = geoms[-1]

        energy = last.get("DFT_Energy", {}).get("finalEn") \
            or last["Energy"][0]["totalEnergy"][0][0]
        gibbs_energy = last.get("ThermalData", {}).get("GibbsFreeEnergy", energy)

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

    def _make_calculator(self, xyz_file, basename, ncpu, ram, charge=0, mult=1, working_dir=None):
        calc = Calculator(
            basename=basename,
            working_dir=str(working_dir or self.working_dir)
        )

        structure = Structure.from_xyz(xyz_file)
        structure.charge = int(charge)
        structure.mult = int(mult)
        calc.structure = structure

        calc.input.ncores = int(ncpu)
        calc.input.memory = int(ram)
        return calc

    def _print_progress(self, pct):
        bar = int(pct / 4)
        sys.stdout.write(f"\r[{'#' * bar}{'.' * (25 - bar)}] {pct:3d}%")
        sys.stdout.flush()

    def _monitor_logfile(self, logfile, stop_event, max_cycles=25):
        last_cycle = 0
        while not stop_event.is_set():
            if logfile.exists():
                text = logfile.read_text(errors="ignore")
                matches = re.findall(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)", text)
                if matches:
                    current = int(matches[-1])
                    if current != last_cycle:
                        last_cycle = current
                        pct = min(100, int((current / max_cycles) * 100))
                        self._print_progress(pct)
            time.sleep(0.2)
        self._print_progress(100)
        print("")

    def _norm(self, kw):
        return None if kw is None else str(kw).strip()

    def _apply_keywords(self, calc, step):
        method = step.get("method", "").upper()
        if hasattr(Dft, method):
            func_keyword = getattr(Dft, method)
        else:
            raise ValueError(f"Unknown DFT functional '{method}'")

        basis = self._norm(step.get("basis"))
        disp_kw = self._norm(step.get("dispersion"))
        grid_kw = self._norm(step.get("grid"))
        ri_kw = self._norm(step.get("ri"))

        task_name = step.get("task", "optimise")
        task = Task.OPT if task_name == "optimise" else Task.SP

        keyword_pieces = [func_keyword.keyword]
        for raw in (basis, disp_kw, grid_kw, ri_kw):
            if raw:
                keyword_pieces.append(raw)
        keyword_pieces.append(task.keyword)

        if step.get("freq", False) and task_name == "optimise":
            keyword_pieces.append("FREQ")

        calc.input.add_simple_keywords(" ".join(keyword_pieces))

        if "ram" in step:
            calc.input.maxcore = int(step["ram"])
        if "ncpu" in step:
            calc.input.nprocs = int(step["ncpu"])

        scf_kwargs = {}
        ref = step.get("reference")
        if ref:
            scf_kwargs["reference"] = self._norm(ref)

        for k, v in step.get("scf", {}).items():
            key = k.lower()
            if key == "diis":
                scf_kwargs["diis"] = DIIS(enabled=v) if isinstance(v, bool) else DIIS(**v)
            elif key == "shift":
                scf_kwargs["shift"] = Shift(scf=v) if isinstance(v, (float, int)) else Shift(**v)
            elif key == "scfmode":
                scf_kwargs["scfmode"] = v
            elif key == "maxiter":
                scf_kwargs["maxiter"] = int(v)
            else:
                scf_kwargs[key] = v
        if scf_kwargs:
            calc.input.add_blocks(BlockScf(**scf_kwargs))

        solv = step.get("solvent")
        if isinstance(solv, dict):
            model = solv.get("model", "").lower()
            eps = float(solv.get("epsilon", 80.0))
            solv_name = solv.get("name", "").upper() if solv.get("name") else None
            use_draco = bool(solv.get("draco", False))

            if model == "cpcm":
                kwargs = {"epsilon": eps}
                if solv_name:
                    kwargs["solvent"] = solv_name
                if use_draco:
                    kwargs["draco"] = True
                calc.input.add_blocks(BlockCpcm(**kwargs))
            elif model == "smd":
                calc.input.add_simple_keywords("SMD")

        logger.debug(f"BACKEND RECEIVED SOLVENT BLOCK: {step.get('solvent')}")

    def optimise(self, xyz_file, step, log_path=None):
        job_name = Path(xyz_file).stem
        job_dir = Path(self.working_dir) / job_name
        job_dir.mkdir(parents=True, exist_ok=True)

        calc = self._make_calculator(
            xyz_file,
            basename=job_name,
            ncpu=step.get("ncpu", 4),
            ram=step.get("ram", 4000),
            charge=step.get("charge", 0),
            mult=step.get("multiplicity", 1),
            working_dir=job_dir
        )
        
        self._apply_keywords(calc, step)
        
        calc.write_input()

        logfile = Path(calc.working_dir) / f"{calc.basename}.out"
        stop_event = threading.Event()
        t = threading.Thread(target=self._monitor_logfile, args=(logfile, stop_event))
        t.start()

        calc.run()  # run ORCA via OPI
        stop_event.set()
        t.join()

        calc.create_jsons()

        output = calc.get_output()
        output.parse()
        return self._wrap_result(calc, output)
    
    #TODO add progress tracking to SPE.

    def single_point(self, xyz_file, step, log_path=None):
        job_name = Path(xyz_file).stem
        job_dir = Path(self.working_dir) / job_name
        job_dir.mkdir(parents=True, exist_ok=True)

        calc = self._make_calculator(
            xyz_file,
            basename=job_name,
            ncpu=step.get("ncpu", 4),
            ram=step.get("ram", 4000),
            charge=step.get("charge", 0),
            mult=step.get("multiplicity", 1),
            working_dir=job_dir
        )
        self._apply_keywords(calc, step)
        calc.write_input()
        calc.run()
        calc.create_jsons()

        output = calc.get_output()
        output.parse()
        return self._wrap_result(calc, output)
