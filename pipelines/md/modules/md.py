import os

from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm.unit import picoseconds
from openmm import Vec3
from openmm.unit import nanometer, kelvin, atmosphere, picoseconds, kilojoule, mole, nanometer as nm

import numpy as np
from tqdm import tqdm

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

protein_resnames = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
    "ILE","LEU","LYS","MET","PHE","PRO","SER","THR",
    "TRP","TYR","VAL"
}

class MDWorkflow:
    def __init__(self, backend, config):
        self.backend = backend
        self.config = config
        self.integrator = None
        self.simulation = None
    from openmm import XmlSerializer

    def load_prepared_system(self, directory):
        topology_path = os.path.join(directory, "topology.pdb")
        system_path = os.path.join(directory, "system.xml")

        if not os.path.exists(topology_path) or not os.path.exists(system_path):
            raise FileNotFoundError(
                "Prepared system not found. Run prepare_system first or provide topology.pdb and system.xml"
            )

        pdb = PDBFile(topology_path)

        with open(system_path) as f:
            system = XmlSerializer.deserialize(f.read())

        self.topology = pdb.topology
        self.positions = pdb.positions
        self.system = system

        logger.info(f"Loaded prepared {system_path} from disk")

    @register_task("prepare_system", category='Molecular dynamics',
                   description="Load inputs, cap chains, parameterise ligand and protein, solvate and save final topology.")
    def prepare_system(self):
        self.backend.prepare_system(self.config["prepare_system"])
        self.system = self.backend.system
        self.topology = self.backend.topology
        self.positions = self.backend.positions

    @register_task(
        "setup_simulation",
        category="Molecular dynamics",
        description="Set up integrator and simulation."
    )
    def setup_simulation(self):

        if not hasattr(self, "system"):
            output_dir = self.config.get("prepare_system", {}).get(
                "output_trajectory",
                "output"
            )
            self.load_prepared_system(output_dir)

        if not hasattr(self, "platform"):
            platform_name = self.config.get("setup_simulation", {}).get("platform", "CPU")
            self.platform = Platform.getPlatformByName(platform_name)

        cfg = self.config["setup_simulation"]

        self.integrator = LangevinIntegrator(
            cfg["temperature"] * kelvin,
            1.0 / picoseconds,
            cfg["timestep"] * femtoseconds
        )

        self.simulation = Simulation(
            self.topology,
            self.system,
            self.integrator,
            self.platform
        )

        self.simulation.context.setPositions(self.positions)

    @register_task("minimize")
    def minimize(self):
        from openmm import CustomExternalForce
        from openmm.unit import kilojoule, mole, nanometer
        import numpy as np
        from tqdm import tqdm

        logger.info("Running staged minimisation")

        atoms = list(self.topology.atoms())

        # -------------------------
        # Identify groups
        # -------------------------
        protein_atoms = [a.index for a in atoms if a.residue.name in protein_resnames]
        ligand_atoms = [a.index for a in atoms if a.residue.name == "UNK"]
        lipid_atoms = [a.index for a in atoms if a.residue.name == "POP"]

        logger.debug(f"Protein atoms: {len(protein_atoms)}")
        logger.debug(f"Ligand atoms: {len(ligand_atoms)}")
        logger.debug(f"Lipid atoms: {len(lipid_atoms)}")

        # -------------------------
        # Single restraint forces (add ONCE)
        # -------------------------
        prot_force = CustomExternalForce("k_prot*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        prot_force.addGlobalParameter("k_prot", 0.0)
        prot_force.addPerParticleParameter("x0")
        prot_force.addPerParticleParameter("y0")
        prot_force.addPerParticleParameter("z0")
        lig_force = CustomExternalForce("k_lig*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        lig_force.addGlobalParameter("k_lig", 0.0)
        lig_force.addPerParticleParameter("x0")
        lig_force.addPerParticleParameter("y0")
        lig_force.addPerParticleParameter("z0")

        # Get current positions
        state = self.simulation.context.getState(getPositions=True)
        pos = state.getPositions()

        for idx in protein_atoms:
            p = pos[idx]
            prot_force.addParticle(idx, [p.x, p.y, p.z])
        for idx in ligand_atoms:
            p = pos[idx]
            lig_force.addParticle(idx, [p.x, p.y, p.z])

        # Add forces ONCE before minimisation loop
        prot_idx = self.system.addForce(prot_force)
        lig_idx = self.system.addForce(lig_force)

        # Reinitialize context once and preserve positions
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(pos)

        # -------------------------
        # Map restrained particles for diagnostics
        # -------------------------
        restrained_atoms = protein_atoms + ligand_atoms
        particle_to_atom = {idx: atoms[idx] for idx in restrained_atoms}

        # -------------------------
        # Staged minimisation
        # -------------------------
        stages = [
            {"desc": "Relaxing lipids (protein + lig restrained)", "k_prot": 1000, "k_lig": 1000, "iter": 200, "tolerance": 100},
            {"desc": "Relaxing headgroups (protein restrained)", "k_prot": 100, "k_lig": 0, "iter": 200, "tolerance": 10},
            {"desc": "Global relaxation (weak restraints)", "k_prot": 10, "k_lig": 0, "iter": 200, "tolerance": 5},
            {"desc": "Final unrestrained minimisation", "k_prot": 0, "k_lig": 0, "iter": 200, "tolerance": 1},
        ]

        for i, stage in enumerate(stages, 1):
            logger.info(f"Stage {i}: {stage['desc']}")

            # Update restraint strengths via global parameters
            self.simulation.context.setParameter("k_prot", stage["k_prot"] * kilojoule/(mole*nanometer**2))
            self.simulation.context.setParameter("k_lig", stage["k_lig"] * kilojoule/(mole*nanometer**2))
            logger.info(f"Applied restraints: protein={stage['k_prot']}, ligand={stage['k_lig']}")

            # Minimisation loop with chunking
            chunk = 5
            total_iter = stage["iter"]
            tolerance = stage.get("tolerance", 10) * kilojoule/(mole*nanometer)
            with tqdm(total=total_iter, desc=f"Stage {i} minimisation", unit="iter") as pbar:
                steps_done = 0
                while steps_done < total_iter:
                    n = min(chunk, total_iter - steps_done)
                    self.simulation.minimizeEnergy(tolerance=tolerance, maxIterations=n)
                    steps_done += n
                    pbar.update(n)

                    # Periodic diagnostics
                    if steps_done % 200 == 0 or steps_done == total_iter:
                        state = self.simulation.context.getState(getEnergy=True, getForces=True)
                        energy_val = state.getPotentialEnergy().value_in_unit(kilojoule/mole)
                        forces = state.getForces(asNumpy=True) / (kilojoule/mole/nanometer)
                        max_force = np.max(np.linalg.norm(forces, axis=1))

                        logger.info(
                            f"[Stage {i} | {steps_done}/{total_iter}] "
                            f"Energy: {energy_val:.1f} kJ/mol | Max force: {max_force:.1f}"
                        )

            # Diagnostics restricted to restrained atoms
            state = self.simulation.context.getState(getEnergy=True, getForces=True)
            forces = state.getForces(asNumpy=True) / (kilojoule/mole/nanometer)
            force_norms = np.linalg.norm(forces, axis=1)
            restrained_forces = {idx: force_norms[idx] for idx in restrained_atoms}
            max_idx = max(restrained_forces, key=restrained_forces.get)
            max_force = restrained_forces[max_idx]
            atom = particle_to_atom[max_idx]

            logger.warning(
                f"Worst restrained atom: {atom.name} in residue {atom.residue.name} "
                f"(resid {atom.residue.id}, chain {atom.residue.chain.id}, system particle index {max_idx})"
            )
            logger.warning(f"Max restrained force: {max_force:.2f} kJ/mol/nm")

            # Optional: warn if unrestrained atoms have very high forces
            unrestrained_indices = np.setdiff1d(np.arange(len(forces)), restrained_atoms)
            if len(unrestrained_indices) > 0:
                unres_forces = force_norms[unrestrained_indices]
                if np.max(unres_forces) > 50000:
                    logger.warning(f"High forces on unrestrained atoms: max {np.max(unres_forces):.1f} kJ/mol/nm")

        # -------------------------
        # Remove restraints safely (after all stages)
        # -------------------------
        self.system.removeForce(lig_idx)
        self.system.removeForce(prot_idx)
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(pos)

        # -------------------------
        # Final check
        # -------------------------
        state = self.simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
        forces = state.getForces(asNumpy=True) / (kilojoule/mole/nanometer)
        max_force = np.max(np.linalg.norm(forces, axis=1))
        logger.info(f"Final max force: {max_force:.2f} kJ/mol/nm")
        if max_force > 10000:
            raise RuntimeError("Minimisation failed: forces still too high")

        from openmm.app import PDBFile
        PDBFile.writeFile(self.topology, state.getPositions(), open("minimized.pdb", "w"))
        logger.info("Minimized PDB written: minimized.pdb")

    @register_task("heat_and_equilibrate", category='Molecular dynamics',
               description="Heating and equilibration.")
    def heat_and_equilibrate(self):
        from openmm import CustomExternalForce, MonteCarloBarostat, NonbondedForce, HarmonicBondForce
        from openmm.unit import kelvin, atmosphere, nanometer, kilojoule, mole, picoseconds
        import numpy as np

        ligand_resname = 'UNK'
        atoms_list = list(self.topology.atoms())
        ligand_atoms = [atom.index for atom in atoms_list if atom.residue.name == ligand_resname]

        # -------------------------------
        # 0. Check ligand parameters
        # -------------------------------
        for force_index in range(self.system.getNumForces()):
            force = self.system.getForce(force_index)
            if isinstance(force, NonbondedForce):
                for i in range(force.getNumParticles()):
                    atom = atoms_list[i]
                    if atom.residue.name == ligand_resname:
                        charge, sigma, epsilon = force.getParticleParameters(i)
                        if sigma <= 0*nanometer or epsilon < 0*kilojoule/mole:
                            logger.warning(f"Ligand atom {i} ({atom.name}) bad sigma/epsilon: {sigma}, {epsilon}")
            elif isinstance(force, HarmonicBondForce):
                for i in range(force.getNumBonds()):
                    p1, p2, length, k = force.getBondParameters(i)
                    atom1, atom2 = atoms_list[p1], atoms_list[p2]
                    if atom1.residue.name == ligand_resname or atom2.residue.name == ligand_resname:
                        if length < 0.05*nanometer or k > 400000*kilojoule/mole/nanometer**2:
                            logger.warning(f"Ligand bond {p1}-{p2} suspicious: length={length}, k={k}")

        # -------------------------------
        # 1. HEATING PHASE (staged)
        # -------------------------------
        heating_cfg = self.config.get("heat_and_equilibrate", {}).get("heating", {})
        num_rounds = heating_cfg.get("num_steps", 8)
        steps_per_round = min(heating_cfg.get("steps_per_round", 5000), 500)
        target_temp = heating_cfg.get("target_temp", 300)
        restraints = heating_cfg.get("restraint_strengths", [2.0, 1.5, 1.0, 0.5, 0.25, 0.1, 0.05, 0.0])

        logger.info(f"Heating with {num_rounds} rounds up to {target_temp} K")

        # backbone + ligand atoms to restrain
        backbone_atoms = [atom.index for atom in atoms_list
                        if atom.residue.name in protein_resnames and atom.name in {"N","CA","C","O"}]
        restrain_atoms = backbone_atoms + ligand_atoms

        # -------------------------------
        # Create restraint force (single reinit)
        # -------------------------------
        restraint_force = CustomExternalForce("k_heating*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint_force.addGlobalParameter("k_heating", restraints[0] * kilojoule/(mole*nanometer**2))
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")

        state = self.simulation.context.getState(getPositions=True)
        positions = state.getPositions()
        for idx in restrain_atoms:
            pos = positions[idx]
            restraint_force.addParticle(idx, [pos.x, pos.y, pos.z])

        # Debug: check initial displacement before adding to system
        max_disp = 0.0
        for i in range(restraint_force.getNumParticles()):
            particle_index, params = restraint_force.getParticleParameters(i)
            x0, y0, z0 = params
            pos = positions[particle_index]
            disp = np.linalg.norm([pos.x - x0, pos.y - y0, pos.z - z0])
            max_disp = max(max_disp, disp)
        logger.info(f"[DEBUG] Max initial displacement before adding restraint: {max_disp:.6f} nm")

        # Add to system and reinitialize context once
        self.system.addForce(restraint_force)
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(positions)
        self.simulation.context.setVelocitiesToTemperature(10*kelvin)

        # Debug: check displacements after reinit
        state = self.simulation.context.getState(getPositions=True)
        positions_post = state.getPositions()
        max_disp_post = 0.0
        for i in range(restraint_force.getNumParticles()):
            particle_index, params = restraint_force.getParticleParameters(i)
            x0, y0, z0 = params
            pos = positions_post[particle_index]
            disp = np.linalg.norm([pos.x - x0, pos.y - y0, pos.z - z0])
            max_disp_post = max(max_disp_post, disp)
        logger.info(f"[DEBUG] Max displacement after reinit: {max_disp_post:.6f} nm")

        # Debug: check initial forces
        state = self.simulation.context.getState(getForces=True)
        forces = state.getForces(asNumpy=True) / (kilojoule / mole / nanometer)
        # max force on restrained atoms only
        restrained_indices = [restraint_force.getParticleParameters(i)[0] for i in range(restraint_force.getNumParticles())]
        max_force = np.max(np.linalg.norm(forces[restrained_indices], axis=1))
        logger.info(f"[DEBUG] Max force on restrained atoms before stepping: {max_force:.2f} kJ/mol/nm")

        # -------------------------------
        # Heating loop with step-level debug
        # -------------------------------
        original_dt = self.integrator.getStepSize()
        small_dt = 0.00025 * picoseconds  # 0.25 fs initially
        step_chunk = 50

        for i in range(num_rounds):
            temp = 10 + (target_temp - 10) * (i + 1) / num_rounds
            k_value = restraints[i] if i < len(restraints) else 0.0
            self.integrator.setTemperature(temp * kelvin)
            self.simulation.context.setParameter("k_heating", k_value * kilojoule/(mole*nanometer**2))
            self.integrator.setStepSize(small_dt if i < 3 else original_dt)

            steps_done = 0
            with tqdm(total=steps_per_round, desc=f"Heating round {i+1}", unit="steps") as pbar:
                while steps_done < steps_per_round:
                    n = min(step_chunk, steps_per_round - steps_done)
                    try:
                        self.simulation.step(n)
                    except Exception as e:
                        # Force check immediately on crash
                        state = self.simulation.context.getState(getForces=True)
                        forces = state.getForces(asNumpy=True) / (kilojoule/mole/nanometer)
                        max_force = np.max(np.linalg.norm(forces, axis=1))
                        top5_idx = np.argsort(np.linalg.norm(forces, axis=1))[-5:]
                        logger.error(f"Crash during heating round {i+1} at step {steps_done}")
                        logger.error(f"Max force before crash: {max_force:.2f} kJ/mol/nm")
                        for idx in top5_idx:
                            atom = atoms_list[idx]
                            fval = np.linalg.norm(forces[idx])
                            logger.error(f"[DEBUG] Atom {atom.name} (res {atom.residue.name}, idx {idx}) force: {fval:.2f} kJ/mol/nm")
                        raise e

                    steps_done += n
                    pbar.update(n)

                    # Step-level force debug
                    state = self.simulation.context.getState(getForces=True)
                    forces = state.getForces(asNumpy=True) / (kilojoule/mole/nanometer)
                    max_force = np.max(np.linalg.norm(forces, axis=1))
                    if max_force > 3000:
                        top5_idx = np.argsort(np.linalg.norm(forces, axis=1))[-5:]
                        logger.warning(f"High force detected: {max_force:.2f} kJ/mol/nm at step {steps_done}")
                        for idx in top5_idx:
                            atom = atoms_list[idx]
                            fval = np.linalg.norm(forces[idx])
                            logger.warning(f"[DEBUG] Atom {atom.name} (res {atom.residue.name}, idx {idx}) force: {fval:.2f} kJ/mol/nm")

        # -------------------------------
        # Remove heating restraint safely
        # -------------------------------
        self.system.removeForce(self.system.getForceIndex(restraint_force))
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(state.getPositions())
        self.simulation.context.setVelocitiesToTemperature(target_temp * kelvin)
        logger.info("Heating complete; restraints removed")

        # -------------------------------
        # 2. EQUILIBRATION PHASE
        # -------------------------------
        equil_cfg = self.config.get("heat_and_equilibrate", {}).get("equilibration", {})
        ensemble = equil_cfg.get("ensemble", "NPT").upper()
        steps = equil_cfg.get("steps", 100000)
        equil_temp = equil_cfg.get("target_temp", target_temp)

        logger.info(f"Starting equilibration: {ensemble} ensemble for {steps} steps at {equil_temp} K")

        # Add barostat if needed
        if ensemble == "NPT" and not any(isinstance(self.system.getForce(i), MonteCarloBarostat) 
                                        for i in range(self.system.getNumForces())):
            pressure = equil_cfg.get("pressure", 1.0) * atmosphere
            frequency = equil_cfg.get("barostat", {}).get("frequency", 25)
            barostat = MonteCarloBarostat(pressure, equil_temp * kelvin, frequency)
            self.system.addForce(barostat)
            self.simulation.context.reinitialize(preserveState=True)
            self.simulation.context.setPositions(state.getPositions())
            self.simulation.context.setVelocitiesToTemperature(equil_temp * kelvin)
            logger.info(f"Added MonteCarloBarostat: {pressure} atm, {equil_temp} K")

        self.integrator.setTemperature(equil_temp * kelvin)

        # Equilibration loop
        step_chunk = 1000
        steps_done = 0
        with tqdm(total=steps, desc="Equilibration", unit="steps") as pbar:
            while steps_done < steps:
                n = min(step_chunk, steps - steps_done)
                self.simulation.step(n)
                steps_done += n
                pbar.update(n)

        # Diagnostics
        state = self.simulation.context.getState(getEnergy=True, getPositions=True)
        box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
        volume_nm3 = (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]).value_in_unit(nanometer**3)
        mass = sum([self.system.getParticleMass(i) for i in range(self.system.getNumParticles())])
        mass_g = mass.value_in_unit(dalton) * 1.66054e-24
        volume_mL = volume_nm3 * 1e-21
        density = mass_g / volume_mL
        logger.info(f"Equilibration complete. Box volume: {volume_nm3:.3f} nm³, density: {density:.3f} g/mL")
        if abs(density - 1.0) > 0.1:
            logger.warning("Density deviates from expected ~1.0 g/mL. Check for bubbles or equilibration issues.")


    @register_task("production", category='Molecular dynamics', description="Run production simulation.")
    def production(self):

        cfg = self.config["production"]
        ns = cfg.get("length_ns", 1)
        ensemble = cfg.get("ensemble", "NPT").upper()
        target_temp = cfg.get("target_temp", self.integrator.getTemperature().value_in_unit(kelvin))
        pressure = cfg.get("pressure", 1.0) 

        # Get current positions, velocities, and box vectors
        state = self.simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
        positions = state.getPositions()
        velocities = state.getVelocities()
        box_vectors = state.getPeriodicBoxVectors()

        #### To center in the box 
        # Convert positions to plain floats in nm
        coords = np.array([
            [
                p.x.value_in_unit(nanometer) if hasattr(p.x, "value_in_unit") else float(p.x),
                p.y.value_in_unit(nanometer) if hasattr(p.y, "value_in_unit") else float(p.y),
                p.z.value_in_unit(nanometer) if hasattr(p.z, "value_in_unit") else float(p.z)
            ]
            for p in positions
        ], dtype=float)

        # Get backbone atom indices
        backbone_atoms = [atom.index for atom in self.topology.atoms() if atom.name in {"N", "CA", "C", "O"}]

        # Compute com
        com = np.mean(coords[backbone_atoms], axis=0)

        # Compute solvent box lengths
        box_lengths = np.array([
            box_vectors[0][0].value_in_unit(nanometer) if hasattr(box_vectors[0][0], "value_in_unit") else float(box_vectors[0][0]),
            box_vectors[1][1].value_in_unit(nanometer) if hasattr(box_vectors[1][1], "value_in_unit") else float(box_vectors[1][1]),
            box_vectors[2][2].value_in_unit(nanometer) if hasattr(box_vectors[2][2], "value_in_unit") else float(box_vectors[2][2])
        ], dtype=float)

        # shift vs com
        shift = box_lengths / 2 - com

        # Apply shift to positions
        for i in range(len(positions)):
            x, y, z = coords[i] + shift
            positions[i] = Vec3(x, y, z) * nanometer

        # Update context
        self.simulation.context.setPositions(positions)
        self.simulation.context.setVelocities(velocities)
        self.simulation.context.setPeriodicBoxVectors(*box_vectors)
        self.integrator.setTemperature(target_temp * kelvin)
        self.simulation.currentStep = 0

        # Attach barostat
        has_barostat = any(isinstance(f, MonteCarloBarostat) for f in
                        [self.system.getForce(i) for i in range(self.system.getNumForces())])
        if ensemble == "NPT" and not has_barostat:
            barostat = MonteCarloBarostat(pressure * atmosphere, target_temp * kelvin, 25)
            self.system.addForce(barostat)
            self.simulation.context.reinitialize(preserveState=True)

            # Restore positions and velocities
            self.simulation.context.setPositions(positions)
            self.simulation.context.setVelocities(velocities)
            self.simulation.context.setPeriodicBoxVectors(*box_vectors)
            logger.info(f"Added MonteCarloBarostat: {pressure} atm, {target_temp} K")

        # Compute total steps
        timestep = self.integrator.getStepSize()
        timestep_ns = timestep.value_in_unit(picoseconds) / 1000.0
        steps_per_ns = int(1.0 / timestep_ns)
        total_steps = int(ns * steps_per_ns)

        # Output control
        split_ns = cfg.get("output_split_ns", ns) 
        split_steps = int(split_ns * steps_per_ns)
        n_segments = int(np.ceil(ns / split_ns))

        output_freq = cfg.get("output_frequency", 1000)
        output_trajectory_base = cfg.get("output_trajectory", "output.dcd")
        output_log_base = cfg.get("output_logfile", "output.log")
        output_dir = os.path.dirname(output_trajectory_base)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        logger.info(f"Running Production: {ns} ns total ({total_steps} steps), "
                    f"split every {split_ns} ns into {n_segments} segments. "
                    f"Writing frames every {output_freq} steps (~{output_freq * timestep.value_in_unit(femtoseconds)/1000:.1f} ps).")

        # Loop over N splits specified in the yaml
        for seg in range(1, n_segments + 1):
            seg_start_step = self.simulation.currentStep
            seg_end_step = min(seg_start_step + split_steps, total_steps)
            seg_steps = seg_end_step - seg_start_step

            traj_file = os.path.join(output_dir, f"trajectory_{seg:03d}.dcd")
            log_file = os.path.join(output_dir, f"log_{seg:03d}.txt")
            chk_file = os.path.join(output_dir, f"restart_{seg:03d}.chk")

            # Attach reporters
            self.simulation.reporters.clear()
            self.simulation.reporters.append(DCDReporter(traj_file, output_freq))
            self.simulation.reporters.append(StateDataReporter(
                log_file, output_freq,
                step=True, temperature=True, potentialEnergy=True, totalEnergy=True,
                density=True, progress=True, remainingTime=True, speed=True,
                totalSteps=seg_steps, separator='\t'
            ))

            logger.info(f"→ Segment {seg}/{n_segments}: running {seg_steps} steps "
                        f"({seg_steps/steps_per_ns:.3f} ns) → {traj_file}")

            with tqdm(total=seg_steps, desc=f"Production segment {seg}", unit="steps") as pbar:
                current_step = 0
                while current_step < seg_steps:
                    n = min(output_freq, seg_steps - current_step)
                    self.simulation.step(n)
                    current_step += n
                    pbar.update(n)

            # Save restart at end of segment
            self.simulation.saveCheckpoint(chk_file)
            logger.info(f"Checkpoint saved: {chk_file}")

        logger.info("Production complete.")


                # # Wrap atoms every 10k steps
                # This doesn't work and may not be necessary anyway so commenting out for now.
                # if self.simulation.currentStep % 10000 == 0:
                #     state = self.simulation.context.getState(getPositions=True)
                #     positions = state.getPositions()
                #     for i in range(len(positions)):
                #         x = (positions[i].x / nanometer) % box_lengths[0]
                #         y = (positions[i].y / nanometer) % box_lengths[1]
                #         z = (positions[i].z / nanometer) % box_lengths[2]
                #         positions[i] = Vec3(float(x), float(y), float(z)) * nanometer
                #     self.simulation.context.setPositions(positions)





