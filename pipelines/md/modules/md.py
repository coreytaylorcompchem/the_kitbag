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

        logger.info("Loaded prepared system from disk")

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
        from openmm.unit import kilojoule, mole, nanometer, kelvin
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
        # Gentle perturbation (NOT large jitter)
        # -------------------------
        logger.info("Applying gentle perturbation to lipid atoms")

        def to_nm(val):
            return val.value_in_unit(nanometer) if hasattr(val, "value_in_unit") else float(val)

        positions_array = np.array([[to_nm(p.x), to_nm(p.y), to_nm(p.z)] for p in self.positions])

        noise = np.random.normal(0, 0.02, size=(len(lipid_atoms), 3))  # 0.2 Å
        positions_array[lipid_atoms] += noise

        from openmm import Vec3
        self.positions[:] = [Vec3(*coord) * nanometer for coord in positions_array]

        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(self.positions)
        self.simulation.context.setVelocitiesToTemperature(10 * kelvin)

        # -------------------------
        # Helper: add restraints
        # -------------------------
        def add_restraints(atom_indices, k, stage_id="stage"):
            # unique parameter name per stage
            param_name = f"k_{stage_id}"
            force = CustomExternalForce(f"{param_name}*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            force.addGlobalParameter(param_name, k)
            force.addPerParticleParameter("x0")
            force.addPerParticleParameter("y0")
            force.addPerParticleParameter("z0")

            state = self.simulation.context.getState(getPositions=True)
            pos = state.getPositions()

            for idx in atom_indices:
                p = pos[idx]
                force.addParticle(idx, [p.x, p.y, p.z])

            self.system.addForce(force)
            self.simulation.context.reinitialize(preserveState=True)
            self.simulation.context.setPositions(pos)

            return force

        # -------------------------
        # Staged minimisation to fix bad lipid geometries and other stuff.
        # -------------------------
        stages = [
            {"desc": "Relaxing lipids (protein restrained)", "k": 1000, "atoms": protein_atoms + ligand_atoms, "iter": 500},
            {"desc": "Relaxing headgroups and solvent", "k": 100, "atoms": protein_atoms, "iter": 5000},
            {"desc": "Global relaxation (weak protein restraints)", "k": 10, "atoms": protein_atoms, "iter": 1000},
            {"desc": "Final unrestrained minimisation", "k": None, "atoms": [], "iter": 1000},
        ]

        for i, stage in enumerate(stages, 1):
            logger.info(f"Stage {i}: {stage['desc']}")

            restraint_force = None
            if stage["k"] is not None:
                k_val = stage["k"] * kilojoule / (mole * nanometer**2)
                restraint_force = add_restraints(stage["atoms"], k_val, stage_id=f"{i}")
                logger.info(f"Applied restraints to {len(stage['atoms'])} atoms (k={stage['k']})")

            chunk = 5  
            total_iter = stage["iter"]

            with tqdm(total=total_iter, desc=f"Stage {i} minimisation", unit="iter") as pbar:
                steps_done = 0
                while steps_done < total_iter:
                    n = min(chunk, total_iter - steps_done)

                    self.simulation.minimizeEnergy(
                        tolerance=1 * kilojoule / (mole * nanometer),
                        maxIterations=n
                    )

                    steps_done += n
                    pbar.update(n)

                    # periodic diagnostics
                    if steps_done % 200 == 0 or steps_done == total_iter:
                        state = self.simulation.context.getState(getEnergy=True, getForces=True)
                        energy_val = state.getPotentialEnergy().value_in_unit(kilojoule/mole)
                        forces = state.getForces(asNumpy=True) / (kilojoule / mole / nanometer)

                        max_force = np.max(np.linalg.norm(forces, axis=1))

                        logger.info(
                            f"[Stage {i} | {steps_done}/{total_iter}] "
                            f"Energy: {energy_val:.1f} kJ/mol | Max force: {max_force:.1f}"
                        )

            state = self.simulation.context.getState(getEnergy=True, getForces=True)
            energy_val = state.getPotentialEnergy().value_in_unit(kilojoule/mole)
            forces = state.getForces(asNumpy=True) / (kilojoule / mole / nanometer)

            force_norms = np.linalg.norm(forces, axis=1)
            max_idx = np.argmax(force_norms)
            max_force = force_norms[max_idx]

            atom = atoms[max_idx]
            logger.warning(f"Worst atom: {atom.name} in residue {atom.residue.name} (index {max_idx})")
            logger.warning(f"Max force: {max_force:.2f} kJ/mol/nm")
            logger.info(f"Energy: {energy_val:.2f} kJ/mol")

            if max_force > 50000:
                logger.warning("Extremely high forces persist - check your geometry")

            # remove restraint force
            if restraint_force is not None:
                self.system.removeForce(self.system.getNumForces() - 1)
                self.simulation.context.reinitialize(preserveState=True)

        # -------------------------
        # Final validation
        # -------------------------
        state = self.simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
        forces = state.getForces(asNumpy=True) / (kilojoule / mole / nanometer)
        max_force = np.max(np.linalg.norm(forces, axis=1))

        logger.info(f"Final max force: {max_force:.2f} kJ/mol/nm")

        if max_force > 5000:
            raise RuntimeError("Minimisation failed: forces still too high")

        from openmm.app import PDBFile
        PDBFile.writeFile(self.topology, state.getPositions(), open("minimized.pdb", "w"))
        logger.info("Minimized PDB written: minimized.pdb")

    @register_task("heat_and_equilibrate", category='Molecular dynamics',
               description="Heating and equilibration.")
    def heat_and_equilibrate(self):
        cfg = self.config["heat_and_equilibrate"]

        from openmm import CustomExternalForce, MonteCarloBarostat, NonbondedForce, HarmonicBondForce
        from openmm.unit import kelvin, atmosphere, nanometer, dalton, kilojoule, mole, picoseconds, elementary_charge
        from tqdm import tqdm

        # -------------------------------
        # 0. Ligand parameter and geometry check
        # -------------------------------

        ligand_resname = 'UNK'
        logger.info(f"Checking ligand atoms and bonds (resname={ligand_resname}):")

        atoms_list = list(self.topology.atoms())
        ligand_atoms = [atom.index for atom in atoms_list if atom.residue.name == ligand_resname]

        # --- Nonbonded parameters ---
        for force_index in range(self.system.getNumForces()):
            force = self.system.getForce(force_index)
            if isinstance(force, NonbondedForce):
                for i in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(i)
                    atom = atoms_list[i]
                    if atom.residue.name == ligand_resname:
                        charge_val = charge.value_in_unit(elementary_charge) if hasattr(charge, "value_in_unit") else charge
                        logger.debug(f"Ligand atom {i} ({atom.name}): charge={charge}, sigma={sigma}, epsilon={epsilon}")
                        if abs(charge_val) > 5.0:
                            logger.warning(f"Ligand atom {i} ({atom.name}) has high charge: {charge}")
                        # FIXED UNIT MISMATCH: epsilon is energy/mole
                        if sigma <= 0*nanometer or epsilon < 0*kilojoule/mole:
                            logger.warning(f"Ligand atom {i} ({atom.name}) has bad sigma/epsilon: {sigma}, {epsilon}")

        # --- Bond parameters ---
        for force_index in range(self.system.getNumForces()):
            force = self.system.getForce(force_index)
            if isinstance(force, HarmonicBondForce):
                for i in range(force.getNumBonds()):
                    p1, p2, length, k = force.getBondParameters(i)
                    atom1, atom2 = atoms_list[p1], atoms_list[p2]
                    if atom1.residue.name == ligand_resname or atom2.residue.name == ligand_resname:
                        logger.debug(f"Ligand bond {p1}-{p2} ({atom1.name}-{atom2.name}): length={length}, k={k}")
                        # FIXED UNIT MISMATCH GOD: k has units kJ/mol/nm^2
                        if length < 0.05*nanometer or k > 400000*kilojoule/mole/nanometer**2:
                            logger.warning(f"Ligand bond {p1}-{p2} ({atom1.name}-{atom2.name}) suspicious: length={length}, k={k}")

        # # --- Positions & neighbor clashes ---
        # logger.info("Checking positions of ligand atoms and neighbors:")
        # state = self.simulation.context.getState(getPositions=True)
        # positions = state.getPositions(asNumpy=True)

        # neighbor_atoms = [atom.index for atom in atoms_list if atom.index not in ligand_atoms]
        # clash_cutoff = 0.15*nanometer

        # for l_idx in ligand_atoms:
        #     l_pos = positions[l_idx]
        #     l_pos_nm = [l_pos[0].value_in_unit(nanometer),
        #                 l_pos[1].value_in_unit(nanometer),
        #                 l_pos[2].value_in_unit(nanometer)]
        #     if any(abs(c) > 10.0 for c in l_pos_nm):
        #         logger.warning(f"Ligand atom {l_idx} ({atoms_list[l_idx].name}) extreme position: {l_pos_nm}")
        #     else:
        #         logger.debug(f"Ligand atom {l_idx} position: {l_pos_nm}")

        #     for n_idx in neighbor_atoms:
        #         n_pos = positions[n_idx]
        #         dx = l_pos[0] - n_pos[0]
        #         dy = l_pos[1] - n_pos[1]
        #         dz = l_pos[2] - n_pos[2]
        #         dist = (dx*dx + dy*dy + dz*dz)**0.5
        #         if dist < clash_cutoff:
        #             logger.warning(
        #                 f"Close contact: Ligand atom {l_idx} ({atoms_list[l_idx].name}) "
        #                 f"and atom {n_idx} ({atoms_list[n_idx].name}, {atoms_list[n_idx].residue.name}) "
        #                 f"distance = {dist.value_in_unit(nanometer):.3f} nm"
        #             )

        # -------------------------------
        # 1. HEATING PHASE (staged)
        # -------------------------------
        heating_cfg = cfg.get("heating", {})
        num_rounds = heating_cfg.get("num_steps", 8)  # increased for gentler ramp
        steps_per_round = min(heating_cfg.get("steps_per_round", 5000), 500)
        target_temp = heating_cfg.get("target_temp", 300)
        # stronger initial restraints and more gradual decay
        restraints = heating_cfg.get("restraint_strengths", [2.0, 1.5, 1.0, 0.5, 0.25, 0.1, 0.05, 0.0])

        logger.debug(heating_cfg)

        logger.info(f"Heating with {num_rounds} rounds up to {target_temp} K")

        backbone_atoms = [atom.index for atom in atoms_list
                        if atom.residue.name in protein_resnames and atom.name in {"N","CA","C","O"}]
        ligand_atoms = [atom.index for atom in atoms_list if atom.residue.name == ligand_resname]

        # restraint force
        restraint_force = CustomExternalForce("k_heating*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint_force.addGlobalParameter("k_heating", restraints[0])
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")

        state = self.simulation.context.getState(getPositions=True)
        positions = state.getPositions()

        for atom_idx in backbone_atoms + ligand_atoms:
            pos = positions[atom_idx]
            restraint_force.addParticle(atom_idx, [pos.x, pos.y, pos.z])

        self.system.addForce(restraint_force)
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(positions)
        self.simulation.context.setVelocitiesToTemperature(10*kelvin)

        # Diagnostics BEFORE heating
        state = self.simulation.context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy()
        forces = state.getForces(asNumpy=True)
        forces_kj_per_mol_nm = forces / (kilojoule/mole/nanometer)
        max_force_value = max((f[0]**2 + f[1]**2 + f[2]**2)**0.5 for f in forces_kj_per_mol_nm)
        logger.info(f"Initial potential energy: {energy.value_in_unit(kilojoule/mole):.2f} kJ/mol")
        logger.info(f"Initial max force: {max_force_value:.2f} kJ/mol/nm")

        # Heating loop
        original_dt = self.integrator.getStepSize()
        small_dt = 0.00025 * picoseconds  # 0.25 fs initially
        step_chunk = 50

        for i in range(num_rounds):
            temp = 10 + (target_temp - 10) * (i + 1) / num_rounds
            k_value = restraints[i] if i < len(restraints) else 0.0
            self.integrator.setTemperature(temp * kelvin)
            self.simulation.context.setParameter("k_heating", k_value)

            # Use small_dt for first 3 rounds
            if i < 3:
                self.integrator.setStepSize(small_dt)
            else:
                self.integrator.setStepSize(original_dt)

            steps_done = 0
            with tqdm(total=steps_per_round, desc=f"Heating round {i+1}", unit="steps") as pbar:
                while steps_done < steps_per_round:
                    n = min(step_chunk, steps_per_round - steps_done)
                    try:
                        self.simulation.step(n)
                    except Exception as e:
                        state = self.simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
                        forces = state.getForces(asNumpy=True)
                        max_force = max((f[0]**2 + f[1]**2 + f[2]**2)**0.5 for f in forces / (kilojoule/mole/nanometer))
                        logger.error(f"Crash during heating round {i+1} at step {steps_done}")
                        logger.error(f"Max force before crash: {max_force:.2f} kJ/mol/nm")
                        raise e

                    steps_done += n
                    pbar.update(n)

                    # Monitor forces
                    state = self.simulation.context.getState(getForces=True)
                    max_force = max((f[0]**2 + f[1]**2 + f[2]**2)**0.5
                                    for f in state.getForces(asNumpy=True) / (kilojoule/mole/nanometer))
                    if max_force > 3000:
                        logger.warning(f"High force detected: {max_force:.2f} kJ/mol/nm at step {steps_done}")

        self.integrator.setStepSize(original_dt)
        self.system.removeForce(self.system.getNumForces() - 1)
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(positions)
        logger.info("Heating complete; restraints removed")

        # -------------------------------
        # 2. EQUILIBRATION PHASE
        # -------------------------------
        equil_cfg = cfg.get("equilibration", {})
        ensemble = equil_cfg.get("ensemble", "NPT").upper()
        steps = equil_cfg.get("steps", 100000)
        equil_temp = equil_cfg.get("target_temp", target_temp)

        logger.info(f"Starting equilibration: {ensemble} ensemble for {steps} steps at {equil_temp} K")

        has_barostat = any(isinstance(f, MonteCarloBarostat)
                        for f in [self.system.getForce(i) for i in range(self.system.getNumForces())])
        if ensemble == "NPT" and not has_barostat:
            pressure = equil_cfg.get("pressure", 1.0) * atmosphere
            frequency = equil_cfg.get("barostat", {}).get("frequency", 25)
            barostat = MonteCarloBarostat(pressure, equil_temp * kelvin, frequency)
            self.system.addForce(barostat)
            self.simulation.context.reinitialize(preserveState=True)
            self.simulation.context.setPositions(positions)
            logger.info(f"Added MonteCarloBarostat: {pressure} atm, {equil_temp} K")

        self.integrator.setTemperature(equil_temp * kelvin)

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

        # Compute Center of Mass
        com = np.mean(coords[backbone_atoms], axis=0)

        # Compute solvent box lengths
        box_lengths = np.array([
            box_vectors[0][0].value_in_unit(nanometer) if hasattr(box_vectors[0][0], "value_in_unit") else float(box_vectors[0][0]),
            box_vectors[1][1].value_in_unit(nanometer) if hasattr(box_vectors[1][1], "value_in_unit") else float(box_vectors[1][1]),
            box_vectors[2][2].value_in_unit(nanometer) if hasattr(box_vectors[2][2], "value_in_unit") else float(box_vectors[2][2])
        ], dtype=float)

        # Compute shift vs COM
        shift = box_lengths / 2 - com

        # Apply shift to positions
        for i in range(len(positions)):
            x, y, z = coords[i] + shift
            positions[i] = Vec3(x, y, z) * nanometer

        # Update simulation context
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
        split_ns = cfg.get("output_split_ns", ns)          # new trajectory every N ns
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

            # Attach new reporters for this segment
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





