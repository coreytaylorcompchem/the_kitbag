import os
import json
import sys
import platform

from tqdm import tqdm

from openmm.app import *
from openmm.app import PDBFile
from openmm import *
from openmm import MonteCarloBarostat
from openmm import XmlSerializer
from openmm import CustomExternalForce, MonteCarloBarostat, NonbondedForce, HarmonicBondForce
from openmm.unit import *
from openmm.unit import picoseconds
from openmm import Vec3
from openmm.unit import nanometer, kelvin, atmosphere, picoseconds, kilojoule, mole

import numpy as np
import pandas as pd

import MDAnalysis as mda
from MDAnalysis.transformations import unwrap, center_in_box, wrap

from modules.utils.plotting import save_plot

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

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
        
        # Set GPU precision: "single", "mixed", or "double"
        if self.platform.getName() in ["CUDA", "OpenCL"]:
            self.platform.setPropertyDefaultValue("Precision", "mixed")

        cfg = self.config["setup_simulation"]

        # -----------------------------
        # Hydrogen Mass Repartitioning helper
        # -----------------------------
        def apply_hmr(system, topology, target_h_mass=3.0*dalton):
            for bond in topology.bonds():
                atom1, atom2 = bond[0], bond[1]
                mass1 = system.getParticleMass(atom1.index)
                mass2 = system.getParticleMass(atom2.index)
                
                # Check if one atom is hydrogen and the other is heavy
                if mass1.value_in_unit(dalton) < 1.2 and mass2.value_in_unit(dalton) > 1.5:
                    delta = target_h_mass.value_in_unit(dalton) - mass1.value_in_unit(dalton)
                    system.setParticleMass(atom1.index, target_h_mass)
                    system.setParticleMass(atom2.index, mass2 - delta * dalton)
                elif mass2.value_in_unit(dalton) < 1.2 and mass1.value_in_unit(dalton) > 1.5:
                    delta = target_h_mass.value_in_unit(dalton) - mass2.value_in_unit(dalton)
                    system.setParticleMass(atom2.index, target_h_mass)
                    system.setParticleMass(atom1.index, mass1 - delta * dalton)
            return system

        if cfg.get("hmr", True):
            self.system = apply_hmr(self.system, self.topology)
            logger.info("Hydrogen Mass Repartitioning (HMR) applied")

        # Create integrator with larger timestep
        timestep_fs = cfg["timestep"] * femtoseconds
        if cfg.get("hmr", True):
            timestep_fs *= 2  # safe doubling of timestep
        self.integrator = LangevinIntegrator(cfg["temperature"] * kelvin, 1.0/picoseconds, timestep_fs)

        self.simulation = Simulation(
            self.topology,
            self.system,
            self.integrator,
            self.platform
        )
        self.simulation.context.setPositions(self.positions)

    @register_task("minimize")
    def minimize(self):

        logger.info("Running staged minimisation")

        atoms = list(self.topology.atoms())

        # Minimisation groups

        protein_atoms = [a.index for a in atoms if a.residue.name in protein_resnames]
        ligand_atoms = [a.index for a in atoms if a.residue.name == "UNK"]
        lipid_atoms = [a.index for a in atoms if a.residue.name == "POP"]

        logger.debug(f"Protein atoms: {len(protein_atoms)}")
        logger.debug(f"Ligand atoms: {len(ligand_atoms)}")
        logger.debug(f"Lipid atoms: {len(lipid_atoms)}")


        # Add restraints

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

        # Add forces once before minimisation loop
        prot_idx = self.system.addForce(prot_force)
        lig_idx = self.system.addForce(lig_force)

        # Reinitialise context once and preserve positions
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.context.setPositions(pos)


        # Map restrained particles for diagnostics

        restrained_atoms = protein_atoms + ligand_atoms
        particle_to_atom = {idx: atoms[idx] for idx in restrained_atoms}

        # -------------------------
        # Staged minimisation
        # -------------------------
        stages = [
            {"desc": "Relaxing lipids (protein + lig restrained)", "k_prot": 1000, "k_lig": 1000, "iter": 2000, "tolerance": 100},
            {"desc": "Relaxing headgroups (protein restrained)", "k_prot": 100, "k_lig": 0, "iter": 5000, "tolerance": 10},
            {"desc": "Global relaxation (weak restraints)", "k_prot": 10, "k_lig": 0, "iter": 2000, "tolerance": 5},
            {"desc": "Final unrestrained minimisation", "k_prot": 0, "k_lig": 0, "iter": 1000, "tolerance": 1},
        ]

        for i, stage in enumerate(stages, 1):
            logger.info(f"Stage {i}: {stage['desc']}")

            # Update restraint strengths via global parameters
            self.simulation.context.setParameter("k_prot", stage["k_prot"] * kilojoule/(mole*nanometer**2))
            self.simulation.context.setParameter("k_lig", stage["k_lig"] * kilojoule/(mole*nanometer**2))
            logger.debug(f"Applied restraints: protein={stage['k_prot']}, ligand={stage['k_lig']}")

            # Minimisation loop with chunking
            minim_chunk = 5 # how often to output updates
            total_iter = stage["iter"]
            tolerance = stage.get("tolerance", 10) * kilojoule/(mole*nanometer)
            with tqdm(total=total_iter, desc=f"Stage {i} minimisation", unit="iter") as pbar:
                steps_done = 0
                while steps_done < total_iter:
                    n = min(minim_chunk, total_iter - steps_done)
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

            logger.debug(
                f"Worst restrained atom: {atom.name} in residue {atom.residue.name} "
                f"(resid {atom.residue.id}, chain {atom.residue.chain.id})"
            )
            logger.debug(f"Max restrained force: {max_force:.2f} kJ/mol/nm")

            # Warn if unrestrained atoms have very high forces
            unrestrained_indices = np.setdiff1d(np.arange(len(forces)), restrained_atoms)
            if len(unrestrained_indices) > 0:
                unres_forces = force_norms[unrestrained_indices]
                if np.max(unres_forces) > 50000:
                    logger.warning(f"High forces on unrestrained atoms: max {np.max(unres_forces):.1f} kJ/mol/nm")


        # Remove restraints

        self.system.removeForce(lig_idx)
        self.system.removeForce(prot_idx)
        self.simulation.context.reinitialize(preserveState=True)

        # Final check

        state = self.simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
        forces = state.getForces(asNumpy=True) / (kilojoule/mole/nanometer)
        max_force = np.max(np.linalg.norm(forces, axis=1))
        logger.info(f"Final max force: {max_force:.2f} kJ/mol/nm")
        if max_force > 50000:
            logger.warning(f"Forces are too high - re-prepare your system.")

        output_dir = self.config.get("heat_and_equilibrate", {}).get("output_dir", ".")
        os.makedirs(output_dir, exist_ok=True)

        minimised_pdb_path = os.path.join(output_dir, "minimised_system.pdb")
        PDBFile.writeFile(self.topology, state.getPositions(), open(minimised_pdb_path, "w"))
        
        logger.info(f"Minimised PDB written: {minimised_pdb_path}")

    @register_task("heat_and_equilibrate", category='Molecular dynamics',
               description="Heating and equilibration.")
    def heat_and_equilibrate(self):

        ligand_resname = 'LIG'
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
        # 1. HEATING PHASE
        # -------------------------------
        heating_cfg = self.config.get("heat_and_equilibrate", {}).get("heating", {})
        num_rounds = heating_cfg.get("num_steps", 8)
        steps_per_round = heating_cfg.get("steps_per_round", 5000)
        target_temp = heating_cfg.get("target_temp", 300)
        restraints = heating_cfg.get("restraint_strengths", [2.0, 1.5, 1.0, 0.5, 0.25, 0.1, 0.05, 0.0])

        logger.info(f"Heating with {num_rounds} rounds up to {target_temp} K")

        # backbone + ligand atoms to restrain
        backbone_atoms = [atom.index for atom in atoms_list
                        if atom.residue.name in protein_resnames and atom.name in {"N","CA","C","O"}]
        restrain_atoms = backbone_atoms + ligand_atoms

        # -------------------------------
        # Create restraint force (one reinit)
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

        # Debug: check initial displacement before adding to system TODO: this is broken
        max_disp = 0.0
        for i in range(restraint_force.getNumParticles()):
            particle_index, params = restraint_force.getParticleParameters(i)
            x0, y0, z0 = params
            pos = positions[particle_index]
            disp = np.linalg.norm([pos.x - x0, pos.y - y0, pos.z - z0])
            max_disp = max(max_disp, disp)
        logger.debug(f"Max initial displacement before adding restraint: {max_disp:.6f} nm")

        # Add to system and reinitialise context once
        restraint_index = self.system.addForce(restraint_force)
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
        logger.debug(f"Max displacement after reinit: {max_disp_post:.6f} nm")

        # Debug: check initial forces
        state = self.simulation.context.getState(getForces=True)
        forces = state.getForces(asNumpy=True) / (kilojoule / mole / nanometer)
        # max force on restrained atoms only
        restrained_indices = [restraint_force.getParticleParameters(i)[0] for i in range(restraint_force.getNumParticles())]
        max_force = np.max(np.linalg.norm(forces[restrained_indices], axis=1))
        logger.debug(f"Max force on restrained atoms before stepping: {max_force:.2f} kJ/mol/nm")

        # -------------------------------
        # Heating loop
        # -------------------------------
        original_dt = self.integrator.getStepSize()
        small_dt = 0.00025 * picoseconds  # 0.25 fs initially
        heat_step_chunk = 50 # how often to update and record data (too small does slow this down) TODO: tweak it a bit

        heating_data = {
            "step": [],
            "temperature": [],
            "potential_energy": [],
        }

        for i in range(num_rounds):
            temp = 10 + (target_temp - 10) * (i + 1) / num_rounds
            k_value = restraints[i] if i < len(restraints) else 0.0
            self.integrator.setTemperature(temp * kelvin)
            self.simulation.context.setParameter("k_heating", k_value * kilojoule/(mole*nanometer**2))
            self.integrator.setStepSize(small_dt if i < 3 else original_dt)

            steps_done = 0
            with tqdm(total=steps_per_round, desc=f"Heating round {i+1}", unit="steps") as pbar:
                while steps_done < steps_per_round:
                    n = min(heat_step_chunk, steps_per_round - steps_done)
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
                            logger.error(f"Atom {atom.name} (res {atom.residue.name}, idx {idx}) force: {fval:.2f} kJ/mol/nm")
                        raise e

                    steps_done += n
                    pbar.update(n)

                    # Force debugs in eahc stage
                    state = self.simulation.context.getState(getForces=True)
                    forces = state.getForces(asNumpy=True) / (kilojoule/mole/nanometer)
                    max_force = np.max(np.linalg.norm(forces, axis=1))
                    if max_force > 3000:
                        top5_idx = np.argsort(np.linalg.norm(forces, axis=1))[-5:]
                        logger.debug(f"High force detected: {max_force:.2f} kJ/mol/nm at step {steps_done}")
                        for idx in top5_idx:
                            atom = atoms_list[idx]
                            fval = np.linalg.norm(forces[idx])
                            logger.debug(f"Atom {atom.name} (res {atom.residue.name}, idx {idx}) force: {fval:.2f} kJ/mol/nm")
                    
                    # track diagnostics for plotting

                    state = self.simulation.context.getState(getEnergy=True)
                    temp_inst = self.integrator.getTemperature().value_in_unit(kelvin)
                    pot_energy = state.getPotentialEnergy().value_in_unit(kilojoule/mole)

                    heating_data["step"].append(steps_done + i * steps_per_round)
                    heating_data["temperature"].append(temp_inst)
                    heating_data["potential_energy"].append(pot_energy)


        # Remove heating restraint

        state = self.simulation.context.getState(getPositions=True)

        self.system.removeForce(restraint_index)
        self.simulation.context.reinitialize(preserveState=True)
        # self.simulation.context.setPositions(state.getPositions())
        self.simulation.context.setVelocitiesToTemperature(target_temp * kelvin)

        logger.info("Heating stage complete")

        # -------------------------------
        # 2. EQUILIBRATION PHASE
        # -------------------------------

        equil_cfg = self.config.get("heat_and_equilibrate", {}).get("equilibration", {})
        ensemble = equil_cfg.get("ensemble", "NPT").upper()
        steps = equil_cfg.get("steps", 100000)
        equil_temp = equil_cfg.get("target_temp", target_temp)

        logger.info(f"Starting equilibration: {ensemble} ensemble for {steps} steps at {equil_temp} K")

        # Add barostat
        if ensemble == "NPT" and not any(isinstance(self.system.getForce(i), MonteCarloBarostat) 
                                        for i in range(self.system.getNumForces())):
            pressure = equil_cfg.get("pressure", 1.0) * atmosphere
            frequency = equil_cfg.get("barostat", {}).get("frequency", 25)
            barostat = MonteCarloBarostat(pressure, equil_temp * kelvin, frequency)
            self.system.addForce(barostat)
            self.simulation.context.reinitialize(preserveState=True)
            # self.simulation.context.setPositions(state.getPositions())
            self.simulation.context.setVelocitiesToTemperature(equil_temp * kelvin)
            logger.info(f"Added MonteCarloBarostat: {pressure} atm, {equil_temp} K")

        self.integrator.setTemperature(equil_temp * kelvin)

        lipid_resnames = {"POP"}  # should match topology

        lipid_residues = [
            res for res in self.topology.residues()
            if res.name in lipid_resnames
        ]

        n_lipids = len(lipid_residues)

        if n_lipids == 0:
            logger.warning("No lipids detected - APL will be NaN")

        logger.debug(f"Detected {n_lipids} lipids")

        # Prepare DataFrame
        columns = ["step", "temperature", "potential_energy", "volume", "density",
                "lx", "ly", "lz", "apl", "anisotropy"]
        equil_df = pd.DataFrame(columns=columns)

        # Equilibration loop
        equil_step_chunk = 10
        steps_done = 0

        with tqdm(total=steps, desc="Equilibration", unit="steps") as pbar:
            while steps_done < steps:
                n = min(equil_step_chunk, steps - steps_done)
                self.simulation.step(n)
                steps_done += n
                pbar.update(n)

                # -------------------------
                # Diagnostics
                # -------------------------
                state = self.simulation.context.getState(getEnergy=True, getPositions=True)

                kinetic_energy = state.getKineticEnergy()

                dof = 3 * self.system.getNumParticles() - self.system.getNumConstraints()
                dof -= 3  # remove com
                kB = 0.00831446261815324  # kJ/mol/K

                temp_inst = (2 * kinetic_energy.value_in_unit(kilojoule/mole)) / (dof * kB)
                pot_energy = state.getPotentialEnergy().value_in_unit(kilojoule/mole)

                box = state.getPeriodicBoxVectors(asNumpy=True)

                lx = box[0][0].value_in_unit(nanometer)
                ly = box[1][1].value_in_unit(nanometer)
                lz = box[2][2].value_in_unit(nanometer)

                volume_nm3 = lx * ly * lz

                # APL
                if n_lipids > 0:
                    apl = (lx * ly) / (n_lipids / 2.0)
                else:
                    apl = float("nan")

                mass = sum([self.system.getParticleMass(i) for i in range(self.system.getNumParticles())])
                mass_g = mass.value_in_unit(dalton) * 1.66054e-24
                volume_mL = volume_nm3 * 1e-21
                density = mass_g / volume_mL

                anisotropy = lz / ((lx + ly) / 2.0)

                # Append all diags to DataFrame
                equil_df.loc[len(equil_df)] = [steps_done, temp_inst, pot_energy, volume_nm3, density,
                                        lx, ly, lz, apl, anisotropy]

        # Diagnostics
        state = self.simulation.context.getState(getEnergy=True, getPositions=True)
        box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
        volume_nm3 = (box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]).value_in_unit(nanometer**3)
        mass = sum([self.system.getParticleMass(i) for i in range(self.system.getNumParticles())])
        mass_g = mass.value_in_unit(dalton) * 1.66054e-24
        volume_mL = volume_nm3 * 1e-21
        density = mass_g / volume_mL

        logger.info(f"Equilibration stage complete. Box volume: {volume_nm3:.3f} nm³, density: {density:.3f} g/mL")
        
        if abs(density - 1.0) > 0.1:
            logger.warning("Density deviates from expected ~1.0 g/mL. Check for bubbles or equilibration issues.")

        # Save plots

        output_dir = (
            self.config.get("heat_and_equilibrate", {}).get("output_dir")
        )

        # -------------------------------
        # Save diagnostic plots
        # -------------------------------
        output_dir = self.config.get("heat_and_equilibrate", {}).get("output_dir")
        plot_dir = os.path.join(output_dir, "heat_equil_diagnostic_plots")
        os.makedirs(plot_dir, exist_ok=True)
        logger.info(f"Saving heat/equil diagnostic plots to: {plot_dir}")

        plots = [
            # Heating
            {"x": heating_data["step"], "y": heating_data["temperature"], "ylabel": "Temperature (K)", "title": "Heating Temperature", "file": "heating_temperature.png"},
            {"x": heating_data["step"], "y": heating_data["potential_energy"], "ylabel": "Potential Energy (kJ/mol)", "title": "Heating Energy", "file": "heating_energy.png"},
            # Equilibration
            {"x": equil_df["step"], "y": equil_df["temperature"], "ylabel": "Temperature (K)", "title": "Equilibration Temperature", "file": "equil_temperature.png"},
            {"x": equil_df["step"], "y": equil_df["potential_energy"], "ylabel": "Potential Energy (kJ/mol)", "title": "Equilibration Energy", "file": "equil_energy.png"},
            {"x": equil_df["step"], "y": equil_df["density"], "ylabel": "Density (g/mL)", "title": "Equilibration Density", "file": "equil_density.png"},
            {"x": equil_df["step"], "y": equil_df["volume"], "ylabel": "Volume (nm³)", "title": "Equilibration Volume", "file": "equil_volume.png"},
            {"x": equil_df["step"], "y": {"Lx": equil_df["lx"], "Ly": equil_df["ly"], "Lz": equil_df["lz"]}, "ylabel": "Box length (nm)", "title": "Box Dimensions", "file": "equil_box_dimensions.png", "labels": True},
            {"x": equil_df["step"], "y": equil_df["apl"], "ylabel": "Area per lipid (nm²)", "title": "Area per Lipid (APL)", "file": "equil_apl.png"},
            {"x": equil_df["step"], "y": equil_df["anisotropy"], "ylabel": "Z / XY ratio", "title": "Box Anisotropy", "file": "equil_anisotropy.png"},
        ]

        # save plots
        for p in plots:
            save_plot(
                x=p["x"],
                y=p["y"],
                xlabel="Step",
                ylabel=p["ylabel"],
                title=p["title"],
                filepath=os.path.join(plot_dir, p["file"]),
                labels=p.get("labels", False),
            )
        
        # output equilibrated pdb
        
        equilibrated_pdb_path = os.path.join(output_dir, "equilibrated_system.pdb")
        # psf_path = os.path.join(output_dir, "topology.psf")

        # with open(psf_path, "w") as f:
        #     CharmmPsfFile.writeFile(self.topology, f)

        # logger.info(f"PSF written: {psf_path}")

        # Write raw coordinates
        raw_pdb = os.path.join(output_dir, "_temp_raw.pdb")
        state = self.simulation.context.getState(getPositions=True, enforcePeriodicBox=False)

        with open(raw_pdb, "w") as f:
            PDBFile.writeFile(self.topology, state.getPositions(), f)

        # Load with correct bonding
        with open(equilibrated_pdb_path, "w") as f:
            PDBFile.writeFile(self.topology, state.getPositions(), f)

        os.remove(raw_pdb)

        logger.info(f"Equilibrated PDB written: {equilibrated_pdb_path}")
    
    def reimage_trajectory(self, traj_path, top_path, output_path):

        logger.info(f"Re-imaging trajectory: {traj_path}")

        u = mda.Universe(top_path, traj_path)

        membrane = u.select_atoms("resname POP")
        protein = u.select_atoms("protein")

        if len(membrane) > 0:
            anchor = membrane
            logger.info("Centering on membrane")
        else:
            anchor = protein
            logger.warning("No membrane found, centering on protein")

        transformations = [
            unwrap(u.atoms),  # required for your version
            center_in_box(anchor, wrap=False),
            wrap(u.atoms)
        ]

        u.trajectory.add_transformations(*transformations)

        with mda.Writer(output_path, n_atoms=u.atoms.n_atoms) as W:
            for ts in u.trajectory:
                W.write(u.atoms)

        logger.info(f"Saved wrapped trajectory: {output_path}")
        
    @register_task("production", category='Molecular dynamics', description="Run production simulation.")
    def production(self):

        cfg = self.config["production"]
        ns = cfg.get("length_ns", 1)
        ensemble = cfg.get("ensemble", "NPT").upper()
        target_temp = cfg.get("target_temp", self.integrator.getTemperature().value_in_unit(kelvin))
        pressure = cfg.get("pressure", 1.0)

        # Thresholds for event triggers
        max_temp_dev = 100        # K deviation from target
        max_density_dev = 0.5     # g/mL deviation

        # Get current state
        state = self.simulation.context.getState(getPositions=True, getVelocities=True)
        positions = state.getPositions()
        velocities = state.getVelocities()
        box_vectors = state.getPeriodicBoxVectors()

        # #### Center system in the box
        # coords = np.array([
        #     [
        #         p.x.value_in_unit(nanometer) if hasattr(p.x, "value_in_unit") else float(p.x),
        #         p.y.value_in_unit(nanometer) if hasattr(p.y, "value_in_unit") else float(p.y),
        #         p.z.value_in_unit(nanometer) if hasattr(p.z, "value_in_unit") else float(p.z)
        #     ]
        #     for p in positions
        # ], dtype=float)
        
        # anchor_atoms = [
        #     atom.index for atom in self.topology.atoms()
        #     if atom.residue.name in {"POP"} or atom.residue.name == "UNK" or atom.residue.chain.id == "A"
        # ]
        
        # com = np.mean(coords[anchor_atoms], axis=0)
        # box_lengths = np.array([box_vectors[i][i].value_in_unit(nanometer) for i in range(3)])
        # shift = box_lengths / 2 - com
        # for i in range(len(positions)):
        #     x, y, z = coords[i] + shift
        #     positions[i] = Vec3(x, y, z) * nanometer

        # Update context
        self.simulation.context.setPositions(positions)
        self.simulation.context.setVelocities(velocities)
        self.simulation.context.setPeriodicBoxVectors(*box_vectors)
        self.integrator.setTemperature(target_temp * kelvin)
        self.simulation.currentStep = 0

        # Attach barostat if NPT
        has_barostat = any(isinstance(f, MonteCarloBarostat) for f in
                        [self.system.getForce(i) for i in range(self.system.getNumForces())])
        if ensemble == "NPT" and not has_barostat:
            barostat = MonteCarloBarostat(pressure * atmosphere, target_temp * kelvin, 25)
            self.system.addForce(barostat)
            self.simulation.context.reinitialize(preserveState=True)
            self.simulation.context.setPositions(positions)
            self.simulation.context.setVelocities(velocities)
            self.simulation.context.setPeriodicBoxVectors(*box_vectors)
            logger.info(f"Added MonteCarloBarostat: {pressure} atm, {target_temp} K")

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

        # set how often to output data to log
        log_multiplier = 1 # set to 10, 100, etc. to run expensive evaluations that are dumped to log less frequently.
        energy_log_freq = output_freq * log_multiplier

        logger.info(f"Running PROD simulation for {ns} ns ({total_steps} steps)")
        logger.info(f"Split traj every {split_ns} ns ({n_segments} segments).")
        logger.info(f"Writing frames every {output_freq} steps.")

        # Rolling diagnostic averages to check simulation is running smoothly
        rolling_window = 10  # number of diagnostic points to average over
        temp_history = []
        density_history = []
        energy_history = []

        # Segment loop
        for seg in range(1, n_segments + 1):
            seg_start_step = self.simulation.currentStep
            seg_end_step = min(seg_start_step + split_steps, total_steps)
            seg_steps = seg_end_step - seg_start_step

            traj_file = os.path.join(output_dir, f"trajectory_{seg:03d}.dcd")
            log_file = os.path.join(output_dir, f"log_{seg:03d}.txt")
            chk_file = os.path.join(output_dir, f"restart_{seg:03d}.chk")

            # Attach reporters
            self.simulation.reporters.clear()
            self.simulation.reporters.append(DCDReporter(traj_file, output_freq, enforcePeriodicBox=True))
            self.simulation.reporters.append(StateDataReporter(
                log_file, energy_log_freq,
                step=True, temperature=True, potentialEnergy=True, totalEnergy=True,
                density=True, progress=True, remainingTime=True, speed=True,
                totalSteps=seg_steps, separator='\t'
            ))

            logger.info(f"Segment {seg}/{n_segments}: {seg_steps} steps → {traj_file}")

            with tqdm(total=seg_steps, desc=f"Production seg {seg}", unit="steps") as pbar:
                current_step = 0
                diagnostic_freq = energy_log_freq  # diagnostics check frequency
                while current_step < seg_steps:
                    n = min(output_freq, seg_steps - current_step)
                    self.simulation.step(n)
                    current_step += n
                    pbar.update(n)

                    # Diagnostics only at intervals
                    if current_step % diagnostic_freq == 0:
                        state = self.simulation.context.getState(getEnergy=True)

                        # Potential energy check
                        potential_energy = state.getPotentialEnergy().value_in_unit(kilojoule/mole)
                        if potential_energy > 1e6:  # adjust threshold as needed
                            logger.warning(f"High potential energy detected: {potential_energy:.1f} kJ/mol at step {self.simulation.currentStep}")

                        # Temperature and density
                        kinetic_energy = state.getKineticEnergy().value_in_unit(kilojoule/mole)
                        dof = 3 * self.system.getNumParticles() - self.system.getNumConstraints()
                        dof -= 3  # remove com
                        kB = 0.00831446261815324
                        temp_inst = 2 * kinetic_energy / (dof * kB)
                        box_vectors = state.getPeriodicBoxVectors()
                        lx, ly, lz = [box_vectors[i][i].value_in_unit(nanometer) for i in range(3)]
                        volume = lx * ly * lz
                        mass = sum([self.system.getParticleMass(i) for i in range(self.system.getNumParticles())]).value_in_unit(dalton) * 1.66054e-24
                        density = mass / (volume * 1e-21)

                        if abs(temp_inst - target_temp) > max_temp_dev:
                            logger.warning(f"Temperature deviation: {temp_inst:.1f} K (target={target_temp} K)")
                        if abs(density - 1.0) > max_density_dev:
                            logger.warning(f"Density deviation: {density:.2f} g/mL")

            # Save checkpoint at end of segment
            self.simulation.saveCheckpoint(chk_file)
            logger.info(f"Checkpoint saved: {chk_file}")

            # Re-image trajectory after writing
            wrapped_traj = traj_file.replace(".dcd", "_wrapped.xtc")

            try:
                self.reimage_trajectory(
                    traj_path=traj_file,
                    top_path=os.path.join(output_dir, "equilibrated_system.pdb"),
                    output_path=wrapped_traj
                )
                logger.info(f"Wrapped trajectory written: {wrapped_traj}")
            except Exception as e:
                logger.warning(f"Re-imaging failed for {traj_file}: {e}")

        # Save final checkpoint
        final_chk = os.path.join(output_dir, "restart_final.chk")
        self.simulation.saveCheckpoint(final_chk)
        logger.info(f"Final checkpoint saved: {final_chk}")

        # -------------------------------
        # Save provenance / metadata
        # -------------------------------
        prov_file = os.path.join(output_dir, "provenance.json")
        metadata = {
            "ensemble": ensemble,
            "temperature": target_temp,
            "pressure": pressure,
            "length_ns": ns,
            "timestep_fs": timestep.value_in_unit(femtoseconds),
            "split_ns": split_ns,
            "output_frequency": output_freq,
            "software": {
                "python_version": platform.python_version(),
                "openmm_version": sys.modules['openmm'].__version__,
                "platform": platform.platform()
            },
            "system": {
                "num_particles": self.system.getNumParticles()
            }
        }
        with open(prov_file, "w") as f:
            json.dump(metadata, f, indent=2)
        logger.info(f"Simulation metadata saved: {prov_file}")

        logger.info("Production simulation complete.")


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





