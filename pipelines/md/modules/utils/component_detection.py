# from pipeline.logger import setup_logger

# logger = setup_logger(__name__, debug_mode=True, simple_format=True)

import numpy as np
import MDAnalysis as mda

class ComponentDetector:
    """
    Detects molecular components in a system:
    - receptor (largest protein chain)
    - partner (second protein chain, e.g. VHH)
    - ligand (optional small molecule)
    """

    def __init__(self, universe, ligand_resname="UNK", water_resname="HOH", ion_resnames=None):
        self.u = universe
        self.ligand_resname = ligand_resname
        self.water_resname = water_resname
        self.ion_resnames = ion_resnames or ["NA", "CL", "K", "MG", "CA"]

    def detect(self):
        protein = self.u.select_atoms("protein")
        if len(protein) == 0:
            raise ValueError("No protein atoms found in universe!")

        # -----------------------------
        # STEP 1: detect receptor and partner by chain
        # -----------------------------
        chains = sorted(set(atom.segid for atom in protein))
        receptor = protein.select_atoms(f"segid {chains[0]}") if chains else None
        partner = protein.select_atoms(f"segid {chains[1]}") if len(chains) > 1 else None

        # -----------------------------
        # STEP 2: ligand detection
        # -----------------------------
        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")
        if len(ligand) == 0:
            # fallback: select any non-protein, non-water, non-ion atoms
            ions_str = " ".join(self.ion_resnames)
            ligand = self.u.select_atoms(
                f"not protein and not resname {self.water_resname} and not resname {ions_str}"
            )
            if len(ligand) == 0:
                ligand = None

        # -----------------------------
        # STEP 3: return components
        # -----------------------------
        return {
            "receptor": receptor,
            "partner": partner,
            "ligand": ligand,
            "has_partner": partner is not None,
            "has_ligand": ligand is not None
        }