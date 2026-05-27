from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

import numpy as np
import MDAnalysis as mda

class ComponentDetector:

    def __init__(
        self,
        universe,
        ligand_resname="UNK",
        water_resname="HOH",
        ion_resnames=None
    ):
        self.u = universe
        self.ligand_resname = ligand_resname
        self.water_resname = water_resname
        self.ion_resnames = ion_resnames or [
            "NA", "CL", "K", "MG", "CA"
        ]

    def detect(self):

        protein = self.u.select_atoms("protein")

        if len(protein) == 0:
            raise ValueError(
                "No protein atoms found in universe!"
            )

        protein_segids = sorted(set(protein.segids))

        logger.debug(
            f"Detected protein segids: {protein_segids}"
        )

        ligand = self.u.select_atoms(
            f"resname {self.ligand_resname}"
        )

        if len(ligand) == 0:

            ions_str = " ".join(self.ion_resnames)

            lipid_resnames = [
                "POP", "POPC", "POPE",
                "POPG", "DOPC", "DPPC",
                "CHL1"
            ]

            lipid_str = " ".join(lipid_resnames)

            ligand = self.u.select_atoms(
                f"""
                not protein
                and not resname {self.water_resname}
                and not resname {ions_str}
                and not resname {lipid_str}
                """
            )

            if len(ligand) == 0:
                ligand = None

        return {
            "protein_segids": protein_segids,
            "ligand": ligand,
            "ligand_resnames": (
                sorted(set(ligand.resnames))
                if ligand is not None else []
            ),
            "has_ligand": ligand is not None
        }