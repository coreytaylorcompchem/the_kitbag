from dataclasses import dataclass
from typing import Iterable, Set
import numpy as np

import openmm.unit as unit


@dataclass
class FlexibleShellSelector:
    """
    Selects protein atoms for induced-fit refinement based on distance
    to ligand atoms and residue-level constraints.
    """
    cutoff_angstrom: float
    residue_select: Iterable[str] = ()
    backbone_refinement: bool = False
    ligand_resname: str = "LIG"

    BACKBONE_ATOMS = {"N", "CA", "C", "O"}

    def select_flexible_atoms(self, topology, positions) -> Set[int]:
        """
        Returns atom indices that are allowed to move (not restrained).
        """
        cutoff_nm = self.cutoff_angstrom * 0.1
        residue_select = set(self.residue_select)

        # Convert positions → numpy (nm)
        pos = np.array([p.value_in_unit(unit.nanometer) for p in positions])

        ligand_atoms = [
            a for a in topology.atoms()
            if a.residue.name == self.ligand_resname
        ]
        if not ligand_atoms:
            raise RuntimeError("Ligand atoms not found in topology")

        ligand_indices = [a.index for a in ligand_atoms]

        flexible_atoms = set()

        for atom in topology.atoms():
            res = atom.residue

            # Skip ligand itself
            if res.name == self.ligand_resname:
                continue

            # Optional residue-name filtering
            if residue_select and res.name not in residue_select:
                continue

            # Backbone handling
            if not self.backbone_refinement and atom.name in self.BACKBONE_ATOMS:
                continue

            atom_pos = pos[atom.index]

            # Minimum distance to ligand
            min_dist = min(
                np.linalg.norm(atom_pos - pos[i]) for i in ligand_indices
            )

            if min_dist <= cutoff_nm:
                flexible_atoms.add(atom.index)

        return flexible_atoms

    def select_restrained_atoms(self, topology, positions) -> Set[int]:
        """
        Convenience method: returns atoms to restrain.
        """
        flexible = self.select_flexible_atoms(topology, positions)
        all_atoms = {a.index for a in topology.atoms()}
        return all_atoms - flexible
