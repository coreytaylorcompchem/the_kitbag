# from pipeline.logger import setup_logger

# logger = setup_logger(__name__, debug_mode=True, simple_format=True)

import numpy as np

from MDAnalysis.lib.distances import distance_array

class ComponentDetector:
    """
    Detects molecular components in a system:
    - receptor (largest protein)
    - partner (second protein, e.g. VHH)
    - ligand (optional small molecule)
    """

    def __init__(self, universe, ligand_resname="UNK", water_resname="HOH"):
        self.u = universe
        self.ligand_resname = ligand_resname
        self.water_resname = water_resname

    def detect(self):
        protein = self.u.select_atoms("protein")

        # STEP 1: cluster protein into connected components

        coords = protein.positions

        # distance cutoff for connectivity (Å)
        cutoff = 4.5

        # build adjacency
        dist = distance_array(coords, coords)
        adjacency = dist < cutoff

        # simple DFS to find connected components
        visited = set()
        clusters = []

        for i in range(len(coords)):
            if i in visited:
                continue

            stack = [i]
            group = []

            while stack:
                j = stack.pop()
                if j in visited:
                    continue
                visited.add(j)
                group.append(j)

                neighbors = list(np.where(adjacency[j])[0])
                stack.extend(neighbors)

            clusters.append(protein[group])

        # STEP 2: build components
        components = [{
            "atoms": cluster,
            "n_residues": len(cluster.residues)
        } for cluster in clusters]

        components = sorted(components, key=lambda x: x["n_residues"], reverse=True)

        # STEP 3: decide if partner exists
        receptor = components[0]["atoms"] if components else None
        partner = None

        if len(components) > 1:
            # Heuristic: second component must be "large enough"
            size_ratio = components[1]["n_residues"] / components[0]["n_residues"]

            # Only call it a partner if it's substantial (e.g. VHH)
            if size_ratio > 0.2:
                partner = components[1]["atoms"]

        # STEP 4: ligand detection
        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")

        if len(ligand) == 0:
            ligand = self.u.select_atoms(
                f"not protein and not resname {self.water_resname} and not ions"
            )

        if len(ligand) == 0:
            ligand = None

        return {
            "receptor": receptor,
            "partner": partner,
            "ligand": ligand,
            "has_partner": partner is not None,
            "has_ligand": ligand is not None
        }