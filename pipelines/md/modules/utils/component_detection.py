from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

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
        # Build chain → atomgroup mapping
        chain_groups = {
            seg: protein.select_atoms(f"segid {seg}")
            for seg in set(atom.segid for atom in protein)
        }

        # Debug: list all detected chains
        logger.debug("[ComponentDetector] Detected protein chains:")
        for seg, ag in chain_groups.items():
            logger.debug(f"  - segid {seg}: {len(ag.atoms)} atoms, {len(ag.residues)} residues")

        # Sort chains by size (number of atoms)
        sorted_chains = sorted(chain_groups.items(), key=lambda x: len(x[1]), reverse=True)

        logger.debug("[ComponentDetector] Chains sorted by size:")
        for i, (seg, ag) in enumerate(sorted_chains):
            logger.debug(f"  {i+1}. segid {seg}: {len(ag.residues)} residues")

        # Start with all chains as receptor candidates
        receptor_chains = [ag for _, ag in sorted_chains]

        partner = None

        logger.debug("[ComponentDetector] Evaluating partner candidate:")

        if len(sorted_chains) > 1:
            seg_second, second = sorted_chains[1]

            n_res_receptor = len(sorted_chains[0][1].residues)
            n_res_second = len(second.residues)

            ratio = n_res_second / n_res_receptor

            logger.debug(f"  Candidate segid {seg_second}")
            logger.debug(f"    residues: {n_res_second}")
            logger.debug(f"    ratio to receptor: {ratio:.3f}")

            is_vhh_like = 110 <= n_res_second <= 135
            is_reasonable_partner = (
                n_res_second > 80 and ratio > 0.02
            )

            logger.debug(f"    is_vhh_like: {is_vhh_like}")
            logger.debug(f"    is_reasonable_partner: {is_reasonable_partner}")

            if is_vhh_like or is_reasonable_partner:
                partner = second

                # REMOVE partner chain from receptor
                receptor_chains = [ag for ag in receptor_chains if ag is not second]

        # Merge all remaining chains into one AtomGroup
        receptor = sum(receptor_chains[1:], receptor_chains[0]) if receptor_chains else None

        partner = None

        # Heuristic: partner must be "substantial"
        if len(sorted_chains) > 1:
            second = sorted_chains[1][1]

            # ---- Heuristic thresholds ----
            min_atoms = 300        # avoids tiny fragments
            size_ratio = 0.2       # relative to receptor

            # Count residues instead of atoms
            n_res_receptor = len(receptor.residues)
            n_res_second = len(second.residues)

            ratio = n_res_second / n_res_receptor

            is_vhh_like = 110 <= n_res_second <= 135
            is_reasonable_partner = (
                n_res_second > 80 and ratio > 0.02
            )

            if is_vhh_like or is_reasonable_partner:
                partner = second

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

        logger.info("Automatically detecting system components")   
        logger.info("[ComponentDetector] Final assignment:")

        if receptor is not None:
            rec_segids = sorted(set(receptor.segids))
            logger.info(f"  Receptor chains: {rec_segids}")
            logger.info(f"  Total residues: {len(receptor.residues)}")

        if partner is not None:
            part_segids = sorted(set(partner.segids))
            logger.info(f"  Partner chains: {part_segids}")
            logger.info(f"  Residues: {len(partner.residues)}")
        else:
            logger.info("  No protein partner detected")

        if ligand is not None:
            lig_resnames = set(ligand.resnames)
            logger.info(f"  Ligand detected: {lig_resnames} ({len(ligand.atoms)} atoms)")
        else:
            logger.info("  No ligand detected")

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