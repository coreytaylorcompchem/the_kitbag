import MDAnalysis as mda
import networkx as nx
import tqdm
# from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from typing import List
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


class HydrogenBondAnalysisTask:
    """
    Task to calculate both direct and water-mediated hydrogen bonds
    between a binding site and a ligand.
    """

    def __init__(self,
                 topology: str,
                 trajectory: str,
                 binding_site_resids: List[int],
                 ligand_resname: str,
                 start: int = 0,
                 stop: int = -1,
                 step: int = 1,
                 water_resname: str = "HOH"):
        self.topology = topology
        self.trajectory = trajectory
        self.binding_site_resids = binding_site_resids
        self.ligand_resname = ligand_resname
        self.start = start
        self.stop = stop
        self.step = step
        self.water_resname = water_resname

        # Load universe
        self.u = mda.Universe(self.topology, self.trajectory)
        self.direct_paths = []
        self.water_paths = []

    def _compute_hbonds(self, water_mediated: bool = False):
        paths = []
        binding_site_atoms = self.u.select_atoms(
            "protein and resid " + " ".join(str(r) for r in self.binding_site_resids)
        )

        for ts in tqdm.tqdm(self.u.trajectory[self.start:self.stop:self.step],
                            desc="Direct HBonds" if not water_mediated else "Water-mediated HBonds"):
            if water_mediated:
                donors_sel = f"prop mass > 13 and resname {self.ligand_resname} " \
                             f"or (prop mass > 13 and resname {self.water_resname} and around 10 resname {self.ligand_resname}) " \
                             f"or (prop mass > 13 and protein and index {' '.join(map(str, binding_site_atoms.indices))})"
                hydrogens_sel = f"prop mass < 2 and resname {self.ligand_resname} " \
                                f"or (prop mass < 2 and resname {self.water_resname} and around 10 resname {self.ligand_resname}) " \
                                f"or (prop mass < 2 and protein and index {' '.join(map(str, binding_site_atoms.indices))})"
                acceptors_sel = donors_sel
                between = [f"protein or resname {self.ligand_resname}", f"resname {self.water_resname}"]
                path_length = 3
            else:
                donors_sel = f"prop mass > 13 and resname {self.ligand_resname} " \
                             f"or (prop mass > 13 and protein and index {' '.join(map(str, binding_site_atoms.indices))})"
                hydrogens_sel = f"prop mass < 2 and resname {self.ligand_resname} " \
                                f"or (prop mass < 2 and protein and index {' '.join(map(str, binding_site_atoms.indices))})"
                acceptors_sel = donors_sel
                between = ["protein", f"resname {self.ligand_resname}"]
                path_length = 2

            hba = HBA(
                universe=self.u,
                donors_sel=donors_sel,
                hydrogens_sel=hydrogens_sel,
                acceptors_sel=acceptors_sel,
                between=between,
                d_h_cutoff=1.2,
                d_h_a_angle_cutoff=120,
                d_a_cutoff=3.0
            )
            hba.run(start=ts.frame, stop=ts.frame + 1)
            unique_hbs = [list(i) for i in hba.hbonds[:, 1:-2] if i[0] != i[1]]

            G = nx.Graph()
            for edge in unique_hbs:
                donor = f"{self.u.atoms[int(edge[0])].resname}{self.u.atoms[int(edge[0])].resid}"
                acceptor = f"{self.u.atoms[int(edge[2])].resname}{self.u.atoms[int(edge[2])].resid}"
                if self.u.atoms[int(edge[0])].resname != self.water_resname:
                    donor += f": {self.u.atoms[int(edge[0])].name}"
                if self.u.atoms[int(edge[2])].resname != self.water_resname:
                    acceptor += f": {self.u.atoms[int(edge[2])].name}"
                G.add_edge(donor, acceptor)

            sources = [n for n in G.nodes if self.ligand_resname not in n and self.water_resname not in n]
            targets = [n for n in G.nodes if self.ligand_resname in n]

            for target in targets:
                if G.has_node(target):
                    for s in sources:
                        if nx.has_path(G, source=s, target=target):
                            paths_iter = nx.all_shortest_paths(G, source=s, target=target)
                            for p in paths_iter:
                                if len(p) == path_length:
                                    paths.append([list(p), ts.frame])
        return paths

    @register_task(
        "solvent_hbonds",
        category="Analyses",
        description="Compute direct and water-mediated hydrogen bonds."
    )
    def run(self):
        logger.info("Calculating direct hydrogen bonds...")
        self.direct_paths = self._compute_hbonds(water_mediated=False)
        logger.info("Calculating water-mediated hydrogen bonds...")
        self.water_paths = self._compute_hbonds(water_mediated=True)
        logger.info(f"HBond analysis complete. Frames analyzed: {len(self.u.trajectory[self.start:self.stop:self.step])}")
        return {"direct": self.direct_paths, "water_mediated": self.water_paths}
