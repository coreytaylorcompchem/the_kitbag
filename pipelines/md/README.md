# Molecular dynamics pipeline

Perform different types of MD in automated, end-to-end fashion. 

![Molecular dynamics logo](_images/md_logo.png)

## Setup and dependencies

* Mainly Python libs
* [Cuda](https://developer.nvidia.com/cuda-toolkit)
* [Rdkit](https://www.rdkit.org/docs/Install.html)
* [OpenMM](https://docs.openmm.org/latest/userguide/application/01_getting_started.html)
* [PDBFixer](https://github.com/openmm/pdbfixer)
* [Biopython](https://biopython.org/) (may remove this)

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python ../../md_pipeline.py --params ../../configs/md.yaml`

If the calculation runs correctly, you should see outputs like logs, trajectories, etc. appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
Available Tasks:
  [Molecular dynamics]:
    - heat_and_equilibrate: Heating and equilibration.
    - minimize: Initial energy minimization.
    - prepare_system: Load inputs, cap chains, parameterise ligand and protein, solvate and save final topology.
    - production: Run production simulation.
    - setup_simulation: Set up integrator and simulation.
  [Post-proc; graph analyses]:
    - hydration_site_energy: Identify hydration sites and rank them by approximate free energy.
    - network_embedding_analysis: Convert trajectory frames into residue–ligand contact graphs, then perform Node2Vec + t-SNE embedding to visualise network evolution.
    - protein_ligand_communities: Detect cooperative residue clusters (communities) in the protein–ligand interaction network.
    - solvent_hbonds: Compute direct and water-mediated hbonds with ligand.
    - temporal_motif_persistence: Quantify persistence of small recurring motifs (e.g., ligand–water–residue triangles) in the solvent network.
  [Post-proc; traj analyses]:
    - interaction_fingerprint: Compute IFP using ProLIF, generate barcode plot.
    - rmsd_analysis: Compute RMSD of protein bb and ligand.
    - rmsf_analysis: Compute per-residue RMSF for protein (all atoms and Cα only).
Available Workflows:
 - molecular_dynamics: Perform molecular dynamics simulation using specified backend.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* More options for MD simulations.
    * Different constraints (harmonic, restraints, etc.)
    * Control outputs.
    * Restarts.
    * Classic and ML force fields.
    * Loop modelling.
* Automated analyses.
    * Basic trajectory analysis (RMSD/F, hbonds, torsions, solvent lifetimes, etc.)
    * Advanced trajectory anlyses (PCA, water analyses, graph-based workflows). 
* More advanced simulation types.
    * QMMM.
    * Steered MD.
    * Metadynamics.