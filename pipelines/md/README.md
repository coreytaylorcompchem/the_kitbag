# Molecular dynamics pipeline

Perform different types of MD in automated, end-to-end fashion. 

![Molecular dynamics logo](_images/md_logo.png)

## Setup and dependencies

* Mainly Python libs
* [OpenMM](https://docs.openmm.org/latest/userguide/application/01_getting_started.html)
* [OpenFE](https://github.com/OpenFreeEnergy/openfe)
* [MDAnalysis](https://www.mdanalysis.org/pages/installation_quick_start/)
* [Cuda](https://developer.nvidia.com/cuda-toolkit)
* [Rdkit](https://www.rdkit.org/docs/Install.html)
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
  [Free Energy]:
    - openfe_prepare_receptor: Prepare membrane receptor for OpenFE.
    - openfe_create_network: Create OpenFE ligand network for membrane RBFE.    
    - openfe_create_alchemical_network: Create OpenFE RBFE alchemical network.
    - openfe_run_network: Run all OpenFE transformations sequentially.
    - openfe_gather_results: Gather OpenFE results.
  [Molecular dynamics]:
    - prepare_system: Load inputs, cap chains, parameterise ligand and protein, solvate and save final topology.
    - setup_simulation: Set up integrator and simulation.
    - minimize: Configurable staged minimisation of protein systems.
    - heat_and_equilibrate: Configurable staged heating and NPT equilibration with diagnostic plots.
    - production: Run production NPT simulation.
  [Post-proc; free energy]:
    - mmpbsa: Compute MM/GB(PB)SA binding free energies at beginning and end of trajectory (WIP).
  [Post-proc; graph analyses]:
    - hydration_site_energy: Identify hydration sites and rank them by approximate free energy.
    - network_embedding_analysis: Network-based visualisation of protein–ligand contact graphs.
    - protein_ligand_communities: Detect cooperative residue clusters (communities) in the protein–ligand interaction network.
    - protein_protein_network_embedding: Network-based visualisation of protein–protein contact graphs.
    - solvent_hbonds: Compute direct and water-mediated hbonds with ligand.
    - temporal_motif_persistence: Detect persistence of small recurring solvent-mediated motifs.
  [Post-proc; kinetics]:
    - msm_analysis: MSM using TICA + clustering + MLE.
  [Post-proc; traj analyses]:
    - cluster_analysis: Cluster trajectory (KNN) and output centroid structures.
    - free_energy_landscape: Compute 2D free energy landscape (PCA).
    - interaction_fingerprint: Compute and plot protein-ligand IFPs (Prolif).
    - protein_protein_interaction_fingerprint: Compute and plot protein-protein IFPs (Prolif).
    - rmsd_analysis: Compute and plot RMSDs by component (protein, ligand, antibody, etc.)
    - rmsf_analysis: Compute per-residue RMSF for protein (all atoms and Cα).
Available Workflows:
 - molecular_dynamics: Perform molecular dynamics simulation using specified backend.
 - rbfe: Prepare and build OpenFE RBFE calculations.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## Notes

* General MD
    * At the moment, the MD code outputs restarts every 10 ns, which is set in the yaml with `output_split_ns`. To restart a failed simulation, you just need to run the pipeline again with your yaml and it will automatically detect from the yaml parameter `length_ns` where to restart the simulation. 
    * If you want to run a longer simulation from a completed job, increase `length_ns`, run the pipeline again and it will automatically detect the appropriate restart point.
* Free energy 
    * Use the OpenFE env for running simulations they provide during installation. The code base is being revised fairly regularly so dependencies will change with newer versions.
    * The ligands in the sdf file you provide ``must be aligned``.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* More options for MD simulations.
    * Classic and ML force fields.
    * Loop modelling.
* More advanced simulation types.
    * QMMM.
    * Steered MD.
    * Metadynamics.