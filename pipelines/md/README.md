# Molecular dynamics pipeline

Perform different types of MD in automated, end-to-end fashion. 

## Setup and dependencies

* Mainly Python libs
* [Cuda](https://developer.nvidia.com/cuda-toolkit)
* [Rdkit](https://www.rdkit.org/)
* [OpenMM](https://openmm.org/)
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
Importing workflow module: workflows.md_workflow
INFO 
Available Tasks:
INFO 
Available Workflows:
INFO  - molecular_dynamics: Perform MD.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

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