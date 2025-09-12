# Docking pipeline

In this repo are workflows to dock compounds with various use cases (unbiased, constrained core, etc.). 

## Setup and dependencies

* gnina
* propka
* xtb
* Most Python libs.

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python run_data_pipeline --params ../../configs/vanilla_docking.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `output/` directory should appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
Available Backends:
 - gnina

Available Tasks:
 - standardize_ligand: Standardize ligand from SMILES.
 - generate_conformers: Generate RDKit conformers.
 - cluster_conformers: Cluster and select conformers.
 - optimize_with_xtb: Optimize conformers using GFN1-xTB.
 - save_final_conformers: Save final conformers to SDF.
 - convert_to_pdbqt: Convert final conformers to PDBQT.
 - dock: Run docking using Gnina backend.

Available Workflows:
 - constrained_docking: Dock with core constraints.
 - ensemble_docking: Ensemble docking with Gnina.
 - vanilla_docking: Preparation and docking using Gnina.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

* Add more docking backends ([Uni-dock](https://github.com/dptech-corp/Uni-Dock) and ML-based backends like Diffdock and Equibind).
* Add other docking modes (constrained core, ensemble, etc.)
* Run docking on multiple targets. 