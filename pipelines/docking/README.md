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
 - standardize_ligand: Prep: standardise ligand from SMILES.
 - generate_conformers: Prep: generate RDKit conformers.
 - cluster_conformers: Prep: cluster and select conformers.
 - optimize_with_xtb: Prep: optimise conformers using GFN1-xTB.
 - save_final_conformers: Prep: save final conformers to sdf.
 - convert_to_pdbqt: Prep: convert final conformers to PDBQT for docking.
 - dock: Docking: run docking.

Available Workflows:
 - constrained_docking: Prepare, dock and score with core constraints.
 - ensemble_docking: Prepare, dock and score emsemble of proteins.
 - vanilla_docking: Prepare, dock and score with no constraints.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

* Add more docking backends ([Uni-dock](https://github.com/dptech-corp/Uni-Dock) and ML-based backends like Diffdock and Equibind).
* Add other docking modes (constrained core, ensemble, etc.)
* Run docking on multiple targets. 