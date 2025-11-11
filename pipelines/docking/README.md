# Docking pipeline

![Docking logo](_images/docking_logo.png)

In this repo are workflows to dock compounds with various use cases (unbiased, constrained core, etc.). 

## Setup and dependencies

* [Rdkit](https://www.rdkit.org/docs/Install.html)
* [gnina](https://github.com/gnina/)
* [Unidock](https://github.com/dptech-corp/Uni-Dock)
* [propka](https://propka.readthedocs.io/en/latest/)
* [xtb](https://xtb-docs.readthedocs.io/en/latest/setup.html)
* [fpocket](https://github.com/Discngine/fpocket)
* [LightGBM](https://lightgbm.readthedocs.io/en/latest/Installation-Guide.html)
* Most Python libs.

See `environment.yml` for the rest. The list is indicative only, as this repo is under development.

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python run_data_pipeline --params ../../configs/vanilla_docking.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `output/` directory should appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
 Available Tasks:
  [Docking]:
    - active_learn_docking: Active learning loop over docking.
       ↳ Backends: gnina
    - dock: Dock with backend specified.
       ↳ Backends: gnina, unidock
    - induced_fit_docking: Dock, minimise nearby residues and re-dock.
       ↳ Backends: gnina
  [Ligand preparation]:
    - cluster_conformers: Cluster conformers by specified energy and RMSD criteria.
       ↳ Backends: gnina, unidock
    - convert_to_pdbqt: Convert ligands to pdbqt for docking.
       ↳ Backends: gnina, unidock
    - generate_conformers: Generate multiple feasible conformers from ligands.
       ↳ Backends: gnina, unidock
    - save_final_conformers: Save conformers that meet energy and RMSD criteria.
       ↳ Backends: gnina, unidock
    - standardise_ligand: Prepare and generate 3D coords for each ligand.
       ↳ Backends: gnina, unidock
  [Receptor preparation]:
    - prepare_receptor_pdbqt: Prepare the receptor for docking. Format: pdbqt
       ↳ Backends: gnina, unidock
Available Workflows:
 - active_learning_docking: prepare, dock and score with Active Learning loop.
 - constrained_docking: Prepare, dock and score with core constraints.
 - induced_fit_docking: Perform induced fit docking.
 - multi_pocket_docking: Dock ligands into top N pockets detected by fpocket.
 - multi_structure_docking: Dock ligands into multiple PDB structures.
 - vanilla_docking: Prepare, dock and score with no constraints.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## Performance

[Click here for a comparison of Gnina vs Unidock](https://docs.google.com/spreadsheets/d/1b8IvSBlcC0kmzoTG9YV1Gm6nOiDCGrEKMdmMHDxkT7c/edit?gid=0#gid=0). The winner in this case is clearly **Gnina** but I'm not 100% sure if I have optimal parameters for Unidock.

All docking methods run in parallel. I've set the maximum number of CPUs up to 20 - 2 to be safe, but maybe I'll add some options to the yaml to tune this.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* Add ML-based docking backends like [Diffdock](https://github.com/gcorso/DiffDock) and [Equibind](https://github.com/HannesStark/EquiBind).