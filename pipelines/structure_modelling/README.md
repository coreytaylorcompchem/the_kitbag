# Structure modelling pipeline

Modelling of proteins and peptides, mainly.

![Structure modelling logo](_images/structure_modelling_logo.jpg)

## Setup and dependencies
### Structure modelling

* [Rdkit](https://www.rdkit.org/docs/Install.html)
* [Biopython](https://biopython.org/wiki/Download)
* [MODELLER](https://salilab.org/modeller/download_installation.html) (default backend)
* [OpenMM](https://docs.openmm.org/latest/userguide/application/01_getting_started.html)
* [PyRosetta](https://www.pyrosetta.org/downloads) 
  * NOTE: 1.7 Gb download

### Structure prediction

* [Boltz-2](https://github.com/jwohlwend/boltz)
* [ColabFold](https://github.com/sokrypton/ColabFold)
* [Chai-1](https://github.com/chaidiscovery/chai-lab)
* [OpenFold](https://github.com/aqlaboratory/openfold)

See `environment.yml` for the rest. The list is indicative only, as this repo is under development.

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python ../../run_structure_modelling_pipeline.py --params ../../configs/loop_modelling_and_fixes.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `output/loop_modelling` directory should appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
Available Backends:
 - base
 - boltz
 - chai
 - colabfold
 - modeller
 - openfold
 - peptidebuilder
 - pyrosetta
 - structureinference
 - xtb
Available Tasks:
  [Peptide modelling]:
    - analyse_peptide_flexibility: Compute peptide flexibility metrics.
    - build_peptide: Build peptides from a sequence file.
    - generate_peptide_conformers: Generate and minimise (MMFF) peptide conformers.
    - plot_peptide_flexibility: Plot flexibility (Boltzmann-weighted).
  [Protein modelling]:
    - fix_loops: Model missing loops.
    - fix_residues: Fix incomplete residues.
    - refine_loops: Refine modelled loops (energy criteria).
  [Structure modelling]:
    - predict_structures: Run structure prediction tools (Chai, Boltz, etc.) in parallel.
    - rank_predictions: Rank AI-predicted structures by confidence metrics.

Available Workflows:
 - peptide_modelling: Build, minimise peptides and calculate metrics.
 - protein_modelling: Fix, model, and refine protein structure.
 - structure_prediction: Parallel AI-based structure prediction

```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* Model membrane proteins + membrane.
* Model proteins with the *Fold ecosystem once I have access to a decent GPU or two.