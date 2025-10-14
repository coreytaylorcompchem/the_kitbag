# Structure modelling pipeline

Modelling of proteins and peptides, mainly.

![Structure modelling logo](_images/structure_modelling_logo.jpg)

## Setup and dependencies

* [Rdkit](https://www.rdkit.org/docs/Install.html)
* [Biopython](https://biopython.org/wiki/Download)
* [MODELLER](https://salilab.org/modeller/download_installation.html) (default backend)
* [OpenMM](https://docs.openmm.org/latest/userguide/application/01_getting_started.html)
* [PyRosetta](https://www.pyrosetta.org/downloads) 
  * NOTE: 1.7 Gb download
* [ColabFold](https://github.com/sokrypton/ColabFold)

See `environment.yml` for the rest. The list is indicative only, as this repo is under development.

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python ../../run_structure_modelling_pipeline.py --params ../../configs/loop_modelling_and_fixes.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `output/loop_modelling` directory should appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
Available Tasks:
  [Peptide modelling]:
    - build_peptide: Build peptides from sequence file
    - minimise_peptide: Minimize peptide structures
  [Protein modelling]:
    - fix_loops: Model missing loops.
    - fix_residues: Fix incomplete residues.
    - refine_loops: Refine loop structures.

Available Workflows:
 - peptide_modelling: Build and minimise peptides.
 - protein_modelling: Fix, model, and refine protein structure.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* Model membrane proteins + membrane.
* Model proteins with the *Fold ecosystem once I have access to a decent GPU or two.