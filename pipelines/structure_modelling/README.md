# Structural modelling pipeline

Modelling of proteins and peptide.

![Structure modelling logo](_images/structure_modelling_logo.jpg)

## Setup and dependencies

* Rdkit
* PyRosetta
* ColabFold
* PDBFixer
* 

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python ../../run_structure_modelling_pipeline.py --params ../../configs/aaa.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `output/physchem` directory should appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
Available Tasks:
  [aaa]:
    - aaa
Available Workflows:
 - aaa
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* Model membrane proteins + membrane.