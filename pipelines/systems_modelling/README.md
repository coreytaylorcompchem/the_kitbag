# Systems modelling pipeline

PK modelling for brain penetration, liver PK and others.

![Systems modelling logo](_images/systems_modelling_logo.jpg)

## Setup and dependencies

* [Rdkit](https://www.rdkit.org/docs/Install.html)
* [Scipy](https://scipy.org/install/) - solve ODEs.
* [Pytorch](https://pytorch.org/get-started/locally/) - instantiate Torch-based ML models.

See `environment.yml` for the rest. The list is indicative only, as this repo is under development.

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python ../../run_systems_modelling_pipeline.py --params ../../configs/pbpk_modelling.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `output/loop_modelling` directory should appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
Available Backends:
Available Tasks:
  [PBPK]:
    - pbpk_analysis: Analyse PBPK simulation results.
    - pbpk_model_assembly: Assemble PBPK model
    - pbpk_parameter_prediction: Predict PBPK parameters from SMILES
  [PBPK - brain:plasma]:
    - pbpk_simulation: Run 2-compartment plasma–brain PBPK simulation with selectable exchange model.

Available Workflows:
 - pbpk_workflow: Predict PBPK parameters, assemble model, simulate, and analyse.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* Model membrane proteins + membrane.
* Model proteins with the *Fold ecosystem once I have access to a decent GPU or two.