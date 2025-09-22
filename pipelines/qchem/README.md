# Quantum chemistry pipeline

Here's where we can run some quantum chemistry calculations. 

![Quantum chemistry logo](_images/qchem.png)

## Setup and dependencies

You'll need `psi4` and several others. Will create an env.yaml in time. Further explanations on what the other bits do - soon.

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

Just to verify all the code and software is running and installed correctly, there's a test run in `projects/test` Run the code from there with:

`python ../../main.py water.xyz ../../configs/config_default_workflow.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `results` directory appear.

## So what can it do?

Not much - yet.

The framework to build workflows is there, so it should be a matter of implementing the specific calculations so workflows can actually run them. The framework is currently implemented in such a way that it will eventually be possible to stack calculations with checkpointing at each step so if there's a psi4 crash, it won't need to re-run steps. So you could run a quick opt with XTB, run some other structural steps, another optimisation with a DFT method then more advanced calculations on the DFT geometry (e.g. sapt0).

The code is implemented with automatic checkers to see what's actually available. Run `python check_available_tasks_workflows_backends.py` for a list.

```
Available Tasks:
  [EDA]:
    - qtaim: Do AIM (Atoms In Molecules) analysis.
    - sapt0: Do a sapt0 calculation (psi4 only)
  [Energy]:
    - single_point: Performs a single-point energy calculation using the selected backend.
  [Infra]:
    - pipeline: Allows us to stack calculations in a workflow
  [Molecular properties]:
    - mesp_map: Calculates a Molecular Electrostatic Potential Map
  [Optimization]:
    - optimise: Performs geometry optimization using the selected backend.
  [PES exploration]:
    - torsion_scan: Do a torsion scan of a bond.

Available Workflows:
 - staged_workflow: Performs single or multi-stage calculations.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

* Implement XTB and Orca API to give further options.
* Add all the common calculations (torsion scans, mesp maps, etc.)
