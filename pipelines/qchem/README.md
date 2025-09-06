# Quantum chemistry pipeline

Here's where we can run some quantum chemistry calculations. 

## Setup and dependencies

You'll need `psi4` and several others. Will create an env.yaml in time. Further explanations on what the other bits do - soon.

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

Just to verify all the code and software is running and installed correctly, there's a test run in `projects/test` Run the code from there with:

`python ../../main.py water.xyz ../../configs/config_default_workflow.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `results` directory appear.

## So what can it do?

Not much - yet.

The frame workwork to build workflows is there, so it should be a matter of implmenting the specific calculations so workflows can actually run them. The framework is currently implmented in such a way that it will eventually be possible to stack calculations with checkpointing at each step so if there's a psi4 crash, it won't need to re-run steps. So you could run a quick opt with XTB, run some other structural steps, another optimisation with a DFT method then more advanced calculations on the DFT geometry (e.g. sapt0).

The code is implemented with automatic checkers to see what's actually available. Run `python check_available_tasks_workflows_backends.py` for a list.

```
python check_available_tasks_workflows_backends.py 
Available Backends:
 - psi4
 - xtb

Available Tasks:
 - mesp_map: Calculates a Molecular Electrostatic Potential Map
   ↳ Backends: None
 - optimise: Performs geometry optimization using the selected backend.
   ↳ Backends: psi4, xtb
 - pipeline: Allows us to stack calculations in a workflow
   ↳ Backends: None
 - sapt0: Do a sapt0 calculation (psi4 only)
   ↳ Backends: None
 - torsion_scan: Do a torsion scan of a bond.
   ↳ Backends: None

Available Workflows:
 - advanced: optimize → mesp_map → torsion_scan → sapt0
 - advanced_multistage_workflow: Runs opt → torsion → opt → MESP with different backends.
 - default: Performs structure optimisation only.
```

## TODO

* Implement XTB and Orca API to give further options.
* Add all the common calculations (torsion scans, mesp maps, etc.)
