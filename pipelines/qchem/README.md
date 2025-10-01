# Quantum chemistry pipeline

Here's where we can run some quantum chemistry calculations. 

![Quantum chemistry logo](_images/qchem.png)

## Setup and dependencies

* Python > 3.11 (OPI needs it)
* [Rdkit](https://www.rdkit.org/docs/Install.html)
* [psi4](https://psicode.org/installs/v110zero/)
* [Orca](https://orcaforum.kofo.mpg.de/app.php/portal) - 6.1 or above (big download ~500Mb from their forum)
* [Orca Python Interface (OPI)](https://github.com/faccts/opi)
* [Xtb](https://xtb-docs.readthedocs.io/en/latest/setup.html) - binary is well supported, python-xtb not
* [Multiwfn](http://sobereva.com/multiwfn/) - this software is great and the dude who maintains it is worth supporting.
* [libxm4](https://packages.debian.org/search?keywords=libxm4) - for Multiwfn

Will create an env.yaml in time. Further explanations on what the other bits do - soon.

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

Just to verify all the code and software is running and installed correctly, there's a test run in `projects/test` Run the code from there with:

`python ../../run_qchem_pipeline --params ../../configs/pipeline_single_stage.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `results` directory appear.

## So what can it do?

The code is implemented with automatic checkers to see what's actually available. Run `python list_available_tasks_and workflows.py` for a list.

```
Available Backends:
 - multiwfn
 - orca
 - psi4
 - xtb

Available Tasks:
    [EDA]:
    - sapt0: Do a sapt0 calculation (psi4 only)
  [Energy]:
    - single_point: Performs a single-point energy calculation using the selected backend.
    - torsion_scan: Do a torsion scan of a bond.
  [Optimization]:
    - optimise: Performs geometry optimization.
  [Property]:
    - calc_pka: Calculate site pKa values from HA/A- free energies.
    - mesp: Calculates Molecular Electrostatic Potential (MESP) and outputs cube files for visualization.
    - site_pka: Enumerate acidic sites, generate protonated and deprotonated geometries for pKa calculations.
  [Wave function analysis]:
    - basin: Performs real-space analysis (basins, DIs).
    - nci: Performs NCI calculation, outputs cube files.
    - qtaim: Performs QTAIM analysis on .molden formatted file.

Available Workflows:
 - staged_workflow: Performs single or multi-stage calculations.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* Implement XTB and Orca API to give further options.
* For pKa, to make it actually useful;
  * Add solvent to calculation (Orca only, Psi4 doesn't handle solvent well).
  * Add empirical correction to solvent pKas.
    * Will use solvent pKa calculation to generate dataset from literature values.
  * Generate multiple energetic conformers to feed into pKa calculation to get the average pKa.
    * Common tactic; pKa is geometry dependent. 
* Add all the common calculations (torsion scans, mesp maps, etc.)
