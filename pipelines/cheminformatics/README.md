# Cheminformatics pipeline

Perform all sorts of cheminf calculations from basic physchem filtering to more advanced analyses. The overall goal is to squeeze down huge lists of molecules from data extraction to something more manageable and representative for virtual screening.

## Setup and dependencies

* Mainly Python libs
* Rdkit

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python ../../run_cheminformatics_pipeline.py --params ../../configs/chembl_pipeline_physchem.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `output/physchem` directory should appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
INFO 
Available Tasks:
INFO  - physchem_filtering: Physchem filtering with mandatory and optional cutoffs + conformer generation
INFO 
Available Workflows:
INFO  - physchem_filtering: Parallel physchem filtering using RDKit descriptors
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

* More advanced filtering:
    * Shape-based and pharmacophore filtering
    * Tox filters (e.g. PAINS)
    * Clustering/PCA
    * Fragment/scaffold filtering
* More advanced analyses
    * MMPs
    * SAR analyses
    * Scaffold analyses
    * Use of ML models for ADME prediction
    * Free-Wilson analysis