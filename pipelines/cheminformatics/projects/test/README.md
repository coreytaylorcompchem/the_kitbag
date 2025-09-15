# Test runs of pipeline

Here are test cases for the following use cases:

## Bioactivities

* Retrieve and clean bioactivities for a single target (EP4)

`python ../../run_data_pipeline.py --params ../../configs/chembl_pipeline_single_target.yaml` 

* Retrieve and clean bioactivities for a multiple targets (EP1 - EP4)

`python ../../run_data_pipeline.py --params ../../configs/chembl_pipeline_multi_targets.yaml`

## ADME/Tox

* Retrieve and clean values for LogD, microsomal stability (CLint), Solubility, Caco2 (A->B Papp) and PPB.

`python ../../run_data_pipeline.py --params ../../configs/chembl_pipeline_adme.yaml` 

* Retrieve and clean bioactivities for a tox targets (Eurofins44)

`python ../../run_data_pipeline.py --params ../../configs/chembl_pipeline_tox_targets.yaml`