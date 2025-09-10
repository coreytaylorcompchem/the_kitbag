# Test runs of pipeline

Here are test cases for the following use cases:

* Retrieve and clean bioactivities for a single target (EP4)

`python ../../run_data_pipeline.py --params ../../configs/chembl_pipeline_single_target.yaml` 

* Retrieve and clean bioactivities for a multiple targets (EP1 - EP4)

`python ../../run_data_pipeline.py --params ../../configs/chembl_pipeline_multi_targets.yaml`

* Retrieve and clean bioactivities for a tox targets (Eurofins44)

`python ../../run_data_pipeline.py --params ../../configs/chembl_pipeline_tox_targets.yaml` 
