# Retrieve all Class A GPCR bioactivities

Here we download all Class A GPCR data from Chembl, clean and process it for further use.

## Bioactivities

`cd gpcrs`

* Retrieve and clean bioactivities for all Class A GPCRs

`python ../../../run_data_pipeline.py --params ../../../configs/chembl_pipeline_protein_class_a_gpcr.yaml` 

* Retrieve and clean bioactivities for all Class B GPCRs

`python ../../../run_data_pipeline.py --params ../../../configs/chembl_pipeline_protein_class_b_gpcr.yaml`

