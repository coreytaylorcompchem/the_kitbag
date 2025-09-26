# Test runs of pipeline

Here are test cases for the following use cases (see yamls for details on parameters used):

## Filtering

### Lipinski, physchem and toxic / reactive SMARTS

`python ../../run_cheminformatics_pipeline.py --params ../../configs/chembl_pipeline_tox_reactive.yaml`

### Lipinski, physchem and tox / reactive SMARTS + ML ADME predictors

NOTE: at present, even though the filtering in the yaml will be done, only the ADME results will be returned.

`python ../../run_cheminformatics_pipeline.py --params ../../configs/chembl_pipeline_tox_reactive_adme.yaml`

### Lipinski, physchem and tox / reactive SMARTS + scaffold novelty

`python ../../run_cheminformatics_pipeline.py --params ../../configs/chembl_pipeline_physchem_scaffold_novelty.yaml`

## Library generation

### Focused separate libraries based on detected scaffolds

`python ../../run_cheminformatics_pipeline.py --params ../../configs/chembl_pipeline_generate_focused_lib_by_scaffold.yaml`

### Combinatorial library from R-group enumeration.

`python ../../run_cheminformatics_pipeline.py --params ../../configs/chembl_pipeline_generate_focused_lib_combinatorial.yaml`

### Reaction-based enumeration (Hartenfeller reactions x Enamine fragment libraries)

`python ../../run_cheminformatics_pipeline.py --params ../../configs/chembl_pipeline_generate_rxn_based_enumeration.yaml`