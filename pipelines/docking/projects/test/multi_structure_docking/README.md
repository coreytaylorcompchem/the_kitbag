# Test runs of pipeline - multi-structure docking

Here we see examples where we dock ligands exhastively into multiple PDB structures. 
The outputs show results from all docking runs as well as generate a couple of summary csvs and plots. 

## Process

`python ../../../run_docking_pipeline.py --params ../../../configs/multi_structure_docking_gnina.yaml`
`python ../../../run_docking_pipeline.py --params ../../../configs/multi_structure_docking_unidock.yaml`