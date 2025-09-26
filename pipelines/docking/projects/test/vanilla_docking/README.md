# Test runs of pipeline - vanilla docking

Here we see examples where we dock ligands into a single protein with the binding pocket determined from the crystal ligand coords. 
The outputs show results from all docking runs (SDFs or pdbqt, depending on the backend used). 

## Process

`python ../../../run_docking_pipeline.py --params ../../../configs/vanilla_docking_gnina.yaml`
`python ../../../run_docking_pipeline.py --params ../../../configs/vanilla_docking_unidock.yaml`