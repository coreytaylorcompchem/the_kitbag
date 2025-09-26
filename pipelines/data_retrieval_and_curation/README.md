# Data curation pipeline

Here's where we can curate and combine affinity, ADME, etc, data from multiple sources. 

![Data curation logo](_images/data_curation_logo.png)

## Setup and dependencies

Mainly Python libs and the [Chembl webresource_client](https://github.com/chembl/chembl_webresource_client).

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python run_data_pipeline --params ../../configs/chembl_workflow_bioactivity_single_target.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `output/single_target` directory should appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
Available Tasks:
  [ADME]:
    - clean_adme_data: Check and standardise ADME data.
    - retrieve_chembl_adme_data: Retrieve CHEMBL ADME data.
  [Bioactivity]:
    - annotate_bioactivity_pactivity: Compute p(readout)) and add to retrieval results.
    - clean_bioactivities: Check and standardise bioactivities.
    - retrieve_chembl_bioactivities: Retrieve bioactivity data from CHEMBL.
    - retrieve_compound_smiles: Retrieve SMILES from downloaded compound data.
    - retrieve_protein_class_target_list: Retrieve UniProt IDs for protein target class.
  [PDB]:
    - align_structures: Align all PDBs to the first retrieved structure.
    - retrieve_pdbs: Retrieve all PDBs for a UniProt ID.
    - standardise_pdbs: Standardise PDBs (Chain A + ligand, remove solvent).
  [Post-processing]:
    - merge_bioactivity: Merge data from different sources (WIP).

Available Workflows:
 - chembl_adme_data: Retrieve, standardise and collate ADME data - ChEMBL.
 - chembl_multi_target: Retrieve, standardise and collate bioactivities for multiple targets - CHEMBL.
 - chembl_bioactivity_single_target: Retrieve, standardise and collate bioactivities for a single target - CHEMBL.
 - chembl_tox_targets: Retrieve, standardise and collate bioactivities for tox-relevant targets - CHEMBL.
 - pdb_multi_target: Process PDB structures for multiple UniProt targets.
 - pdb_single_target: Process PDBs for a single UniProt target.
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* Implement more advanced curation (at the moment, ADME retrieves all readouts / units).
* Add more parsers for alternative data sources.
* Workflows to combine data from different sources.