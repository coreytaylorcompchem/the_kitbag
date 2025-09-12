# Data curation pipeline

Here's where we can curate and combine affinity, ADME, etc, data from multiple sources. 

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
 - retrieve_chembl_adme_data: Retrieve CHEMBL ADME data.
 - retrieve_chembl_bioactivities: Retrieve bioactivity data from CHEMBL.
 - clean_bioactivities: Check and standardise bioactivities.
 - retrieve_compound_smiles: Retrieve SMILES from downloaded compound data.
 - annotate_bioactivity_pactivity: Compute p(readout)) and add to retrieval results.
 - clean_adme_data: Check and standardise ADME data.
 - merge_bioactivity: Merge data from different sources (WIP).

Available Workflows:
 - chembl_adme_data: Retrieve, standardise and collate ADME data - ChEMBL.
 - chembl_multi_target: Retrieve, standardise and collate bioactivities for multiple targets - CHEMBL.
 - chembl_bioactivity_single_target: Retrieve, standardise and collate bioactivities for a single target - CHEMBL.
 - chembl_tox_targets: Retrieve, standardise and collate bioactivities for tox-relevant targets - CHEMBL.

```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

* Implement more advanced curation (at the moment, ADME retrieves all readouts / units).
* Add more parsers for alternative data sources.
* Workflows to combine data from different sources.