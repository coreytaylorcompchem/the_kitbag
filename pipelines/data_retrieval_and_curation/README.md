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

* Download and clean (Chembl only):
  * Bioactivities for single targets
  * Bioactivities for multiple targets
  * Bioactivities for all tox targets on the Eurofins Safety44 panel.
  * ADME assay results for LogD, human microsomal/heps stability, solubility, Caco-2 permeability and hPPB.

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

* Implement more advanaced curation (e.g. Caco-2 current retrieves all readouts).
* Add more parsers for alternative data sources.
* Workflows to combine data from different sources.