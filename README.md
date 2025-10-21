![The Kitbag logo](_images/the_kitbag.jpg)

<!-- <img src="_images/the_kitbag.jpg" alt="The Kitbag logo" width="400" align="left" style="margin-right: 20px;" /> -->

# The Kitbag
## End-to-end and modular Python pipelines for computational chemistry

Here's where all my end-to-end comp chem software pipelines live. Generally knitting together a bunch of stuff using Python, yamls and workflows. The aim is to make these able to run as independent pipelines or interoperate. Anyone who has ideas for how to best do that, drop me a message.

## How does it work?

At the moment, the pipelines are separately stored and run in `pipelines`:

* [Data curation and retrieval](https://github.com/coreytaylorcompchem/the_kitbag/tree/main/pipelines/data_retrieval_and_curation)
* [Structure modelling](https://github.com/coreytaylorcompchem/the_kitbag/tree/main/pipelines/structure_modelling)
* [Systems modelling](https://github.com/coreytaylorcompchem/the_kitbag/tree/main/pipelines/systems_modelling)
* [Cheminformatics](https://github.com/coreytaylorcompchem/the_kitbag/tree/main/pipelines/cheminformatics)
* [Docking](https://github.com/coreytaylorcompchem/the_kitbag/tree/main/pipelines/docking)
* [Machine learning](https://github.com/coreytaylorcompchem/the_kitbag/tree/main/pipelines/machine_learning)
* [Molecular dynamics](https://github.com/coreytaylorcompchem/the_kitbag/tree/main/pipelines/md)
* [Quantum chemistry](https://github.com/coreytaylorcompchem/the_kitbag/tree/main/pipelines/qchem)

I haven't gotten around to putting together a `env.yml` for them, so you'll be chasing Python deps for now.

## What does / will it do?

The code is organised into sections that reflect various computational chem software for various computational chem pipelines;

* Data retrieval and curation
    * Retrieve data from public sources (Chembl, Pubchem, etc.)
    * Curate them for further calculations.
    * Smoosh them together for various purposes
* Structure modelling
    * Homology and loop modelling
    * Fixing structures
    * Modelling of proteins and short peptides
    * Minimisation
* Systems modelling
    * PBPK modelling
* Cheminformatics 
    * Various data-driven use cases for:
        * Hit discovery
        * Lead opt
        * Multi-target
* Docking
    * Modelling
    * Running the docking
    * Post-processing
    * Benchmarking
* Machine learning
    * Model training
        * Pose ranker
        * ADME models
        * Affinity models
    * AI tools
        * Molecular generation.
* MD
    * System prep, setup, run.
    * Post-processing
    * Automated trajectory analysis
* Quantum chemistry
    * Setup, structure opt, SPE, etc.
    * Post-processing
    * Benchmarking

## TODO

* Y'know, build the pipelines.
* GPU code will be used as much as possible, so you'll need CUDA.
* TDD will in-principle be used (this will always be a TODO won't it).
