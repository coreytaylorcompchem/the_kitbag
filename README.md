Here's where all my end-to-end comp chem software pipelines live. Generally knitting together a bunch of stuff using Python. The aim is to make these able to run as independent pipelines or interoperate. Anyone who has ideas for how to best do that, drop me a message.

The code is organised into sections that reflect various computational chem software for various computational chem pipelines;

* Docking
    * Modelling
    * Running the docking (Gnina, generally)
    * Post-processing
    * Benchmarking
* MD
    * System prep, setup, run.
    * Post-processing
    * Automated analysis
* Quantum chemistry
    * Setup, structure opt, SP
    * Post-processing
    * Benchmarking
* Cheminformatics
    * Not sure how to do this one yet

TODO

* Y'know, build the pipelines.
* GPU code will be used as much as possible, so you'll need CUDA.
* TDD will in-principle be used (this will always be a TODO won't it).
