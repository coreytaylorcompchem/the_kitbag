# Code assessment pipeline

![Code assessment logo](_images/code_assessment_logo.png)

In this repo is code to assess the performance of the other pipelines (WIP). 

## Setup and dependencies

* aaa

See `environment.yml` for the rest. The list is indicative only, as this repo is under development.

## Running calculations

I use a yaml/workflow system. Examples for each are in `configs/*yaml` and `workflows/*py`.

See a few workflow runs in `projects/test` Run the code from there with, for example:

`python run_data_pipeline --params ../../configs/vanilla_docking.yaml > log.log`

If the calculation runs correctly, you should see output (`log.log`) and a `output/` directory should appear.

## So what can it do?

`python list_available_tasks_and_workflows.py`

```
 Available Tasks:
Available Workflows:
 
```

If you want to register new workflows, you'll also need to do so with metadata to describe what it does.

## TODO

See [issues](https://github.com/coreytaylorcompchem/the_kitbag/issues) for a more complete list.

* aaa