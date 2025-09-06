import argparse
from smd import SMDPipeline
from qmmm import QMMMPipeline
from metadynamics import MetaDynamicsPipeline
from md import MDPipeline

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--mode", choices=["default", "smd", "qmmm", "metad"], default="default")
    args = parser.parse_args()

    if args.mode == "smd":
        workflow = SMDPipeline.from_yaml(args.config)
    elif args.mode == "qmmm":
        workflow = QMMMPipeline.from_yaml(args.config)
    elif args.mode == "metadynamics":
        workflow = MetaDynamicsPipeline.from_yaml(args.config)
    else:
        workflow = MDPipeline.from_yaml(args.config)

    workflow.run()

if __name__ == "__main__":
    main()
