"""
Periodic code efficiency assessment entry point.
Now supports command-line arguments for flexibility.
"""

import argparse
from modules.main import run_assessment


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run static and performance assessment on a codebase."
    )
    parser.add_argument(
        "--code-dir",
        type=str,
        default="src",
        help="Path to the directory containing the pipeline code to assess.",
    )
    parser.add_argument(
        "--workflows-dir",
        type=str,
        default="workflows",
        help="Path to directory containing YAML workflow definitions.",
    )
    parser.add_argument(
        "--profile-entrypoint",
        type=str,
        default="src/run_pipeline.py",
        help="Path to the pipeline entrypoint script to profile.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="reports",
        help="Directory where reports will be saved.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_assessment(
        code_dir=args.code_dir,
        workflows_dir=args.workflows_dir,
        profile_entrypoint=args.profile_entrypoint,
        output_dir=args.output_dir,
    )
