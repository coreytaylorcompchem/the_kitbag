import os
import importlib
import sys

def import_modules_from_dir(directory: str, package: str):
    """
    Dynamically import all modules from a given directory/package.

    Args:
        directory (str): Absolute path to the directory.
        package (str): Dotted import path (e.g., 'workflows').
    """
    for filename in os.listdir(directory):
        if filename.endswith(".py") and not filename.startswith("__"):
            module_name = filename[:-3]
            full_module_name = f"{package}.{module_name}"
            importlib.import_module(full_module_name)
