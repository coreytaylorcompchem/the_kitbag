import importlib.util
import sys
import logging

# Configure a logger or use your existing logger setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def check_dependencies(required_modules):
    """
    Checks if the required modules are installed.

    Args:
        required_modules (dict): Mapping of module names to install instructions

    Returns:
        missing (list): List of missing module names
    """
    missing = []

    for module, install_instr in required_modules.items():
        if importlib.util.find_spec(module) is None:
            logger.warning(f"Missing dependency: '{module}'. Install via: {install_instr}")
            missing.append((module, install_instr))
    
    return missing


def fail_if_missing(missing_modules):
    if missing_modules:
        print("\n🚫 Missing critical dependencies:")
        for mod, instr in missing_modules:
            print(f"  ❌ {mod} — Install via: {instr}")
        print("\n💡 Please install the missing dependencies before running the pipeline.")
        sys.exit(1)
