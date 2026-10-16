import importlib
import pkgutil

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

from pipeline.optional_dependencies import OPTIONAL_DEPENDENCIES

# Registry for workflows
_workflow_registry = {}

def register_workflow(name, description=None):
    def decorator(func):
        _workflow_registry[name.lower()] = {
            "func": func,
            "description": description
        }
        return func
    return decorator

def get_workflow(name):
    entry = _workflow_registry.get(name.lower())
    return entry["func"] if entry else None

def list_workflows():
    return list(_workflow_registry.keys())

def get_workflow_metadata(name):
    return _workflow_registry.get(name.lower(), {})

def load_all_workflows():

    for _, name, _ in pkgutil.iter_modules(__path__):

        if name.startswith("_"):
            continue

        module_name = f"{__name__}.{name}"

        logger.debug(
            f"Discovered workflow module: {module_name}"
        )

        try:
            importlib.import_module(module_name)

        except ModuleNotFoundError as exc:

            missing_package = exc.name.split(".")[0]

            if missing_package in OPTIONAL_DEPENDENCIES:

                logger.warning(
                    f"Skipping workflow module '{module_name}': "
                    f"optional dependency '{exc.name}' is unavailable."
                )

                continue

            # Genuine import/programming error
            logger.error(
                f"Failed to import workflow module '{module_name}': "
                f"missing module '{exc.name}'."
            )

            raise
# Load all workflows at import time
load_all_workflows()