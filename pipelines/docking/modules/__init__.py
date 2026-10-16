import importlib
import pkgutil

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

from pipeline.optional_dependencies import OPTIONAL_DEPENDENCIES

def load_all_tasks():

    for _, name, _ in pkgutil.iter_modules(__path__):

        if name.startswith("_"):
            continue

        module_name = f"{__name__}.{name}"

        try:
            importlib.import_module(module_name)

        except ModuleNotFoundError as exc:

            missing_package = exc.name.split(".")[0]

            if missing_package in OPTIONAL_DEPENDENCIES:

                logger.debug(
                    f"Skipping task module '{module_name}': "
                    f"optional dependency '{exc.name}' is unavailable."
                )

                continue

            # This wasn't an optional dependency.
            # It could be a programming/import error.
            raise


load_all_tasks()
