import importlib
import pkgutil
from importlib.util import find_spec

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


OPTIONAL_MODULE_DEPENDENCIES = {
    "rbfe": [
        "openfe",
        "cinnabar",
        "openff",
    ],
}


def _missing_dependencies(module_name):
    deps = OPTIONAL_MODULE_DEPENDENCIES.get(module_name, [])
    missing = []

    for dep in deps:
        if find_spec(dep) is None:
            missing.append(dep)

    return missing


def load_all_tasks(
    selected_modules=None,
    include_optional=False,
):
    """
    Import task modules so their @register_task decorators run.

    Parameters
    ----------
    selected_modules : set[str] | list[str] | None
        If provided, only these module names are imported.
        Example: {"md"} or {"rbfe"}.

    include_optional : bool
        If False, modules with missing optional dependencies are skipped.
    """

    if selected_modules is not None:
        selected_modules = set(selected_modules)

    package_name = __name__

    for _, name, is_pkg in pkgutil.iter_modules(__path__):
        if is_pkg:
            continue

        if name.startswith("_"):
            continue

        if selected_modules is not None and name not in selected_modules:
            continue

        missing = _missing_dependencies(name)

        if missing and not include_optional:
            logger.info(
                f"Skipping optional task module modules.{name}; "
                f"missing dependencies: {', '.join(missing)}"
            )
            continue

        full_module_name = f"{package_name}.{name}"

        logger.debug(f"Importing task module: {full_module_name}")
        importlib.import_module(full_module_name)