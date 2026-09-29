import importlib
import pkgutil

import modules


def load_all_tasks():
    """
    Import task modules explicitly during pipeline startup.

    Do not run this automatically when the modules package
    is imported, because utility imports such as
    modules.utils.transforms must remain side-effect free.
    """
    for _, name, _ in pkgutil.iter_modules(
        modules.__path__
    ):
        if name.startswith("_"):
            continue

        importlib.import_module(
            f"{modules.__name__}.{name}"
        )