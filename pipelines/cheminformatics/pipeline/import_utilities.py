import importlib
import pkgutil
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def import_modules_recursively(base_dir: str | Path, base_package: str):
    base_dir = Path(base_dir)
    if not base_dir.exists():
        logger.error(f"Base directory does not exist: {base_dir}")
        return

    def _recursive_import(package_path: Path, package_name: str):
        for finder, name, is_pkg in pkgutil.iter_modules([str(package_path)]):
            full_name = f"{package_name}.{name}"
            try:
                importlib.import_module(full_name)
                logger.info(f"Imported module {full_name}")
            except Exception as e:
                logger.error(f"Failed to import module '{full_name}': {e}")

            if is_pkg:
                _recursive_import(package_path / name, full_name)

    _recursive_import(base_dir, base_package)
