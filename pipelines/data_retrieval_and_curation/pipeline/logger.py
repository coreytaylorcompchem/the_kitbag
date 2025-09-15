import logging
import sys
from pathlib import Path

def setup_logger(
    name: str = None,
    log_level: int = logging.INFO,
    log_to_file: bool = False,
    log_file: str = "app.log",
    debug_mode: bool = False,
    simple_format: bool = False,
    minimalist_format: bool = False
) -> logging.Logger:
    """
    Set up and return a custom logger.
    """
    logger = logging.getLogger(name)
    if logger.hasHandlers():
        return logger

    level = logging.DEBUG if debug_mode else log_level
    logger.setLevel(level)

    if minimalist_format:
        formatter = logging.Formatter('%(levelname)s %(message)s')
    elif simple_format:
        formatter = logging.Formatter('%(levelname)s - %(name)s - %(message)s')
    else:
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(name)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Optional file handler
    if log_to_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger
