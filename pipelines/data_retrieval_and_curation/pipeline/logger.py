import logging
import sys
from pathlib import Path

# ANSI color codes
RESET = "\x1b[0m"
BLACK = "\x1b[30m"
RED = "\x1b[31m"
GREEN = "\x1b[92m"
YELLOW = "\x1b[33m"
BLUE = "\x1b[34m"
MAGENTA = "\x1b[35m"
CYAN = "\x1b[36m"
WHITE = "\x1b[37m"
BOLD = "\x1b[1m"

LEVEL_COLORS = {
    logging.DEBUG: BOLD + CYAN,
    logging.INFO: BOLD + GREEN,
    logging.WARNING: BOLD + YELLOW,
    logging.ERROR: RED,
    logging.CRITICAL: BOLD + RED
}

class ColorFormatter(logging.Formatter):
    def format(self, record):
        color = LEVEL_COLORS.get(record.levelno, RESET)
        record.msg = f"{color}{record.msg}{RESET}"
        return super().format(record)

def setup_logger(
    name: str = None,
    log_level: int = logging.INFO,
    log_to_file: bool = False,
    log_file: str = "app.log",
    debug_mode: bool = False,
    simple_format: bool = False,
    minimalist_format: bool = False,
    spartan_format: bool = False
) -> logging.Logger:
    """
    Set up and return a custom logger with optional color coding.
    """
    logger = logging.getLogger(name)
    if logger.hasHandlers():
        return logger

    level = logging.DEBUG if debug_mode else log_level
    logger.setLevel(level)

    if spartan_format:
        fmt_str = '%(message)s'
    elif minimalist_format:
        fmt_str = '%(levelname)s %(message)s'
    elif simple_format:
        fmt_str = '%(levelname)s - %(name)s - %(message)s'
    else:
        fmt_str = '%(asctime)s - %(levelname)s - %(name)s - %(message)s'
        datefmt = '%Y-%m-%d %H:%M:%S'

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    if spartan_format or minimalist_format or simple_format:
        ch.setFormatter(ColorFormatter(fmt_str))
    else:
        ch.setFormatter(ColorFormatter(fmt_str, datefmt=datefmt))
    logger.addHandler(ch)

    # Optional file handler (no color)
    if log_to_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        if spartan_format or minimalist_format or simple_format:
            fh.setFormatter(logging.Formatter(fmt_str))
        else:
            fh.setFormatter(logging.Formatter(fmt_str, datefmt=datefmt))
        logger.addHandler(fh)

    return logger
