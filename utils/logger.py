"""
Logging Utility
"""

import logging
import sys
from pathlib import Path


def setup_logger(
    name: str = "solvent_analysis",
    level: int = logging.INFO,
    log_file: str | None = None,
    console: bool = True,
) -> logging.Logger:
    """
    Configures and returns a logger.

    Args:
        name: Name of the logger.
        level: Logging level (e.g., logging.INFO, logging.DEBUG).
        log_file: Optional path to a log file. If provided, logs will be written here.
        console: Whether to output logs to the console.

    Returns:
        logging.Logger: A configured logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Remove any existing handlers to avoid duplicate logs
    logger.handlers.clear()

    # Create a formatter for log messages
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler (stdout)
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # File handler
    if log_file:
        # Ensure the log directory exists
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_logger(name: str = "solvent_analysis") -> logging.Logger:
    """
    Retrieves a logger. If it hasn't been configured, sets it up with defaults.

    Args:
        name: Name of the logger.

    Returns:
        logging.Logger: The requested logger instance.
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        # Apply default configuration if no handlers are present
        setup_logger(name)
    return logger


class LogMixin:
    """A mixin class that provides a convenient logger property to its subclasses."""

    @property
    def logger(self) -> logging.Logger:
        """Returns a logger instance specific to the class using this mixin."""
        if not hasattr(self, "_logger"):
            # Create a logger named after the class for better context
            class_name = self.__class__.__name__
            self._logger = get_logger(f"solvent_analysis.{class_name}")
        return self._logger