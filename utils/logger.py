"""
日志工具
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
    设置日志记录器

    Args:
        name: 日志记录器名称
        level: 日志级别
        log_file: 日志文件路径（可选）
        console: 是否输出到控制台

    Returns:
        logging.Logger: 配置好的日志记录器
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # 清除现有的处理器
    logger.handlers.clear()

    # 创建格式化器
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # 控制台处理器
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # 文件处理器
    if log_file:
        # 确保日志目录存在
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_logger(name: str = "solvent_analysis") -> logging.Logger:
    """
    获取日志记录器

    Args:
        name: 日志记录器名称

    Returns:
        logging.Logger: 日志记录器
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        # 如果没有处理器，使用默认设置
        setup_logger(name)
    return logger


class LogMixin:
    """日志混合类"""

    @property
    def logger(self) -> logging.Logger:
        """获取类日志记录器"""
        if not hasattr(self, "_logger"):
            class_name = self.__class__.__name__
            self._logger = get_logger(f"solvent_analysis.{class_name}")
        return self._logger
