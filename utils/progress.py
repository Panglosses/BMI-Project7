"""
进度条工具
"""

import sys
import time
from typing import Optional


class ProgressBar:
    """简单进度条"""

    def __init__(
        self,
        total: int,
        prefix: str = "",
        suffix: str = "",
        length: int = 50,
        fill: str = "█",
        print_end: str = "\r",
    ):
        """
        Args:
            total: 总任务数
            prefix: 前缀文本
            suffix: 后缀文本
            length: 进度条长度
            fill: 填充字符
            print_end: 打印结束字符
        """
        self.total = total
        self.prefix = prefix
        self.suffix = suffix
        self.length = length
        self.fill = fill
        self.print_end = print_end
        self.start_time = time.time()
        self.current = 0

    def update(self, iteration: int):
        """更新进度"""
        self.current = iteration
        percent = f"{100 * (iteration / float(self.total)):.1f}"
        filled_length = int(self.length * iteration // self.total)
        bar = self.fill * filled_length + "-" * (self.length - filled_length)

        # 计算剩余时间
        elapsed = time.time() - self.start_time
        if iteration > 0:
            eta = (elapsed / iteration) * (self.total - iteration)
            eta_str = f"ETA: {eta:.1f}s"
        else:
            eta_str = "ETA: --"

        sys.stdout.write(
            f"\r{self.prefix} |{bar}| {percent}% "
            f"({iteration}/{self.total}) {eta_str} {self.suffix}"
        )
        sys.stdout.flush()

    def increment(self):
        """增加一个进度"""
        self.update(self.current + 1)

    def finish(self):
        """完成进度条"""
        self.update(self.total)
        print()  # 换行

    @staticmethod
    def iterate(iterable, total: Optional[int] = None, **kwargs):
        """
        迭代器包装器

        Args:
            iterable: 可迭代对象
            total: 总数量（如果可迭代对象没有len()）
            **kwargs: 传递给ProgressBar的参数

        Yields:
            迭代元素
        """
        if total is None:
            try:
                total = len(iterable)
            except TypeError:
                raise ValueError("必须提供total参数或可迭代对象支持len()")

        progress = ProgressBar(total, **kwargs)
        for i, item in enumerate(iterable):
            yield item
            progress.update(i + 1)
        progress.finish()