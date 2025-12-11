"""
Progress Bar Utility
"""

import sys
import time


class ProgressBar:
    """A simple progress bar for tracking task completion."""

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
            total: Total number of iterations/tasks.
            prefix: Text to display before the progress bar.
            suffix: Text to display after the progress bar and ETA.
            length: Character length of the progress bar.
            fill: Character used to fill the completed portion.
            print_end: End character for prints (e.g., '\\r' for updating in place).
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
        """Updates the progress bar display for the given iteration."""
        self.current = iteration
        percent = f"{100 * (iteration / float(self.total)):.1f}"
        filled_length = int(self.length * iteration // self.total)
        bar = self.fill * filled_length + "-" * (self.length - filled_length)

        # Calculate estimated time remaining
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
        """Increments the progress by one step and updates the display."""
        self.update(self.current + 1)

    def finish(self):
        """Completes the progress bar, ensuring it shows 100%."""
        self.update(self.total)
        print()  # Move to the next line

    @staticmethod
    def iterate(iterable, total: int | None = None, **kwargs):
        """
        A wrapper that yields items from an iterable while displaying a progress bar.

        Args:
            iterable: The object to iterate over.
            total: Total number of items (required if iterable has no __len__).
            **kwargs: Additional arguments passed to the ProgressBar constructor.

        Yields:
            Each item from the iterable.
        """
        if total is None:
            try:
                total = len(iterable)
            except TypeError:
                raise ValueError("A 'total' argument must be provided if the iterable has no __len__.")

        progress = ProgressBar(total, **kwargs)
        for i, item in enumerate(iterable):
            yield item
            progress.update(i + 1)
        progress.finish()