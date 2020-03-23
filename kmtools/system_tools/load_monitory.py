import os
import signal
import threading
import time

import psutil


class LoadMonitor(threading.Thread):
    """Monitor CPU utilization of child processes and kill those processes when CPU utilization
    falls below a certain threshold.
    """

    def __init__(
        self, pid: int, poll_interval: float = 60, min_utilization: float = 10, *args, **kwargs
    ) -> None:
        """Initialize a `LoadMonitor` instance.

        Args:
            pid: Process id of the parent process. CPU utilization will be monitored for
                the children of this process.
            poll_interval: Number of seconds to sleep between CPU utilization checks.
            min_utilization: The minimum average utilization (in percent) for the child processes.
            *args: Positional arguments to pass up to the parent class.
            **kwargs: Keyword arguments to pass up to the parent class.
        """
        super().__init__(*args, **kwargs)
        self.pid = pid
        self.poll_interval = poll_interval
        self.min_utilization = min_utilization
        self._stop_event = threading.Event()

    def run(self) -> None:
        proc = psutil.Process(os.getpid())
        prev_utilization = 100.0
        while not self.stopped():
            child_procs = list(proc.children(recursive=True))

            utilizations = []
            for child_proc in child_procs:
                try:
                    utilizations.append(child_proc.cpu_percent(interval=1))
                except psutil.NoSuchProcess:
                    continue

            if utilizations:
                cur_utilization = sum(utilizations) / len(utilizations)
                if (
                    cur_utilization < self.min_utilization
                    and prev_utilization < self.min_utilization
                ):
                    self._kill_all_processes(child_procs)
                    self.stop()
                prev_utilization = cur_utilization
                time.sleep(self.poll_interval)
            else:
                time.sleep(self.poll_interval)

    def _kill_all_processes(self, proc_list):
        for proc in proc_list:
            try:
                os.kill(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    def stop(self):
        """Stop the execution of the thread."""
        self._stop_event.set()

    def stopped(self) -> bool:
        """Returns `True` if the thread has been signalled to stop."""
        return self._stop_event.is_set()
