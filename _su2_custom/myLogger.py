# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 17:26:52 2025

@author: BriceM

Trying my hand at building my own logger that can replace print statements that
will log and print to the console. This will currently only work for my own scripts

Update: I severely udnerstimated how well I understood logging
"""
# fdlogger_tee.py
import os, sys, time, threading
from pathlib import Path
from types import SimpleNamespace

class FDTeeLogger:
    """
    FD-level capture that *mirrors* output to both console and a log file.
    Works across Python/C extensions and most subprocesses.
    
    USAGE: 
        logger = FDLogger(filepath="logs", filename="mylogfile.log")
        function_outputs = logger.log(FunctionCall, 2, 3)
    """
    def __init__(self, filepath=".", filename="run.log", encoding="utf-8"):
        """
        
        Parameters
        ----------
        filepath : str, optional
            filepath to log file. The default is ".".
        filename : str, optional
            name of log file. The default is "run.log".
        encoding : str, optional
            The default is "utf-8".

        Returns
        -------
        None.

        """
        self.log_dir = Path(filepath); self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_path = self.log_dir / filename
        self.encoding = encoding
        self._saved = SimpleNamespace(out=None, err=None)
        self._threads = []
        self._pipes = []  # [(read_fd, write_fd, orig_fd, name), ...]
        self._stop_evt = threading.Event()

    def _spawn_tee(self, name, target_fd, log_file_handle):
        # Create a pipe; dup write end onto target_fd (1 or 2)
        rfd, wfd = os.pipe()
        orig = os.dup(target_fd)
        os.dup2(wfd, target_fd)
        os.close(wfd)

        def _pump():
            # Read from pipe, write both to original fd and log file
            with os.fdopen(rfd, "rb", buffering=0) as r:
                while not self._stop_evt.is_set():
                    chunk = r.read(8192)
                    if not chunk:
                        break
                    # Write to original console fd
                    os.write(orig, chunk)
                    # And to log
                    log_file_handle.buffer.write(chunk)  # raw bytes
                    log_file_handle.flush()

            os.close(orig)

        t = threading.Thread(target=_pump, name=f"FDTee-{name}", daemon=True)
        t.start()
        self._threads.append(t)
        self._pipes.append((rfd, orig, target_fd, name))

    def _start(self):
        self._stop_evt.clear()
        self._log_fh = open(self.log_path, "ab", buffering=0)
        # tee stdout (1) and stderr (2)
        self._spawn_tee("stdout", 1, self._log_fh)
        self._spawn_tee("stderr", 2, self._log_fh)

    def _stop(self):
        # Close the write ends by restoring original fds; then join threads
        # (Simplest way: dup the saved originals back—threads will see EOF)
        for rfd, orig_fd, target_fd, _ in self._pipes:
            try:
                os.dup2(orig_fd, target_fd)
            except OSError:
                pass
        self._stop_evt.set()
        for t in self._threads:
            t.join(timeout=1.0)
        self._threads.clear()
        self._pipes.clear()
        try:
            self._log_fh.close()
        except Exception:
            pass

    def log(self, func, *args, **kwargs):
        start = time.time()
        with open(self.log_path, "a", encoding=self.encoding) as ftxt:
            ftxt.write(f"\n===== START {func.__name__} {time.ctime(start)} =====\n")
            ftxt.flush()
        self._start()
        try:
            return func(*args, **kwargs)
        finally:
            self._stop()
            end = time.time()
            with open(self.log_path, "a", encoding=self.encoding) as ftxt:
                ftxt.write(f"===== END {func.__name__} {time.ctime(end)} (elapsed: {end-start:.2f}s) =====\n")
                ftxt.flush()

            
    