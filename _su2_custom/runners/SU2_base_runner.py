# -*- coding: utf-8 -*-
from __future__ import annotations

import os
import time
import subprocess
from datetime import timedelta


def _format_duration(seconds: float) -> str:
    """Format seconds as H:MM:SS (timedelta handles >24h)."""
    return str(timedelta(seconds=int(seconds)))



def base_CFD_run(cfg_filename: str, cfg_filepath: str, save_output: bool = False) -> None:
    '''
    Runs the configuration file cfg_filename located within cfg_filepath through SU2.

    Parameters
    ----------
    cfg_filename : str
        Full filename of the configuration file: Ex. "filename.cfg". 
    cfg_filepath : str
        Full filepath to the configuration file (Ends with folder containing .cfg file).
    save_output : bool, optional
        Enables the log file writing when True, saves the log to "cfg_filepath". The default is False.

    Returns
    -------
    None
        

    '''
    
    run_cmd = f"%SU2_RUN%\\SU2_CFD {cfg_filename}"
    start = time.perf_counter()

    if save_output:
        log_file = (
            os.path.join(cfg_filepath, cfg_filename[:-4] + ".log")
            if cfg_filepath
            else cfg_filename[:-4] + ".log"
        )
        # Open the log file in write mode; append total time at the end.
        with open(log_file, "w", encoding="utf-8", errors="replace") as f:
            process = subprocess.Popen(
                run_cmd,
                shell=True,  # %SU2_RUN% expansion on Windows
                cwd=cfg_filepath or None,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,  # line-buffered
            )

            assert process.stdout is not None  # For type checkers
            for line in process.stdout:
                elapsed = _format_duration(time.perf_counter() - start)
                print(f"[{elapsed}] {line}", end="", flush=True)
                f.write(line)

            exit_code = process.wait()
            total = _format_duration(time.perf_counter() - start)
            f.write(f"\nTotal elapsed: {total}\n")

        print(f"\nProcess finished with exit code {exit_code} (total {total})")
        print(f"Output was saved to {log_file}")

    else:
        process = subprocess.Popen(
            run_cmd,
            shell=True,  # %SU2_RUN% expansion on Windows
            cwd=cfg_filepath or None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        assert process.stdout is not None
        for line in process.stdout:
            elapsed = _format_duration(time.perf_counter() - start)
            print(f"[{elapsed}] {line}", end="", flush=True)

        exit_code = process.wait()
        total = _format_duration(time.perf_counter() - start)
        print(f"\nProcess finished with exit code {exit_code} (total {total})")


if __name__ == "__main__":
    # Example usage (adjust paths as needed)
    config_filename = "3D_Swirler_test_RANS.cfg"
    config_filepath = "C:\\Users\\BriceM\\Documents\\SU2 CFD Data\\3D_Tests" #os.getcwd()
    basic_CFD_run(config_filename, config_filepath, save_output=True)