# -*- coding: utf-8 -*-
import subprocess 


def basic_CFD_run(cfg_filename, cfg_filepath, save_output=False):
    run_cmd = f"%SU2_RUN%\\SU2_CFD {cfg_filename}"

    if save_output:
        log_file = cfg_filepath + '\\' + cfg_filename[:-4] + '.log'
    
        # Open the log file in write mode
        with open(log_file, "w", encoding="utf-8") as f:
            # Start the process
            process = subprocess.Popen(
                run_cmd,
                shell=True,                 # Needed for %SU2_RUN% expansion
                cwd=cfg_filepath,        # Set working directory to where the cfg file is
                stdout=subprocess.PIPE,     
                stderr=subprocess.STDOUT,   # Merge stderr into stdout
                text=True,                  # Decode bytes -> str automatically
                bufsize=1                   # Line-buffered
            )
        
            # Stream output and log simultaneously
            for line in process.stdout:
                print(line, end="")   # Show live in console
                f.write(line)         # Save to file
        
            # Wait for process to complete
            exit_code = process.wait()
        
        print(f"\nProcess finished with exit code {exit_code}")
        print(f"Output was saved to {log_file}")
        
    else:
        # Start the process
        process = subprocess.Popen(
            run_cmd,
            shell=True,                 # Needed for %SU2_RUN% expansion
            cwd=cfg_filepath,        # Set working directory to where the cfg file is
            stdout=subprocess.PIPE,     
            stderr=subprocess.STDOUT,   # Merge stderr into stdout
            text=True,                  # Decode bytes -> str automatically
            bufsize=1                   # Line-buffered
        )
    
        # Stream output and log simultaneously
        for line in process.stdout:
            print(line, end="")   # Show live in console
    
        # Wait for process to complete
        exit_code = process.wait()
    
    print(f"\nProcess finished with exit code {exit_code}")
            

if __name__ == "__main__":
    config_filename = "inv_NACA0012.cfg"
    config_filepath = "C:\\Users\\BriceM\\Documents\\Modules\\SU2_win64\\bin\\QuickStart"
    basic_CFD_run(config_filename, config_filepath, True)