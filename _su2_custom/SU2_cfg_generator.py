# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 22:17:18 2025

@author: BriceM
"""
import traceback 
from pathlib import Path


def dict_to_cfg(cfg_dict: dict, cfg_filename: str, cfg_filepath: str = "") -> bool:
    '''
    Converts a dictionary cfg_dict to a configuration file and saves it to cfg_filename in the
    directory cfg_filepath. If the filepath does not exist, it will create it.
    Returns whether it was successful or not.

    Parameters
    ----------
    cfg_dict : dict
        Dictionary containing all configuration parameters.
    cfg_filename : str
        Filename of the desired cfg file output (needs .cfg at the end), NO FILEPATH.
    cfg_filepath : str 
        Default is empty, constains as a str the desired final location of the config file.
    Returns
    -------
    success : bool
        True if no error occurs, False if exception is caught.

    '''
    try: 
        cfg_filepath = cfg_filepath + "\\" if cfg_filepath != "" else ""
        out_path = Path(cfg_filepath + cfg_filename)
        # Ensure parent directory exists (no-op if '.' or already present).
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with out_path.open("w", encoding="utf-8", newline="\n") as f:
            for key, value in cfg_dict.items():
                # Cast to str to avoid TypeErrors on non-strings.
                f.write(f"{str(key)}= {str(value)}\n")

        print(f"SU2 configuration file '{out_path}' created successfully.")
        success = True
    except Exception as e:
        print("\nError Occured")
        print(e)
        traceback.print_exc()
        success = False
    finally:
        # handle other things if needed
        return success 
       
        
   
if __name__ == "__main__":
    from base_config_params import cfg_params as params 
    filename = "test_config.cfg"
    filepath = "Configs"
    dict_to_cfg(params, filename, filepath)
    