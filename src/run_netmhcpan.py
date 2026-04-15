import numpy as np
import sys
import subprocess
from src.pan_predictor import PanPredictor
from src.single_predictor import SinglePredictor
import pandas as pd
from pathlib import Path
import os
from src.config import PATH_TO_MAPPING, PATH_TO_SUPERGROUPS, TMPDIR, N_CORES
#from src.config import PATH_TO_NETMHCPAN




def run_netmhcpan(prediction_mode: str,
                  path_to_netmhcpan: str,
                  df: pd.DataFrame,
                  epitope_colname: str,

                  supergroup_column: str = None, #for panprediction
                  hla_column: str = None        #for single prediction
                  ):
    
    if not (os.path.exists(TMPDIR) or os.path.isdir(TMPDIR)):
        subprocess.run(f'mkdir {TMPDIR}', shell = True)

    if not (os.path.exists(PATH_TO_MAPPING) or os.path.exists(PATH_TO_SUPERGROUPS)):
        raise ValueError('datasets/mapping or datasets/supergroups files not found. Ensure that original repo structure is not violated or configure pathways in config.py')
    
    if not N_CORES:
        raise ValueError('Number of cores to compute is undefined. Please, configure N_CORES value in config.py')
    
    
    result = subprocess.run(f'{path_to_netmhcpan}/netMHCpan -h', shell = True, capture_output=True, text = True)
    if result.returncode != 0:
        raise ValueError('Could not initialize netMHCpan. Adjust PATH_TO_NETMHCPAN in config.py and ensure that your netMHCpan is properly configured')
    
    if prediction_mode == 'pan':

        predictor = PanPredictor(
            path_to_netmhcpan=path_to_netmhcpan,
            path_to_supergroups=PATH_TO_SUPERGROUPS,
            tmpdir=TMPDIR,
            n_cores=N_CORES,
        )

        if supergroup_column:
            df = predictor.pan_hla_matching_prediction(df, epitope_colname, supergroup_column)
        else:
            df = predictor.pan_hla_matching_prediction(df, epitope_colname)



    elif prediction_mode == 'single':

        predictor = SinglePredictor(
            path_to_netmhcpan = path_to_netmhcpan,
            path_to_mapping=PATH_TO_MAPPING,
            tmpdir=TMPDIR,
            n_cores=N_CORES
        )

        df = predictor.match_hla(df, hla_column)
        df = predictor.predict_affinity_dataframe(df, epitope_colname, 'matched_hla')

    else:
        raise ValueError(f'Unknown mode {prediction_mode}. Please select either pan or single')
    
    return(df)




if __name__ == '__main__':
    run_netmhcpan()
    




    


