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
    
    """Entry point that dispatches peptide–HLA affinity prediction in a chosen mode.

    Validates environment prerequisites (tmp directory, dataset files, core count,
    netMHCpan binary) before instantiating the appropriate predictor and running
    predictions on the supplied DataFrame.

    Args:
        prediction_mode (str): Either 'pan' (pan-HLA search across all supergroups)
            or 'single' (prediction against a user-supplied specific allele).
        path_to_netmhcpan (str): Path to netMHCpan directory.
        df (pd.DataFrame): Input DataFrame containing at minimum the epitope column
            and, depending on mode, an HLA or supergroup column.
        epitope_colname (str): Name of the column containing peptide sequences.
        supergroup_column (str, optional): For 'pan' mode — name of the column
            containing HLA supergroup labels to narrow the allele search.
            Defaults to None (full pan-HLA search).
        hla_column (str, optional): For 'single' mode — name of the column
            containing raw HLA allele strings to be normalised and predicted against.
            Defaults to None.

    Returns:
        pd.DataFrame: The input DataFrame with prediction result columns appended
            (columns vary by mode; see PanPredictor and SinglePredictor for details).

    Raises:
        ValueError: If dataset files are missing, N_CORES is not configured,
            netMHCpan cannot be initialised, or an unknown prediction_mode is given.
    """

    
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
    




    


