import pandas as pd
import numpy as np
import subprocess
import re
import tempfile
import os
from multiprocessing import Pool


class SinglePredictor:

    """Predicts MHC-I peptide binding affinity for a known, user-supplied HLA allele.

    Handles allele name normalisation via a nomenclature mapping file before
    calling netMHCpan, and supports parallel batch predictions over a DataFrame.

    Attributes:
        PATH_TO_NETMHCPAN (Path): Path to the netMHCpan executable directory.
        PATH_TO_MAPPING (Path): Path to the allele nomenclature mapping file.
        TMPDIR (Path): Directory for temporary files.
        n_cores (int): Number of parallel worker processes.
        ALPHABET (set[str]): Set of valid single-letter amino acid codes.
    """


    def __init__(self, path_to_netmhcpan, path_to_mapping, tmpdir, n_cores=4):

        """Initializes SinglePredictor with paths and runtime configuration.

        Args:
            path_to_netmhcpan (Path | str): Path to the netMHCpan directory.
            path_to_mapping (Path | str): Path to the allele nomenclature mapping file.
            tmpdir (Path | str): Directory for temporary files.
            n_cores (int, optional): Number of parallel worker processes.
        """


        self.PATH_TO_NETMHCPAN = path_to_netmhcpan
        self.PATH_TO_MAPPING = path_to_mapping
        self.TMPDIR = tmpdir
        self.n_cores = n_cores


        self.ALPHABET  = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}


    def validate_peptide(self, epitope):
        """Validates that a peptide sequence contains only canonical aminoacid characters.

        Args:
            epitope (str): The peptide sequence to validate.

        Returns:
            bool: True if the sequence consists solely of standard aminoacid
                single-letter codes, False otherwise.
        """

        return bool(re.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", epitope))

    def predict_affinnity(self, epitope, hla):
        """Runs netMHCpan for a single peptide–HLA pair and parses the binding output.

        Writes the peptide to a temporary file, invokes netMHCpan in peptide mode
        with binding affinity output (-BA), then extracts Score_BA, %Rank_BA, and
        Aff(nM) from stdout. Temporary files are always cleaned up.

        Args:
            epitope (str): A validated peptide sequence.
            hla (str): HLA allele in netMHCpan format (e.g. 'HLA-A0201').

        Returns:
            tuple[float, float, float]: A tuple of
                (Score_BA, %Rank_BA, Aff(nM)).

        Raises:
            ValueError: If epitope validation fails or binding data cannot be
                extracted from netMHCpan output.
        """

        if not self.validate_peptide(epitope):
            raise ValueError(f'Epitope validation for {epitope} failed')

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', dir=self.TMPDIR, delete=False) as f:
            f.write(epitope)
            tmp_input = f.name
        tmp_output = tmp_input.replace('.txt', '_affinity.xls')

        try:
            res = str(subprocess.run(
                f'{self.PATH_TO_NETMHCPAN}/netMHCpan -a {hla} -p {tmp_input} -xls -xlsfile {tmp_output} -BA -mode 0',
                shell=True, check=True, capture_output=True
            ).stdout.decode('utf-8'))

            match = re.search(r'Score_BA\s+([\d.]+)\s+%Rank_BA\s+([\d.]+)\s+Aff\(nM\)\s+([\d.]+)(?:\s+<= [A-Z]+)?', res)

            if not match:
                lines = res.split('\n')
                for line in lines:
                    if epitope in line:
                        numbers = re.findall(r'\d+(?:\.\d+)?', line)
                        if len(numbers) >= 3:
                            score_ba = float(numbers[-3])
                            percent_rank_ba = float(numbers[-2])
                            affinity_nm = float(numbers[-1])
                            return (score_ba, percent_rank_ba, affinity_nm)

                raise ValueError(f'Could not extract binding data for epitope {epitope}')

            score_ba = float(match.group(1))
            percent_rank_ba = float(match.group(2))
            affinity_nm = float(match.group(3))

            return (score_ba, percent_rank_ba, affinity_nm)

        finally:
            for path in [tmp_input, tmp_output]:
                if os.path.exists(path):
                    os.remove(path)

    @staticmethod
    def _worker(args):
        """Multiprocessing-compatible static worker for a single peptide–HLA prediction.

        Instantiates a fresh SinglePredictor (required for subprocess safety) and
        calls predict_affinnity. Returns NaN values on any exception.

        Args:
            args (tuple): A tuple of
                (epitope, hla, path_to_netmhcpan, tmpdir).

        Returns:
            tuple[float, float, float]: A tuple of
                (Score_BA, %Rank_BA, Aff(nM)),
                or (nan, nan, nan) on failure.
        """
        epitope, hla, path_to_netmhcpan, tmpdir = args
        try:
            self = SinglePredictor(path_to_netmhcpan, None, tmpdir)
            return self.predict_affinnity(epitope, hla)
        except:
            return (np.nan, np.nan, np.nan)

    def predict_affinity_dataframe(self, df, epitope_colname, hla_colname):

        """Runs affinity prediction for all epitope–HLA pairs in a DataFrame in parallel.

        Zips the epitope and HLA columns into worker arguments, distributes them
        across a multiprocessing pool, and appends results as new columns.

        Args:
            df (pd.DataFrame): Input DataFrame containing epitope and HLA columns.
            epitope_colname (str): Name of the column containing peptide sequences.
            hla_colname (str): Name of the column containing HLA allele strings
                (should already be normalised via match_hla).

        Returns:
            pd.DataFrame: The input DataFrame with three new columns appended:
                - '%Rank Score_BA' (float): Binding affinity BA score.
                - '%Rank_BA' (float): Percentile rank of binding affinity.
                - 'Aff(nM)' (float): Predicted binding affinity in nanomolar.
        """


        args = [
            (epitope, hla, self.PATH_TO_NETMHCPAN, self.TMPDIR)
            for epitope, hla in zip(df[epitope_colname], df[hla_colname])
        ]

        with Pool(processes=self.n_cores) as pool:
            results = pool.map(SinglePredictor._worker, args)

        scores, ranks, affs = zip(*results)

        df['%Rank Score_BA'] = scores
        df['%Rank_BA'] = ranks
        df['Aff(nM)'] = affs

        return df

    def match_hla(self, df, hla_colname):
        """Normalises HLA allele names in a DataFrame column using a mapping file.

        Reads a semicolon-delimited mapping file where each line maps a canonical
        allele name to a set of synonymous aliases. For each query allele, checks
        whether it is a canonical key or a known alias, then writes the normalised
        netMHCpan-compatible allele name (asterisk removed) to a new column.
        Unresolved alleles are set to NaN.

        Args:
            df (pd.DataFrame): Input DataFrame containing the HLA column.
            hla_colname (str): Name of the column containing raw HLA allele strings.

        Returns:
            pd.DataFrame: The input DataFrame with a new 'matched_hla' column
                containing normalised allele strings, or NaN where no match is found.
        """


        mapping = {}
        matched_hlas = []

        with open(self.PATH_TO_MAPPING, mode='r') as mp:
            for i in mp:
                line = i.split(';')
                mapping[line[0]] = set(line[1].rstrip('\n').split(','))

        mapping_keys = set(mapping.keys())

        for query in df[hla_colname].tolist():
            if query in mapping_keys:
                matched_hlas.append(''.join(query.split('*')))
            else:
                found = False
                for h in mapping:
                    if query in mapping[h]:
                        matched_hlas.append(''.join(h.split('*')))
                        found = True
                        break
                if not found:
                    matched_hlas.append(np.nan)

        df['matched_hla'] = matched_hlas
        return df
