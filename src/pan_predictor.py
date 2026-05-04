import pandas as pd
import numpy as np
import subprocess
import re
import tempfile
import os
from multiprocessing import Pool




class PanPredictor:

    """Predicts MHC-I peptide binding affinity across all HLA supergroups using netMHCpan.

    Implements a two-round pan-HLA prediction strategy: first identifies the best
    HLA supergroup via representative alleles, then finds the best individual allele
    within that supergroup. Supports multiprocessing for batch predictions. 

    Attributes:
        PATH_TO_NETMHCPAN (Path): Path to the netMHCpan executable directory.
        PATH_TO_SUPERGROUPS (Path): Path to the HLA supergroups definition file.
        TMPDIR (Path): Directory for temporary input/output files.
        n_cores (int): Number of parallel worker processes.
        _supergroups (dict | None): Cached supergroup data after first load.
    """


    def __init__(self, path_to_netmhcpan, path_to_supergroups, tmpdir, n_cores=4):

        """Initializes PanPredictor with paths and runtime configuration.

        Args:
            path_to_netmhcpan (Path | str): Path to the netMHCpan executable directory.
            path_to_supergroups (Path | str): Path to the HLA supergroups file.
            tmpdir (Path | str): Directory for temporary files.
            n_cores (int, optional): Number of parallel worker processes
        """



        self.PATH_TO_NETMHCPAN = path_to_netmhcpan
        self.PATH_TO_SUPERGROUPS = path_to_supergroups
        self.TMPDIR = tmpdir
        self.n_cores = n_cores
        self._supergroups = None


    def _load_supergroups(self):

        supergroups = {}
        with open(self.PATH_TO_SUPERGROUPS, mode='r') as f:
            for string in f:
                sg_label, representative, other_alleles = string.split(';')
                supergroups[sg_label] = {'representative': representative}
                if other_alleles.strip() == '0':
                    supergroups[sg_label]['alleles'] = None
                else:
                    supergroups[sg_label]['alleles'] = other_alleles.strip().split(',')

        self._supergroups = supergroups
        return supergroups

    def validate_peptide(self, epitope):
        """Validates that a peptide sequence contains only canonical amino acid characters.

        Args:
            epitope (str): The peptide sequence to validate.

        Returns:
            bool: True if the sequence consists solely of standard amino acid
                single-letter codes, False otherwise.
        """

        return bool(re.match(r"^[ACDEFGHIKLMNPQRSTVWYX]+$", epitope))

    def predict_affinnity(self, epitope, hla):

        """Runs netMHCpan for a single peptide–HLA pair and parses the binding output.

        Writes the peptide to a temporary file, invokes netMHCpan in peptide mode
        with binding affinity output (-BA), then extracts Score_BA, %Rank_BA, and
        Aff(nM) from stdout. Temporary files are always cleaned up.

        Args:
            epitope (str): A validated peptide sequence.
            hla (str): HLA allele string in netMHCpan format (e.g. 'HLA-A0201').

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
        """Multiprocessing-compatible static worker that predicts affinity for one row.

        Instantiates a fresh PanPredictor and
        calls _predict_row. Returns NaN values on any exception to avoid killing
        the worker pool.

        Args:
            args (tuple): A tuple of
                (epitope, helper_val, path_to_netmhcpan, path_to_supergroups, tmpdir).

        Returns:
            tuple[str | float, float, float, float]: A tuple of
                (best_allele, best_score_ba, best_%Rank_BA, best_Aff_nM),
                or (nan, nan, nan, nan) on failure.
        """
        epitope, helper_val, path_to_netmhcpan, path_to_supergroups, tmpdir = args
        try:
            predictor = PanPredictor(path_to_netmhcpan, path_to_supergroups, tmpdir)
            return predictor._predict_row(epitope, helper_val)
        except:
            return (np.nan, np.nan, np.nan, np.nan)

    def _predict_row(self, epitope, helper_val):
        """Core prediction logic for a single epitope with optional supergroup hint.

        Performs either one or two rounds of netMHCpan calls depending on whether
        a valid supergroup hint is provided:

        - No hint / invalid hint: Round 1 screens all supergroup representatives to
        find the best supergroup; Round 2 screens all alleles within that supergroup.
        - Valid hint: Screens only alleles within the specified supergroup (one round).

        The best allele is selected by jointly minimising %Rank_BA and Aff(nM).

        Args:
            epitope (str): The peptide sequence to predict.
            helper_val (str | None): A supergroup label to restrict the search, or
                None to perform full pan-HLA search.

        Returns:
            tuple[str | float, float | None, float, float]: A tuple of
                (best_allele, best_score_ba, best_%Rank_BA, best_Aff_nM),
                or (nan, nan, nan, nan) if no valid prediction is found.
        """

        supergroups = self._load_supergroups()

        best_rank = 1000
        best_aff = float('inf')
        best_allele = None
        best_score_ba = None

        # Case 1 - no helper value or invalid --> 2 rounds of prediction
        if (not helper_val) or (helper_val not in supergroups):

            best_supergroup = None
            best_supergroup_aff = float('inf')
            best_supergroup_rank = 1000
            best_supergroup_score_ba = None

            # 1st round - finding best supergroup by representative allele
            for sg in supergroups:
                repr = supergroups[sg]['representative']
                if '*' in repr:
                    repr = ''.join(repr.split('*'))
                try:
                    score_ba, percent_rank_ba, affinity_nm = self.predict_affinnity(epitope, repr)
                    if percent_rank_ba < best_supergroup_rank and affinity_nm < best_supergroup_aff:
                        best_supergroup_rank = percent_rank_ba
                        best_supergroup_aff = affinity_nm
                        best_supergroup = sg
                        best_supergroup_score_ba = score_ba
                except:
                    continue

            if best_supergroup is None:
                return (np.nan, np.nan, np.nan, np.nan)

            # 2nd round - finding best allele within best supergroup
            if not supergroups[best_supergroup]['alleles']:
                # singlet: representative is the only allele
                best_rank = best_supergroup_rank
                best_aff = best_supergroup_aff
                best_score_ba = best_supergroup_score_ba
                best_allele = supergroups[best_supergroup]['representative']

                if '*' in best_allele:
                    best_allele = ''.join(best_allele.split('*'))

            else:
                for a in supergroups[best_supergroup]['alleles']:
                    a = ''.join(a.split('*')) if '*' in a else a  # BUG FIX: outside try/except
                    try:
                        score_ba, percent_rank_ba, affinity_nm = self.predict_affinnity(epitope, a)
                        if percent_rank_ba < best_rank and affinity_nm < best_aff:
                            best_rank = percent_rank_ba
                            best_aff = affinity_nm
                            best_allele = a
                            best_score_ba = score_ba
                    except:
                        continue

        # Case 2 - valid helper value --> 1 round on all alleles in that supergroup
        else:
            repr = supergroups[helper_val]['representative']
            others = supergroups[helper_val]['alleles']
            alleles = others if others else [repr]

            for a in alleles:
                a = ''.join(a.split('*')) if '*' in a else a  # BUG FIX: outside try/except
                try:
                    score_ba, percent_rank_ba, affinity_nm = self.predict_affinnity(epitope, a)
                    if percent_rank_ba < best_rank and affinity_nm < best_aff:
                        best_rank = percent_rank_ba
                        best_aff = affinity_nm
                        best_allele = a
                        best_score_ba = score_ba
                except:
                    continue

        return (best_allele, best_score_ba, best_rank, best_aff)

    def pan_hla_matching_prediction(self, df, epitope_column, helper_column=None):

        """Runs pan-HLA affinity prediction on all rows of a DataFrame in parallel.

        Distributes rows across worker processes using multiprocessing.Pool.
        Appends prediction results as new columns to the input DataFrame.

        Args:
            df (pd.DataFrame): Input DataFrame containing at least the epitope column
                and optionally a supergroup helper column.
            epitope_column (str): Name of the column containing peptide sequences.
            helper_column (str, optional): Name of the column containing supergroup
                labels used to narrow the allele search. Defaults to None.

        Returns:
            pd.DataFrame: The input DataFrame with four new columns appended:
                - 'Allele' (str): Best-matching HLA allele.
                - '%Rank Score_BA' (float): Binding affinity score for that allele.
                - '%Rank_BA' (float): Percentile rank of binding affinity.
                - 'Aff(nM)' (float): Predicted binding affinity in nanomolar.
        """


        args = [
            (
                row[epitope_column],
                row[helper_column] if helper_column else None,
                self.PATH_TO_NETMHCPAN,
                self.PATH_TO_SUPERGROUPS,
                self.TMPDIR,
            )
            for _, row in df.iterrows()
        ]

        with Pool(processes=self.n_cores) as pool:
            results = pool.map(PanPredictor._worker, args)

        alleles, scores, ranks, affs = zip(*results)

        df['Allele'] = alleles
        df['%Rank Score_BA'] = scores
        df['%Rank_BA'] = ranks
        df['Aff(nM)'] = affs

        return df
