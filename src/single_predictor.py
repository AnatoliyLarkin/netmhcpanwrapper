import pandas as pd
import numpy as np
import subprocess
import re
import tempfile
import os
from multiprocessing import Pool


class SinglePredictor:

    def __init__(self, path_to_netmhcpan, path_to_mapping, tmpdir, n_cores=4):
        self.PATH_TO_NETMHCPAN = path_to_netmhcpan
        self.PATH_TO_MAPPING = path_to_mapping
        self.TMPDIR = tmpdir
        self.n_cores = n_cores


        self.ALPHABET  = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}


    def validate_peptide(self, epitope):
        return bool(re.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", epitope))

    def predict_affinnity(self, epitope, hla):
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
        epitope, hla, path_to_netmhcpan, tmpdir = args
        try:
            self = SinglePredictor(path_to_netmhcpan, None, tmpdir)
            return self.predict_affinnity(epitope, hla)
        except:
            return (np.nan, np.nan, np.nan)

    def predict_affinity_dataframe(self, df, epitope_colname, hla_colname):
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
