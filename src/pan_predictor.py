import pandas as pd
import numpy as np
import subprocess
import re
import tempfile
import os
from multiprocessing import Pool




class PanPredictor:

    def __init__(self, path_to_netmhcpan, path_to_supergroups, tmpdir, n_cores=4):
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
        epitope, helper_val, path_to_netmhcpan, path_to_supergroups, tmpdir = args
        try:
            predictor = PanPredictor(path_to_netmhcpan, path_to_supergroups, tmpdir)
            return predictor._predict_row(epitope, helper_val)
        except:
            return (np.nan, np.nan, np.nan, np.nan)

    def _predict_row(self, epitope, helper_val):
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
