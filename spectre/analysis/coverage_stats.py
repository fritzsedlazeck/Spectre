import numpy as np


class CoverageStatistics(object):
    # should have "get" and "set" func but just use them directly
    def __init__(self):
        # {"average": np.NaN, "std_dev": np.NaN, "min": np.NaN, "max": np.NaN}
        # Stats
        self.chromosome_len = 0
        self.chromosome_name = ""
        self.average = np.NaN
        self.std_dev = np.NaN
        self.min = np.NaN
        self.max = np.NaN
        self.median = np.NaN
        self.genome_wide_median = np.NaN
        self.genome_wide_mean = np.NaN
        self.coverage_overwrite = np.NaN
        self.empty = True

    def print(self):
        print_me = f'Statistics of chromosome {self.chromosome_name}' \
                   f'  chromosome length: {self.chromosome_len}\n' \
                   f'  average coverage: {self.average}\n' \
                   f'  median coverage: {self.median}\n' \
                   f'  genome wide median coverage: {self.genome_wide_median}\n' \
                   f'  standard deviation: {np.round(self.std_dev, 3)}\n' \
                   f'  min, max coverage: {np.round(self.min, 3)}, {np.round(self.max, 3)}\n' \
                   f'  coverage overwrite: {self.coverage_overwrite}'
        return print_me


class CoverageData(object):
    # should have "get" and "set" func but just use them directly
    def __init__(self):
        # raw data
        self.coverage_raw = np.NaN
        self.positions = np.NaN
        self.coverage_log2 = np.NaN         # log2(x)
        self.normalized_cov = np.NaN        # x/median
        self.normalized_cov_ploidy = np.NaN   # x/median * 2 -> for diploid
        self.empty = True
