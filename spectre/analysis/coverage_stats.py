import numpy as np
import logging as logger


class CoverageStatistics(object):
    # should have "get" and "set" func but just use them directly
    def __init__(self, as_dev=False):
        # {"average": np.NaN, "std_dev": np.NaN, "min": np.NaN, "max": np.NaN}
        # Stats
        logger.basicConfig(level=logger.DEBUG) if as_dev else logger.basicConfig(level=logger.INFO)
        self.logger = logger
        self.chromosome_len = 0
        self.chromosome_name = ""
        self.average = np.NaN
        self.std_dev = np.NaN
        self.min = np.NaN
        self.max = np.NaN
        self.median = np.NaN

    def print(self):
        print_me = f'Statistics of chromosome {self.chromosome_name}' \
                   f'  chromosome length: {self.chromosome_len}\n' \
                   f'  average coverage: {self.average}\n' \
                   f'  median coverage: {self.median}\n' \
                   f'  standard deviation: {np.round(self.std_dev, 3)}\n' \
                   f'  min, max coverage: {np.round(self.min, 3)}, {np.round(self.max, 3)}\n'
        logger.info(print_me)


class CoverageData(object):
    # should have "get" and "set" func but just use them directly
    def __init__(self):
        # raw data
        self.coverage_raw = np.NaN
        self.positions = np.NaN
        self.coverage_log2 = np.NaN         # log2(x)
        self.normalized_cov = np.NaN        # x/median
        self.normalized_cov_ploidy = np.NaN   # x/median * 2 -> for diploid
