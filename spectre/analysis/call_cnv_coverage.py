import numpy as np
from analysis.cnv_candidate import CNVCandidate
import pandas as pd
import logging as logger


class CNVCall(object):
    def __init__(self, as_dev=False):
        self.as_dev = as_dev
        logger.basicConfig(level=logger.DEBUG) if as_dev else logger.basicConfig(level=logger.INFO)
        self.logger = logger

    def cnv_coverage(self, chromosome_coverage_data, bin_size, chromosome_name, sample_origin="",
                     lower_bound: float = None,
                     upper_bound: float = None):
        # get chromosome_coverage_data ("self.genome_analysis"), for each chromosome, "cov_data:
        # * .position: list size n
        #   .coverage_raw = np.NaN
        #   .positions = np.NaN
        #   .coverage_log2 = np.NaN         # log2(x)
        #   .normalized_cov = np.NaN        # x/median
        # * .normalized_cov_ploidy = np.NaN   # x/median * 2 -> for diploid
        candidate_list = []
        # local
        min_run = 10
        run = 0
        cnv_pos_cand_list = []
        cnv_cov_cand_list = []
        run_current_pos = 0
        run_start = 0
        cnv_type = ""
        current_cnv_type = ""
        # lower_bound, upper_bound = 1.5, 2.5  # values outside the CN2
        if not isinstance(chromosome_coverage_data.normalized_cov_ploidy, np.ndarray) or \
                not isinstance(chromosome_coverage_data.positions, np.ndarray):
            self.logger.debug(type(chromosome_coverage_data.normalized_cov_ploidy))
            self.logger.debug(type(chromosome_coverage_data.positions))
            self.logger.error(f'No data is available, check that the mosdepth.region.bed.gz file is not empty '
                              f'or that the selected chromosome (--only-chr <CHR>) exists')

        for (cov, pos) in zip(chromosome_coverage_data.normalized_cov_ploidy, chromosome_coverage_data.positions):
            # start/continue run
            if not np.isnan(cov):
                # Determine CNV Type
                if cov > upper_bound or cov < lower_bound:
                    cnv_type = "DUP" if cov > upper_bound else "DEL"
                    # start run
                    if run_start == 0:
                        current_cnv_type = cnv_type
                        run_start = pos
                    else:
                        # continue run if not in the chromosome end
                        if pos - run_current_pos < bin_size + 1 and current_cnv_type == cnv_type:
                            run += 1
                            cnv_pos_cand_list.append(pos)
                            cnv_cov_cand_list.append(cov)
                        # break run end of chromosome
                        else:
                            if len(cnv_pos_cand_list) > 1:
                                cnv_candidates = CNVCandidate(sample_origin, self.as_dev)
                                cnv_candidates.push_candidates(chromosome_name, cnv_pos_cand_list, cnv_cov_cand_list,
                                                               current_cnv_type)

                                if len(cnv_pos_cand_list) >= min_run:
                                    candidate_list.append(cnv_candidates)
                            # restart run
                            run = 1
                            cnv_pos_cand_list = [pos]
                            cnv_cov_cand_list = [cov]
                            run_start = pos
                            current_cnv_type = cnv_type
                # break run
                else:
                    if len(cnv_pos_cand_list) > 1:
                        cnv_candidates = CNVCandidate(sample_origin, self.as_dev)
                        cnv_candidates.push_candidates(chromosome_name, cnv_pos_cand_list, cnv_cov_cand_list,
                                                       current_cnv_type)

                        if len(cnv_pos_cand_list) >= min_run:
                            candidate_list.append(cnv_candidates)
                    # restart run
                    run = 1
                    cnv_pos_cand_list = []
                    cnv_cov_cand_list = []
                    run_start = pos
                    current_cnv_type = ""
                run_current_pos = pos

        return candidate_list
