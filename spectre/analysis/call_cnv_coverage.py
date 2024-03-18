import numpy as np
from spectre.classes.cnv_candidate import CNVCandidate
from spectre.util import logger


class CNVCall(object):
    def __init__(self, as_dev=False):
        self.as_dev = as_dev
        self.logger = logger.setup_log(__name__, self.as_dev)
        self.candidate_list = []
        self.min_run = 10  # default = 10; minimum number of bins to be considered a CNV

    def append_candidate(self, sample_origin, bin_size, chromosome_name, cnv_pos_cand_list, cnv_cov_cand_list,
                         current_cnv_type):
        if len(cnv_pos_cand_list) >= self.min_run:
            cnv_candidates = CNVCandidate(sample_origin=sample_origin, bin_size=bin_size)
            cnv_candidates.push_candidates(chromosome_name, cnv_pos_cand_list, cnv_cov_cand_list, current_cnv_type)
            self.candidate_list.append(cnv_candidates)

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
        self.candidate_list = []
        # local
        run = 0
        cnv_pos_cand_list = []
        cnv_cov_cand_list = []
        # cnv_raw_cov_cand_list = []
        run_current_pos = 0
        run_start = 0
        cnv_type = ""
        current_cnv_type = ""
        # lower_bound, upper_bound = 1.5, 2.5  # values outside the CN2
        if not isinstance(chromosome_coverage_data.normalized_cov_ploidy, np.ndarray) or \
                not isinstance(chromosome_coverage_data.positions, np.ndarray):
            self.logger.debug(type(chromosome_coverage_data.normalized_cov_ploidy))
            self.logger.debug(type(chromosome_coverage_data.positions))
            self.logger.warning(f'No data is available, check that the mosdepth.region.bed.gz file is not empty '
                              f'or that the selected chromosome (--only-chr <CHR>) exists')

        for (cov, pos, raw_cov) in zip(chromosome_coverage_data.normalized_cov_ploidy,
                                       chromosome_coverage_data.positions, chromosome_coverage_data.coverage_raw):
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
                            cnv_pos_cand_list.append(int(pos))
                            cnv_cov_cand_list.append(float(cov))
                            # cnv_raw_cov_cand_list.append(float(raw_cov))
                        # break run end of chromosome
                        else:
                            self.append_candidate(sample_origin, bin_size, chromosome_name, cnv_pos_cand_list,
                                                  cnv_cov_cand_list, current_cnv_type)
                            # restart run
                            cnv_pos_cand_list = [int(pos)]
                            cnv_cov_cand_list = [float(cov)]
                            # cnv_raw_cov_cand_list = [float(raw_cov)]
                            run_start = pos
                            current_cnv_type = cnv_type
                # break run
                else:
                    self.append_candidate(sample_origin, bin_size, chromosome_name, cnv_pos_cand_list,
                                          cnv_cov_cand_list, current_cnv_type)
                    # restart run
                    cnv_pos_cand_list = []
                    cnv_cov_cand_list = []
                    # cnv_raw_cov_cand_list = []
                    run_start = pos
                    current_cnv_type = ""
                run_current_pos = pos
        # Append last candidate if possible to the candidate list
        self.append_candidate(sample_origin, bin_size, chromosome_name, cnv_pos_cand_list,
                              cnv_cov_cand_list, current_cnv_type)
        return self.candidate_list
