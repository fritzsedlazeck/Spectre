import numpy as np
import logging as logger


class CNVCandidate(object):
    def __init__(self, sample_origin="", as_dev=False):
        self.chromosome = ""
        self.start = 0
        self.end = 0
        self.size = 0
        self.pos = []
        self.cov = []
        self.type = ""
        self.cn_status = np.nan
        self.gt = "./."
        self.median_cov_norm = np.nan
        self.size_kb = np.nan
        self.het_score = float()  # makes sense from 0-1
        self.as_dev = as_dev
        logger.basicConfig(level=logger.DEBUG) if as_dev else logger.basicConfig(level=logger.INFO)
        self.logger = logger
        self.statistics = {}
        self.support_cnv_calls = {}
        self.id = ""
        self.sample_origin = sample_origin  # e.g. hg002
        self.merged_sample_references = set()
        self.median_cov_norm = 0.0

    def push_candidates(self, chromosome, cnv_cand_list, cnv_cov_cand_list, cnv_type):
        self.chromosome = chromosome
        self.pos = cnv_cand_list
        self.start = int(self.pos[0])
        self.end = int(self.pos[-1])
        self.size = self.end - self.start
        self.size_kb = int((self.end - self.start) / 1000)
        self.cov = cnv_cov_cand_list
        [self.cn_status, self.median_cov_norm] = self.set_copy_number_status(cnv_cov_cand_list)
        self.set_gt()
        self.type = cnv_type

    def add_candidates(self, cnv_cand_pos, cnv_cand_cov, cnv_id, cnv_merged_ids):
        self.pos = list(self.pos) + list(cnv_cand_pos)
        self.start = int(self.pos[0])
        self.end = int(self.pos[-1])
        self.size = self.end - self.start
        self.size_kb = int((self.end - self.start) / 1000)
        self.cov = list(self.cov) + list(cnv_cand_cov)
        self.merged_sample_references.add(cnv_id)
        self.merged_sample_references |= cnv_merged_ids  # "|" joins sets

    def set_gt(self):
        if self.cn_status == 0:
            self.gt = "1/1"
        elif self.cn_status == 1 or self.cn_status == 3:
            self.gt = "0/1"
        elif self.cn_status == 2:
            self.gt = "0/0"
        else:
            self.gt = "./."

    def median_coverage_candidates_merged(self):
        if len(self.cov) > 0:
            self.median_cov_norm = np.nanmedian(self.cov)
        else:
            self.median_cov_norm = 0.0

    def set_id(self, id: str):
        self.id = f"Spectre.{self.type}.{id}"

    def reinitialize_candidate_values(self) -> None:
        """
        Used after loading data from .spc file
        :return: None
        """
        self.cn_status = self.set_copy_number_status(self.cov)[0]  # set copy number status
        self.median_coverage_candidates_merged()  # median of normalized coverage
        self.merged_sample_references = set(self.merged_sample_references)  # convert list to set

    @staticmethod
    def set_copy_number_status(candidate_coverage):
        median_candidates_coverage = np.nanmedian(candidate_coverage)
        if median_candidates_coverage == np.inf or median_candidates_coverage == np.nan:
            return [2, 2.0]  # we are working with diploid data
            # one for regular signal ... inf can not be a valid CN state
        # copy number state is given by integer type conversion of the float value median coverage, e.g. 1.51-> 1
        return [int(round(median_candidates_coverage, 0)), median_candidates_coverage]

    # functions required to use this object as a key in a dictionary
    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return hash(self.id) == hash(other.id)
