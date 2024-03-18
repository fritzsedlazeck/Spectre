import numpy as np


class CNVCandidate(object):
    def __init__(self, sample_origin="",bin_size=1000):
        self.chromosome = ""
        self.start = 0
        self.end = 0
        self.size = 0
        self.bin_size = int(bin_size)
        self.pos = []
        self.cov = []
        self.normalization_value = 1.0
        # self.raw_cov = []
        self.type = ""
        self.cn_status = np.nan
        self.gt = "0/0"
        self.median_cov_norm = np.nan
        self.median_raw_cov = np.nan
        self.size_kb = np.nan
        self.is_init = False
        self.filter = "."
        self.quality = "."
        # logging
        self.id = ""
        #
        self.het_score = float()  # makes sense from 0-1
        self.statistics = {}
        self.support_cnv_calls = {}
        self.sample_origin = sample_origin  # e.g. hg002
        self.merged_sample_references = set()
        # SV support from SNFJ data
        self.sv_support = False

    def push_candidates(self, chromosome, cnv_pos_cand_list, cnv_cov_cand_list, cnv_type ):
        self.chromosome = chromosome
        self.pos = cnv_pos_cand_list
        self.start = int(self.pos[0])
        self.end = int(self.pos[-1])
        self.size = self.end - self.start
        self.size_kb = int((self.end - self.start) / self.bin_size)
        self.cov = cnv_cov_cand_list
        # self.raw_cov = cnv_raw_cov_cand_list
        [self.cn_status, self.median_cov_norm] = self.set_copy_number_status(cnv_cov_cand_list)
        self.median_raw_cov = np.nanmedian(self.median_cov_norm * self.normalization_value)
        self.set_gt()
        self.type = cnv_type
        self.is_init = True

    def add_candidates(self, cnv_cand_pos, cnv_cand_cov, cnv_id, cnv_merged_ids, push_front=False):

        if push_front:
            self.pos = list(cnv_cand_pos) + list(self.pos)
            self.cov = list(cnv_cand_cov) + list(self.cov)
        else:
            self.pos = list(self.pos) + list(cnv_cand_pos)
            self.cov = list(self.cov) + list(cnv_cand_cov)
        self.start = int(self.pos[0])
        self.end = int(self.pos[-1])
        self.size = self.end - self.start
        self.size_kb = int((self.end - self.start) / self.bin_size)
        self.median_cov_norm = np.nanmedian(self.cov)
        # self.raw_cov = list(self.raw_cov) + list(cnv_cand_raw_cov)
        self.median_raw_cov = self.median_cov_norm * self.normalization_value
        self.merged_sample_references.add(cnv_id)
        self.merged_sample_references |= cnv_merged_ids  # "|" joins sets

    def add_scaffold_candidate(self, scaf_start:int, scaf_end:int, scaf_cov:list, scaf_pos:list):
        # insert coverage and positions of scaffold into candidate between scaf_start and scaf_end in self.cov and self.pos

        # instert coverage
        self.cov = self.cov[:scaf_start] + scaf_cov + self.cov[scaf_end:]
        # insert positions
        self.pos = self.pos[:scaf_start] + scaf_pos + self.pos[scaf_end:]

        pass

    def set_gt(self):
        if self.cn_status == 0 or self.cn_status == 4:
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
