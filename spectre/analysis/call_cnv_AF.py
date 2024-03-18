import numpy as np
import pysam
from spectre.util import logger


class CNVCall(object):
    def __init__(self, genome_info: dict, user_args):
        self.logger = logger.setup_log(__name__, user_args.as_dev)
        self.as_dev = user_args.as_dev
        self.genome_info = genome_info
        self.output_directory = user_args.out_dir
        self.sample_id = user_args.sample_id
        # snv vcf
        self.snv_vcf = None
        # snv bed
        # self.snv_af_bed_bgzip = None
        # candidates coverage
        self.cnv_candidates_by_coverage = None
        # only user given chromosomes
        self.only_chr_list = user_args.only_chr_list
        # AF to CN state concordance
        # difference between the "expected AF" and "computed AF" for a given CN state
        self.af_diff_buffer_default = 0.1  # will be min(default, 20% min ratio); e.g. if 4:1, then min(0.5, 1/5 * 0.2)
        self.default_proportion_min_af = 0.3  # 30%
        # all AF are MAF meaning if AF > 0.5, then AF = 1 - AF
        self.cn_state_to_af = {
            0: 0,  # complete DEL then we should not have data, init = 0
            1: 0,  # 0 or 1, as we use MAF, then only 0
            3: 1 / 3,  # 2:1 ratio, as we use MAF then 1/3
            4: [1 / 4, 2 / 4],  # 3:1 or 2:2 ratio, then 1/4 or 2/4
            5: [1 / 5, 2 / 5],  # 4:1 or 3:2 ratio, then 1/5 or 2/5
            6: [1 / 6, 2 / 6, 3 / 6],  # 5:1 or 4:2 or 3:3 ratio, then 1/6 or 2/6 or 3/6
        }
        self.user_ploidy = user_args.ploidy

    # CN states and CNV calls
    def get_CNV_type_from_af(self, af, slack):
        maf = af if af < 0.5 else 1 - af  # minor allele frequency
        upper_boundaries = 0.7
        lower_boundaries = 0.3
        if maf == 0.0 or maf == 1.0:
            return "DEL"
        elif upper_boundaries - slack < maf < upper_boundaries + slack or \
            lower_boundaries - slack < maf < lower_boundaries + slack:
            # if maf is beteen 0.7+/-slack (upper_boundaries) or 0.3+/-slack DUP
            return "DUP"
        else:
            return ""

    def check_af_duplication(self, cn, maf, slack) -> bool:
        """
        Determines if the combination of copy number and minor allele frequency corresponds to a duplication.
        :param cn: copy number
        :param maf: minor allele frequency
        :param slack: +/- off set from target allele frequency
        :return: true = Duplication
        """
        dr = cn % 2
        # check for whole chromosome duplication
        if dr == 0:  # duplication residual
            if 0.5 - slack < maf < 0.5 + slack:  # one copy got duplicated
                return True
            else:
                return False

        # calculate upper and lower boundaries
        lower_target_af = (dr / cn)  # CN = 3 -> lower_target_af = 1/3
        upper_target_af = ((cn - dr) / cn)  # CN = 3 -> upper_target_af = 2/3

        if lower_target_af - slack < maf < lower_target_af + slack or \
            upper_target_af - slack < maf < upper_target_af + slack:
            return True
        return np.nan()

    def af_cnv_call_region(self, cnv_candidates, snv_vcf):
        cnv_calls = {}
        vcf_file = pysam.VariantFile(snv_vcf)
        vcf_sample = list(vcf_file.header.samples).pop()  # only the first sample is used if a multisample is given
        for each_chromosome in cnv_candidates.keys():
            self.logger.debug(each_chromosome)
            cnv_calls[each_chromosome] = []
            for each_candidate in cnv_candidates[each_chromosome]:
                genomic_region = f'{each_chromosome}:{each_candidate.start}-{each_candidate.end}'
                # CN state == pliody is a FP as there is no copy number alteration
                # NOTE: should we generalize it ot  cn_state == ploidy?
                if each_candidate.cn_status == 2 and self.user_ploidy == 2:
                    self.logger.debug(f'FP,{genomic_region} => {each_candidate.cn_status},{each_candidate.type}')
                else:
                    # AF is already a float
                    af = np.nanmean([var_record.samples[vcf_sample]["AF"] for var_record in \
                        vcf_file.fetch(region=genomic_region)])
                    if self.af_cn_state_concordance(af, each_candidate.cn_status):
                        cnv_calls[each_chromosome].append(each_candidate)
                    else:
                        self.logger.debug(f'FP,{genomic_region} => {each_candidate.cn_status},{each_candidate.type},{af}')
        return cnv_calls

    def af_cn_state_concordance(self, af, cn_state):
        def af_cn_concordance(_af, _e_af, _af_diff_threshold):
            self.logger.debug([abs(_af - _e_af), _af_diff_threshold])
            if abs(_af - _e_af) <= _af_diff_threshold:
                return 1
            else:
                return 0

        score = 0
        if self.user_ploidy != 2:
            return True
        else:
            if int(cn_state) > 6:
                # cn_state = "6"
                return True
            else:
                expected_af = self.cn_state_to_af[cn_state]
                if isinstance(expected_af, list):
                    for each_expected_af in expected_af:
                        af_diff_threshold = each_expected_af * self.default_proportion_min_af
                        score += af_cn_concordance(af, each_expected_af, af_diff_threshold)
                else:
                    af_diff_threshold = expected_af * self.default_proportion_min_af if cn_state > 1 \
                        else self.af_diff_buffer_default
                    score += af_cn_concordance(af, expected_af, af_diff_threshold)
                # if score == 0 then it is a false positive (FP) call, else we cannot decide, adn we do not decide
            self.logger.debug(f'score: {score}')
        return score != 0
