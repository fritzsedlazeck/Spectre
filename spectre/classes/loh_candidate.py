import itertools
import numpy as np
from spectre.util import logger
from spectre.classes.cnv_candidate import CNVCandidate


class LoHCandidate(CNVCandidate):
    def __init__(self, sample_origin=""):
        # from CNVCandidate
        CNVCandidate.__init__(self, sample_origin)
        # update
        self.type = "LOH"
        self.statistics = {"z-score": {"statistics": "nan", "score": "nan", "pvalue": "nan", "sample_score": "nan"}}
        # loh only
        self.aflist = []
        self.quallist = []
        self.median_cov = np.nan
        self.copy_number = np.nan
        self.pretty_size = ""
        self.number_snvs = 0
        self.gt_alt_prop = 0
        self.gt_alt_count = 0
        self.gt_all_count = 0
        # filtering
        self.loh_pass = False
        # stats
        self.af = np.nan

    def start_candidate_list(self, chromosome, pos, snv_depth, snv_af, snv_qual):
        self.chromosome = chromosome
        self.pos =[pos]
        self.start = pos
        self.end = pos
        self.size = 1
        self.size_kb = 0.0
        self.cov = [snv_depth]
        self.aflist = [snv_af]
        self.quallist = [snv_qual]
        self.is_init = True
        self.number_snvs = 1
        self.gt_alt_count = 1
        self.gt_all_count = 1


    def push_candidates(self, pos, snv_depth, snv_af, snv_qual):
        if pos is not np.nan:
            self.pos.append(pos)
            self.end = pos
            self.size = self.end - self.start
            self.size_kb = int(self.size / 1000)
            self.number_snvs += 1
            self.aflist.append(snv_af)
            self.quallist.append(snv_qual)
            self.gt_alt_count += 1
            self.gt_all_count += 1
        if snv_depth is not np.nan:
            self.cov.append(snv_depth)

    def combine(self, combine_list, cn_neutral, ploidy, snv_vcf_handler, snv_qual_filter):
        # shared elements
        lohc = combine_list[0]
        self.chromosome = lohc.chromosome
        self.start = lohc.start
        self.is_init = True
        self.id = lohc.id
        # merged elements
        self.pos = list(itertools.chain.from_iterable([lohc.pos for lohc in combine_list]))
        self.cov = list(itertools.chain.from_iterable([lohc.cov for lohc in combine_list]))
        self.aflist = list(itertools.chain.from_iterable([lohc.aflist for lohc in combine_list]))
        self.quallist = list(itertools.chain.from_iterable([lohc.quallist for lohc in combine_list]))
        self.number_snvs = sum([lohc.number_snvs for lohc in combine_list])
        # end
        lohc = combine_list[-1]
        self.end = lohc.end
        self.size = self.end - self.start
        self.size_kb = int(self.size / 1000)
        # copy number_snvs
        self.set_copy_number(cn_neutral, ploidy)
        # set gt_alt_prop and getting a bool for updating the
        self.set_gt_alt_prop(snv_vcf_handler, snv_qual_filter)
        #self.set_gt_alt_prop2(snv_vcf_handler, snv_qual_filter, combine_list)

    def set_gt_alt_prop(self, vcf_handler, qsnv):
        # compute gt_alt_prop
        gt_tag = "GT"
        gt_alt = (1, 1)
        loh_region = f'{self.chromosome}:{self.start}-{self.end}'
        gt_alt_sum, gt_all_sum = 0, 0
        for snv_record in vcf_handler.fetch(region=loh_region):
            gt_alt_sum += 1 if (snv_record.samples.values()[0][gt_tag] == gt_alt and \
                snv_record.qual > qsnv) else 0
            gt_all_sum += 1
        self.gt_alt_prop = gt_alt_sum/float(gt_all_sum)

    def set_gt_alt_prop2(self, vcf_handler, qsnv, combine_list):
        # set self.gt_all_count, self.gt_alt_count, self.gt_alt_prop
        gt_tag = "GT"
        gt_alt = (1, 1)
        gt_alt_count = sum([lohc.gt_alt_count for lohc in combine_list])
        gt_all_count = sum([lohc.gt_all_count for lohc in combine_list])
        gt_alt_sum, gt_all_sum = 0, 0
        for cindex in range(0, len(combine_list)-1):
            lohc = combine_list[cindex]
            lohcn = combine_list[cindex + 1]
            loh_region = f'{self.chromosome}:{lohc.end + 1}-{lohcn.start - 1}'
            for snv_record in vcf_handler.fetch(region=loh_region):
                gt_alt_sum += 1 if (snv_record.samples.values()[0][gt_tag] == gt_alt and \
                    snv_record.qual > qsnv) else 0
                gt_all_sum += 1
        gt_alt_count += gt_alt_sum
        gt_all_count += gt_all_sum
        self.gt_all_count = gt_all_count
        self.gt_alt_count = gt_alt_count
        self.gt_alt_prop = gt_alt_count/float(gt_all_count)

    def set_copy_number(self, cn_neutral, ploidy):
        self.median_cov = np.nanmedian(self.cov)
        self.copy_number = round(self.median_cov/cn_neutral*ploidy,2)
        self.median_cov_norm = self.copy_number
        self.cn_status = int(round(self.median_cov/cn_neutral*ploidy,0))

    def print(self, cn_neutral=None, ploidy=None):
        # bed-like format
        self.pretty_dna_span()
        cn = round(np.nanmedian(self.cov)/cn_neutral*ploidy,2) if cn_neutral is not None else "NA"
        return f'{self.chromosome}\t{self.start}\t{self.end}\t{self.id}\t{self.pretty_size}\t{self.size}\t' + \
               f'{np.nanmedian(self.cov)}x;{cn}'

    def vcf_line(self, cn_neutral, ploidy):
        # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  h2009_R8658
        id, ref, qual = ".", ".", "."
        alt = "<LOH>"
        filter = "PASS"
        cn = round(np.nanmedian(self.cov)/cn_neutral*ploidy,2)
        info = f'END={self.end};SVLEN={self.size};SVTYPE=LOH;CN={cn}'
        format = "GT:HO"
        sample = "1/1:NULL"
        return f'{self.chromosome}\t{self.start}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{sample}'

    def pretty_dna_span(self, dna_len=None):
        dna_1kb = 1e3
        dna_1mb = 1e6
        if dna_len is None:
            dna_len = self.size
            self.pretty_size = f'{dna_len}' if dna_len < dna_1kb else (f'{np.round(dna_len/dna_1kb,2)}kb' \
                if dna_len < dna_1mb else f'{np.round(dna_len/dna_1mb,2)}mb')
        else:
            return f'{dna_len}' if dna_len < dna_1kb else (f'{np.round(dna_len/dna_1kb,2)}kb' \
                if dna_len < dna_1mb else f'{np.round(dna_len/dna_1mb,2)}mb')


class MergeCNVLoH(object):
    def __init__(self, chr_list, loh, cnv):
        self.chr_list = chr_list
        self.loh = loh
        self.cnv = cnv
        self.merged_candidates = {}
        self.error = False
        self.logger = logger.setup_log(__name__, True)
        # autorun
        self.merge_candidates()

    def merge_candidates(self):
        if self.loh is None:
            cnv_len = sum([len(self.cnv[contig]) for contig in self.cnv])
            self.logger.debug(f'CNV cand: {cnv_len}')
            self.merged_candidates = self.cnv
        elif len(self.cnv) == 0:
            loh_len = sum([len(self.loh[contig]) for contig in self.loh])
            self.logger.debug(f'LOH cand: {loh_len}')
            self.merged_candidates = self.loh
        else:
            self.logger.info(f'init merge list by chromosome name')
            for each_chr in self.chr_list:
                # init merge list by chromosome
                if each_chr not in self.merged_candidates:
                    self.merged_candidates[each_chr] = []
            for each_chr in self.chr_list:
                self.logger.info(f'MERGE: loh candidates chromosome {each_chr}')
                loh_index = 0
                cnv_index = 0
                failsafe = len(self.loh[each_chr]) + len(self.cnv[each_chr]) + 1
                # no candidates in either
                if len(self.loh[each_chr]) == 0 and len(self.cnv[each_chr]) == 0:
                    self.logger.debug(f'no candidates for chromosome {each_chr}')
                    pass
                # no candidates in loh
                elif len(self.loh[each_chr]) == 0:
                    self.logger.debug(f'no LOH candidates for chromosome {each_chr} => CNV {len(self.cnv[each_chr])}')
                    while cnv_index < len(self.cnv[each_chr]):
                        this_cnv = self.cnv[each_chr][cnv_index]
                        self.merged_candidates[each_chr].append(this_cnv)
                        cnv_index += 1
                        failsafe -= 1
                        if failsafe < 0:
                            assert()
                # no candidates in cnv
                elif len(self.cnv[each_chr]) == 0:
                    self.logger.debug(f'no CNV candidates for chromosome {each_chr} => LOH {len(self.loh[each_chr])}')
                    # ID creator for LoH
                    while loh_index < len(self.loh[each_chr]):
                        this_loh = self.loh[each_chr][loh_index]
                        self.merged_candidates[each_chr].append(this_loh)
                        loh_index += 1
                        failsafe -= 1
                        if failsafe < 0:
                            assert()
                else:
                    self.logger.debug(f'CNV and LOH candidates for chromosome {each_chr}')
                    while loh_index < len(self.loh[each_chr]) or cnv_index < len(self.cnv[each_chr]):
                        # always same chromsome
                        # any is done
                        if loh_index >= len(self.loh[each_chr]) or cnv_index >= len(self.cnv[each_chr]):
                            # loh done, only push cnv
                            if loh_index >= len(self.loh[each_chr]):
                                this_cnv = self.cnv[each_chr][cnv_index]
                                self.merged_candidates[each_chr].append(this_cnv)
                                cnv_index += 1
                            # cnv done, only push loh
                            elif cnv_index >= len(self.cnv[each_chr]):
                                this_loh = self.loh[each_chr][loh_index]
                                self.merged_candidates[each_chr].append(this_loh)
                                loh_index += 1
                            else:
                                self.logger.info(f'MERGE: WARNING during merge, both indices over the limit')

                        else:
                            # loh ahead pos
                            this_loh = self.loh[each_chr][loh_index]
                            this_cnv = self.cnv[each_chr][cnv_index]
                            if this_loh.start >= this_cnv.start:
                                if cnv_index < len(self.cnv[each_chr]):
                                    self.merged_candidates[each_chr].append(this_cnv)
                                    cnv_index += 1
                            # cnv ahead pos
                            elif this_loh.start < this_cnv.start:
                                if loh_index < len(self.loh[each_chr]):
                                    self.merged_candidates[each_chr].append(this_loh)
                                    loh_index += 1
                            else:
                                self.logger.info(f'MERGE: WARNING during merge, none is ahead')
                        failsafe -= 1
                        if failsafe <= 0:
                            self.logger.info(f'MERGE: ERROR loop over the limit')
                            self.error = True
                            break
