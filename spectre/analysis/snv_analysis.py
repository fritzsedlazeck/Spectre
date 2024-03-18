import pysam
import numpy as np
import math
import os
import sys
from spectre.util import logger
import gzip
import copy
from spectre.classes.loh_candidate import LoHCandidate
from spectre.libs.unidip.unidip import UniDip
from spectre.util.vcf_parser import VCFParser
from spectre.util.cnv_id import CNV_ID
from spectre.analysis.call_cnv_AF import CNVCall as CNVAnalysisSNP
from spectre.util.timing_dev import DevTiming

# loh_stats_file = './spectre-cnv/data/loh_stats_1000g_v1.json'
# json_file = open(loh_stats_file)
# loh_stats = json.load(json_file)


class SNVAnalysis(object):
    def __init__(self, snv_file=None, genome=None, args=None):
        # logger
        self.logger = logger.setup_log(__name__, args.as_dev)
        # as_dev in args.as_dev
        self.spectre_args = args
        # params
        self.genome_info = genome
        # self.coverage = coverage
        self.snv_file = snv_file
        self.file_check()
        self.snv_vcf = pysam.VariantFile(self.snv_file)
        # AF for CNVs
        self.snv_af_bed = None
        self.snv_af_df = None
        self.cnv_calls_list_af = None
        # arguments for loh
        if args is not None:
            self.loh_min_snv_perkb = args.loh_min_snv_perkb
            self.loh_min_snv_total = args.loh_min_snv_total
            self.loh_min_region_size = args.loh_min_region_size
            self.ploidy = args.ploidy
        else:
            # Default in case not provided
            self.loh_min_snv_perkb = 10
            self.loh_min_snv_total = 100
            self.loh_min_region_size = 100000.0
        # min values in case someone is funny
        if self.loh_min_snv_perkb < 1:
            self.logger.warning(f'Min acceptable value for "--loh-min-snv-perkb" is 1')
            self.loh_min_snv_perkb = 1
        if self.loh_min_snv_total < 10:
            self.logger.warning(f'Min acceptable value for "--loh-min-snv-total" is 10')
            self.loh_min_snv_total = 10
        if self.loh_min_region_size < 10000.0:
            self.logger.warning(f'Min acceptable value for "--loh-min-region-size" is 10000')
            self.loh_min_region_size = 10000.0
        # CN 2x
        self.perfect_het_depth = []
        self.cn_neutral_coverage_med = np.nan
        self.cn_neutral_coverage_mean = np.nan
        self.cn_neutral_coverage_stdev = np.nan
        # currently fixed
        self.max_abs_dev_2x = 0.2
        # snv
        self.gt_tag="GT"
        self.dp_tag="DP"
        af="AF"    # Clair3
        vaf="VAF"  # DeepVariant
        header=self.snv_vcf.header.formats.keys()
        self.header_has_tags = self.gt_tag in header and self.dp_tag in header and (af in header or vaf in header)
        self.logger.error("Not all tags needed are present in the SNV VCF file. Please remove the SNV file") if \
            not self.header_has_tags else None
        sys.exit(1)
        self.af_tag=af if af in header else vaf
        self.gt_het = (0, 1)
        self.gt_ref = (0, 0)
        self.gt_alt = (1, 1)
        # snv qual
        self.max_deviation_perfect_het_af = 0.12
        self.perfect_het_af = 0.5
        self.min_snv_qual = 5
        self.min_coverage_site = 3
        self.min_hom_alt_af = 0.75
        self.min_het_chr_sites = 0.4  # 40%
        self.min_gt_alt_prop = 0.85  # 85%
        #
        # LoH candidates
        self.loh_candidate_dict = {}
        self.loh_candidate_list = []
        self.merged_loh_candidate_list = []
        # Source for LoH size https://elifesciences.org/articles/70339
        # TODO: add mode for bed file in genes
        self.min_chr_length = 1e6

    def file_check(self):
        vcffile = VCFParser(as_dev=self.spectre_args.as_dev)
        new_snv_file = vcffile.vcf_has_index(self.snv_file)
        self.snv_file = new_snv_file
        
    def snv_multimodal_detect(self):
        self.het_dump()
        self.het_clean() if np.max(self.perfect_het_depth) > 10*np.nanmean(self.perfect_het_depth) else None
        detect_miltimodal = UniDip(self.perfect_het_depth)
        dist_detected = detect_miltimodal.run()
        # if more than 1 then multimodal detected, we expect the first peak is 2x
        if len(dist_detected) > 1:
            self.logger.info(f'Multimodal distributions found ({len(dist_detected)}) => {dist_detected}')
            # we expect that the first one is 2x, then 4x and so on
            dist_se = dist_detected[0]
            start, end = dist_se
            self.cn_neutral_coverage_med = np.nanmedian(self.perfect_het_depth[start:end])
            self.cn_neutral_coverage_mean = np.nanmean(self.perfect_het_depth[start:end])
            self.cn_neutral_coverage_stdev = np.nanstd(self.perfect_het_depth[start:end])

        elif 1 == len(dist_detected):
            self.logger.info(f'No multimodal distributions found ({len(dist_detected)}) => {dist_detected}')
            [dist_se] = dist_detected
            start, end = dist_se
            self.cn_neutral_coverage_med = np.nanmedian(self.perfect_het_depth[start:end])
            self.cn_neutral_coverage_mean = np.nanmean(self.perfect_het_depth[start:end])
            self.cn_neutral_coverage_stdev = np.nanstd(self.perfect_het_depth[start:end])
        else:
            # something happened / possible error, use all
            self.cn_neutral_coverage_med = np.nanmedian(self.perfect_het_depth)
            self.cn_neutral_coverage_mean = np.nanmean(self.perfect_het_depth)
            self.cn_neutral_coverage_stdev = np.nanstd(self.perfect_het_depth)

    def het_dump(self):
        het_dump_file = gzip.open(f'{self.spectre_args.out_dir}/debug/het_dump.txt.gz', 'wt')
        for het in self.perfect_het_depth:
            het_dump_file.write(f'{het}\n')
        het_dump_file.close()

    def het_clean(self):
        # remove top 2% coverage
        p98 = np.percentile(np.array(self.perfect_het_depth), 98)
        self.logger.debug(f'p98 HET filter = {p98}')
        tmp = copy.deepcopy(self.perfect_het_depth)
        self.perfect_het_depth = []
        for het in tmp:
            if het <= p98:
                self.perfect_het_depth.append(het)
        self.logger.debug(f'HET: {len(tmp)} |  p98 HET {len(self.perfect_het_depth)}')

    def _is_gt_het(self, snv):
        return snv[self.gt_tag] == self.gt_het

    def snv_nhet_chromosome(self):
        self.logger.debug(f'Het sites count per chromosome')
        chr_het = {chro: 0 for chro in self.genome_info["chromosomes"]}
        chr_alt = {chro: 0 for chro in self.genome_info["chromosomes"]}
        chr_all = {chro: 0 for chro in self.genome_info["chromosomes"]}
        cov_het = {chro: 0 for chro in self.genome_info["chromosomes"]}
        # go to begginging of file
        self.snv_vcf.seek(0)
        for snv_record in self.snv_vcf.fetch():
            try:
                [snv] = snv_record.samples.values()
            except ValueError:
                [snv] = snv_record.samples.values()[0]
            chr_all[snv_record.chrom] += 1
            if snv_record.qual > self.min_snv_qual and snv[self.dp_tag] > self.min_coverage_site:
                # We only expect a single sample column in the VCF
                # get only HET genotypes
                if self._is_gt_het(snv):
                    if abs(self.perfect_het_af - snv[self.af_tag]) < self.max_deviation_perfect_het_af:
                        chr_het[snv_record.chrom] += 1
                        cov_het[snv_record.chrom] += snv[self.dp_tag]
                else:
                    if snv[self.af_tag] >= self.min_hom_alt_af:
                        chr_alt[snv_record.chrom] += 1
        for chro in self.genome_info["chromosomes"]:
            chr_het_cov = cov_het[chro]/chr_het[chro] if chr_het[chro] > 0 else 0
            all_used = float(chr_alt[chro] + chr_het[chro])
            self.logger.debug(f'{chro}\t{chr_all[chro]}\t{chr_het[chro]}\t{chr_alt[chro]}\t{round(chr_het_cov,2)}x\t' \
                              f'{round(chr_het[chro]/all_used*100,2)}% het pass | ' \
                              f'{round(chr_het[chro]/chr_all[chro]*100,2)}% het all') if all_used > 0 else None
        self.logger.debug(f'Use chromsomes based on proportion of het sites')
        use_this_chr_for_norm = []
        for chro in self.genome_info["chromosomes"]:
            all_used = float(chr_alt[chro] + chr_het[chro])
            if all_used > 0:
                chr_het_cov = cov_het[chro]/chr_het[chro]
                if chr_het[chro]/all_used > self.min_het_chr_sites:
                    self.logger.debug(f'{chro}\t{round(chr_het_cov,2)}')
                    use_this_chr_for_norm.append(chro)
        return use_this_chr_for_norm

    def snv_copy_number_state(self, ploidy=2, return_result=False):
        # Dev timing ---------- #
        dev_timing = DevTiming(timing_process="SNV Copy number state compute", name=__name__)
        dev_timing.start()
        # -------------------- #
        self.ploidy = ploidy
        # dev only
        use_chr_het_sites = self.snv_nhet_chromosome()  # self.genome_info["chromosomes"]
        # Params and filters
        # we need a list of the "perfect heterozygous" meaning those with SF close to 50% AF
        # how close, +- 5%
        if len(self.perfect_het_depth) != 0:
            self.perfect_het_depth = []
        # go to begginging of file
        self.snv_vcf.seek(0)
        for snv_record in self.snv_vcf.fetch():
            # We only expect a single sample column in the VCF
            try:
                [snv] = snv_record.samples.values()
            except ValueError:
                [snv] = snv_record.samples.values()[0]
            if use_chr_het_sites is not None:
                if snv_record.chrom in use_chr_het_sites:
                    # Only high quality calls, threshold suggested by Medhat
                    if snv_record.qual > self.min_snv_qual:
                        # get only HET genotypes
                        if self._is_gt_het(snv):
                            if abs(self.perfect_het_af - snv[self.af_tag]) < self.max_deviation_perfect_het_af:
                                self.perfect_het_depth.append(snv[self.dp_tag])
            else:
                self.logger.error(f'use_chr_het_sites = {use_chr_het_sites}')
        self.logger.debug(f'n het sites used for normalization = {len(self.perfect_het_depth)}')
        # detect multi-modal dist
        self.snv_multimodal_detect()
        # Dev timing ---------- #
        dev_timing.end()
        # -------------------- #
        self.logger.info(f'Copy number neutral at ploidy {self.ploidy} is = {self.cn_neutral_coverage_med}x')
        if return_result:
            return self.cn_neutral_coverage_med

    def call_cnv_af_region(self, cnv_calls_list):
        cnv_by_af = CNVAnalysisSNP(genome_info=self.genome_info,
                                   user_args=self.spectre_args)
        self.logger.info("Calculating CNV events based on SNV data")
        self.cnv_calls_list_af = cnv_by_af.af_cnv_call_region(cnv_calls_list, self.snv_file)

    @staticmethod
    def get_score_sigmoid(lohc):
        # https://en.wikipedia.org/wiki/Logistic_function
        # x is the length of the give loh region
        x = lohc.size
        max_score = 60
        midpoint = 1e6 # 1mb has score of 30
        return min(max_score, max_score/(1+math.exp(-1*(x/midpoint-1))))*np.nanmean(lohc.aflist)

    # LoH
    def candidate_loh(self, existing_cnv_ids):
        # Dev timing ---------- #
        dev_timing = DevTiming(timing_process="Candidate LoH selection", name=__name__)
        dev_timing.start()
        # -------------------- #
        loh_save_candidates = False
        # go to start of VCF
        self.snv_vcf.reset()
        # working variables
        max_skip = 2
        min_save = 10
        fails = {"qc":0, "gt":0}
        for contigs in self.snv_vcf.header.contigs.items():
            chro, chro_info = contigs
            if chro_info.length > self.min_chr_length and chro in self.spectre_args.only_chr_list:
                #init dict
                if chro not in self.loh_candidate_dict:
                    self.loh_candidate_dict[chro] = []
                # candidate start
                loh_candidates = LoHCandidate()
                total_sites = 0
                min_alt_sites_cand = 0.70
                not_loh_sites, loh_sites_count = 0, 0
                loh_per_cand, sites_per_cand = 0, 0
                # timitng per chr
                dev_timing_chr = DevTiming(timing_process=f' {chro}', name=__name__)
                dev_timing_chr.start()
                for snv_record in self.snv_vcf.fetch(region=chro):
                    if snv_record.qual > self.min_snv_qual:
                        total_sites += 1
                        [snv] = snv_record.samples.values()
                        dp, gt, af = int(snv[self.dp_tag]), snv[self.gt_tag], float(snv[self.af_tag])
                        if gt == self.gt_alt:
                            loh_sites_count += 1
                            loh_per_cand += 1
                            sites_per_cand += 1
                            if not loh_candidates.is_init:
                                loh_candidates.start_candidate_list(snv_record.contig, int(snv_record.pos),
                                                                    dp, af, float(snv_record.qual))
                            else:
                                if max_skip > not_loh_sites > 0:
                                    loh_sites_count += not_loh_sites
                                not_loh_sites = 0
                                loh_candidates.push_candidates(int(snv_record.pos), dp, af, float(snv_record.qual))
                        else:
                            fails["gt"] += 1
                            not_loh_sites += 1
                            if not_loh_sites >= max_skip:
                                # save if threshold pass
                                if len(loh_candidates.pos) >= min_save and loh_per_cand >= sites_per_cand * min_alt_sites_cand:
                                    loh_candidates.set_copy_number(self.cn_neutral_coverage_med, self.ploidy)
                                    new_ids = CNV_ID.n_id_generator(existing_ids=existing_cnv_ids)
                                    existing_cnv_ids = existing_cnv_ids + new_ids
                                    loh_candidates.set_id(new_ids.pop())
                                    loh_candidates.gt_alt_prop = loh_per_cand/float(sites_per_cand)
                                    #self.loh_candidate_list.append(loh_candidates)
                                    self.loh_candidate_dict[chro].append(loh_candidates)
                                # re-start
                                sites_per_cand, loh_per_cand = 0, 0
                                not_loh_sites = 0
                                loh_candidates = LoHCandidate()
                            else:
                                sites_per_cand += 1
                    else:
                        fails["qc"] += 1
                # one last dump
                if len(loh_candidates.pos) > min_save:
                    loh_candidates.set_copy_number(self.cn_neutral_coverage_med, self.ploidy)
                    new_ids = CNV_ID.n_id_generator(existing_ids=existing_cnv_ids)
                    existing_cnv_ids = existing_cnv_ids + new_ids
                    loh_candidates.set_id(new_ids.pop())
                    loh_candidates.gt_alt_prop = loh_per_cand/float(sites_per_cand)
                    # self.loh_candidate_list.append(loh_candidates)
                    self.loh_candidate_dict[chro].append(loh_candidates)
                # log
                self.logger.debug(f'{fails}')
                if total_sites > 0:
                    self.logger.info(f'Chromosome: {chro}\t'\
                                     f'Total:{total_sites}\tLoH:{loh_sites_count}({np.round((loh_sites_count/total_sites)*100, 2)}%)\t'\
                                     f'Rest:{total_sites-loh_sites_count}({np.round(((total_sites-loh_sites_count)/total_sites)*100, 2)}%)')
                else:
                    self.logger.info(f'Chromosome: {chro}\tNA')
                dev_timing_chr.end(add_sep=False)
        # Dev timing ---------- #
        dev_timing.end()
        # -------------------- #

    def candidate_loh_neighboring(self, current_loh_candidate_dict):
        def merge_dist(size1, size2):
            size_mean = np.mean([size1, size2])
            # compare to 1mb
            _1mb = 1000000
            # weighted 1. by log10 of size
            ## v1
            prop_threhold = 1/(np.log10(size_mean)**1.75) + ((np.log10(_1mb)/np.log10(size_mean))/100)
            size_merge_thr = size_mean*prop_threhold
            return size_merge_thr

        def merge_candidates(loh_candidate_list, this_cycle):
            new_list_candidates = []
            combine_list = []
            combine_list_ids = []
            merge_occurred = False
            # failed_merges = {}
            self.logger.debug(f' === {this_cycle} === ')
            for cand_index in range(0, len(loh_candidate_list)):
                next_candidate = cand_index + 1
                loh_cand = loh_candidate_list[cand_index]
                # first check if last iteration
                if next_candidate >= len(loh_candidate_list):
                    if 0 == len(combine_list):
                        new_list_candidates.append(loh_cand)
                    elif loh_cand.id in combine_list_ids:
                        loh_candidates = LoHCandidate()
                        loh_candidates.combine(combine_list, self.cn_neutral_coverage_med, self.ploidy,
                                               self.snv_vcf, self.min_snv_qual)
                        new_list_candidates.append(loh_candidates)
                    break
                loh_cand_next = loh_candidate_list[next_candidate]
                # candidate distance test
                candidate_neighboring_dist_threshold = merge_dist(loh_cand.size, loh_cand_next.size)
                candidates_dist = loh_cand_next.start - loh_cand.end
                test_candidates_dist = candidates_dist <= candidate_neighboring_dist_threshold
                # candidate proportion of ALTs
                loh_cand_region = f'{loh_cand.chromosome}:{loh_cand.start}-{loh_cand_next.end}'
                # merge_region = f'{loh_cand.chromosome}:{loh_cand.end+1}-{loh_cand_next.start-1}'
                test_gt_alt_prop = self.compute_gt_alt_prop(loh_cand_region) if test_candidates_dist else 0.0
                if test_candidates_dist and test_gt_alt_prop >= self.min_gt_alt_prop:
                    """
                    A === B
                    A,B <- to be linked
                    [A,B] init  | list(Z,A) <- add B
                    """
                    merge_occurred = True
                    if 0 == len(combine_list):
                        combine_list = [loh_cand, loh_cand_next]
                        combine_list_ids = [loh_cand.id, loh_cand_next.id]
                    else:
                        combine_list.append(loh_cand_next)
                        combine_list_ids.append(loh_cand_next.id)
                else:
                    """
                    A ===== B
                    if []    -> add A to final
                    if [A]   -> add A to final
                    if [Z,A] -> combine
                    """
                    if 0 == len(combine_list):
                        new_list_candidates.append(loh_cand)
                    elif 1 == len(combine_list):
                        loh_cand = combine_list.pop()
                        new_list_candidates.append(loh_cand)
                    elif loh_cand.id in combine_list_ids:
                        loh_candidates = LoHCandidate()
                        loh_candidates.combine(combine_list, self.cn_neutral_coverage_med, self.ploidy,
                                               self.snv_vcf, self.min_snv_qual)
                        new_list_candidates.append(loh_candidates)
                    else:
                        pass
                    combine_list = []
                    combine_list_ids = []
            return new_list_candidates, not(merge_occurred)

        # RUN
        new_candidate_list = []
        max_cycles = 500
        # Dev timing ---------- #
        dev_timing = DevTiming(timing_process=f'Candidate LoH neighboring', name=__name__)
        dev_timing.start()
        # -------------------- #
        for contig in current_loh_candidate_dict.keys():
            no_merges = False
            current_loh_candidate_list = current_loh_candidate_dict[contig]
            dev_timing_chro = DevTiming(timing_process=f'{contig}', name=__name__)
            dev_timing_chro.start()
            new_cand = []
            for cycle in range(0, max_cycles):
                if no_merges:
                    break
                if 0 == cycle:
                    prev = len(current_loh_candidate_list)
                    new_cand, no_merges = merge_candidates(current_loh_candidate_list, cycle)
                    self.logger.debug(f'merged {prev} => {len(new_cand)} candidates')
                else:
                    prev = len(new_cand)
                    new_cand, no_merges = merge_candidates(new_cand, cycle)
                    self.logger.debug(f'merged {prev} => {len(new_cand)} candidates')
            new_candidate_list += new_cand
            dev_timing_chro.end()
        dev_timing.end()
        return new_candidate_list

    def candidate_loh_filtering(self):
        # Dev timing ---------- #
        dev_timing = DevTiming(timing_process="Candidate LoH filtering", name=__name__)
        dev_timing.start()
        # -------------------- #
        prety_min = 3e9
        prety_max = 0
        loh_fail = {"total": 0, "size": 0, "nSNV": 0, "gt_snv": 0}
        loh_pass = 0
        max_score = 60
        if 0 != len(self.merged_loh_candidate_list):
            self.loh_candidate_list = self.merged_loh_candidate_list
        self.snv_vcf.seek(0)
        for cand_index in range(0, len(self.loh_candidate_list)):
            loh_cand = self.loh_candidate_list[cand_index]
            # Filtering starts
            # min size
            snv_perkb_bylen = int(loh_cand.size_kb)*self.loh_min_snv_perkb
            if loh_cand.size > self.loh_min_region_size:
                # Min number of SNVs to be used
                if loh_cand.number_snvs > snv_perkb_bylen or loh_cand.number_snvs > self.loh_min_snv_total:
                    # GT ALT threshold
                    if loh_cand.gt_alt_prop >= self.min_gt_alt_prop:
                        # Update the candidate pass filter in the saved list
                        self.loh_candidate_list[cand_index].loh_pass = True
                        self.loh_candidate_list[cand_index].filter = "PASS"
                        self.loh_candidate_list[cand_index].gt = "1/1"
                        self.loh_candidate_list[cand_index].het_score = np.round(np.nanmean(loh_cand.aflist),2)
                        self.loh_candidate_list[cand_index].quality = np.round(np.nanmean(loh_cand.quallist),2)
                        self.loh_candidate_list[cand_index].statistics["z-score"]["score"] = self.get_score_sigmoid(loh_cand)
                        self.loh_candidate_list[cand_index].statistics["z-score"]["sample_score"] = self.get_score_sigmoid(loh_cand)
                        loh_pass += 1
                        # update min and max size loh events
                        if loh_cand.size < prety_min:
                            prety_min = loh_cand.size
                        if loh_cand.size > prety_max:
                            prety_max = loh_cand.size
                    else:
                        self.loh_candidate_list[cand_index].filter = f'ALT GT propotion: {loh_cand.gt_alt_prop}|{self.min_gt_alt_prop}'
                        loh_fail["total"] += 1
                        loh_fail["gt_snv"] += 1
                else:
                    self.loh_candidate_list[cand_index].filter = f'Number of SVs: {loh_cand.number_snvs}|{snv_perkb_bylen}/kb|{loh_cand.size_kb}'
                    loh_fail["total"] += 1
                    loh_fail["nSNV"] += 1
            else:
                self.loh_candidate_list[cand_index].filter = f'LoH region size kb: {loh_cand.size_kb}|{loh_cand.number_snvs}|{snv_perkb_bylen}/kb'
                loh_fail["total"] += 1
                loh_fail["size"] += 1
        if len(self.loh_candidate_list) > 0:
            prety_min = loh_cand.pretty_dna_span(prety_min)
            prety_max = loh_cand.pretty_dna_span(prety_max)
        else:
            prety_min, prety_max = "NA", "NA"
        self.logger.info(f'Total number of loh candidates = {loh_pass + loh_fail["total"]}, after filtering:')
        self.logger.info(f'  PASS={loh_pass} | FAIL={loh_fail} | MIN_LEN={prety_min} MAX_LEN={prety_max}')
        # Dev timing ---------- #
        dev_timing.end()
        # -------------------- #

    def loh_dump(self, version=""):
        dump_file = f'{self.spectre_args.out_dir}/debug/loh_dump{version}.tsv'
        dump_file_bgz = f'{self.spectre_args.out_dir}/debug/loh_dump{version}.tsv.gz'
        loh_dump_file = open(dump_file, 'w')
        for lohc in self.loh_candidate_list:
            loh_dump_file.write(f'{lohc.print(self.cn_neutral_coverage_med, self.ploidy)}\t' \
                                f'{lohc.filter}\t{lohc.het_score}\t{np.round(lohc.gt_alt_prop, 3)}\n')
        loh_dump_file.close()
        pysam.tabix_compress(filename_in=dump_file, filename_out=dump_file_bgz, force=True)
        pysam.tabix_index(filename=dump_file_bgz, force=True, preset="bed")
        if os.path.isfile(dump_file_bgz) and os.stat(dump_file_bgz).st_size != 0:
            os.remove(dump_file)

    def loh_report(self):
        for lohc in self.loh_candidate_list:
            if lohc.loh_pass:
                self.logger.debug(lohc.print(self.cn_neutral_coverage_med, self.ploidy))

    def loh_pass_only(self):
        pass_only = {}
        for chro in self.spectre_args.only_chr_list:
            if chro not in pass_only:
                pass_only[chro] = []
        for lohc in self.loh_candidate_list:
            if lohc.loh_pass:
                pass_only[lohc.chromosome].append(lohc)
        for chro in pass_only:
            self.logger.debug(f'LOH candidates pass {chro}:  {len(pass_only[chro])}')
        return pass_only

    def loh(self, existing_cnv_ids):
        #  > LoH candidates (quick)
        self.candidate_loh(existing_cnv_ids)
        #  > LoH debug
        # self.loh_dump(version="_v1")
        #  > LoH candidates merge
        # self.merged_loh_candidate_list = self.candidate_loh_neighboring(self.loh_candidate_list)
        self.merged_loh_candidate_list = self.candidate_loh_neighboring(self.loh_candidate_dict)
        #  > LoH filtering
        self.candidate_loh_filtering()
        #  > LoH debug
        #self.loh_dump(version="_v2")
        self.loh_dump()

    def compute_gt_alt_prop(self, loh_region):
        # compute gt_alt_prop
        # uses self.snv_vcf, self.min_snv_qual, self.gt_tag, self.gt_alt
        gt_alt_sum, gt_all_sum = 0, 0
        for snv_record in self.snv_vcf.fetch(region=loh_region):
            gt_alt_sum += 1 if (snv_record.samples.values()[0][self.gt_tag] == self.gt_alt and \
                snv_record.qual > self.min_snv_qual) else 0
            gt_all_sum += 1
        return gt_alt_sum/float(gt_all_sum)
