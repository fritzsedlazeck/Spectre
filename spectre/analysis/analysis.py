import os

import numpy as np
import pandas as pd
import pysam

from spectre.analysis.breakpoint_analysis import Breakpoints
from spectre.util import outputWriter
from spectre.util import logger
from spectre.analysis.call_cnv_coverage import CNVCall
from spectre.analysis.cnv_metrics import CNVMetrics
from spectre.analysis.cnv_post_processing import CNV_Post_Processing
from spectre.analysis.coverage_stats import CoverageData
from spectre.analysis.coverage_stats import CoverageStatistics
from spectre.classes.loh_candidate import MergeCNVLoH
from spectre.plots.plot import CNVPlot
from spectre.plots.plot import CoveragePlot
from spectre.util.cnv_id import CNV_ID
from spectre.util.metadata.metadataCollector import FastaRef


class CNVAnalysis(object):
    def __init__(self, coverage_file, outbed, outvcf, sample_id, spectre_args, metadata_ref, genome_info,
                 debug_dir=""):

        self.logger = logger.setup_log(__name__, spectre_args.as_dev)
        # dev/debug + hidden params
        self.as_dev = spectre_args.as_dev
        self.debug_dir = debug_dir
        # input
        self.spectre_args = spectre_args
        self.metadata_path = metadata_ref
        self.metadata = {}
        self.sample_id = sample_id
        self.output_bed = outbed
        self.output_vcf = outvcf
        self.genome_info = genome_info
        # coverage
        self.coverage_file = coverage_file  # s.coverage_file
        # self.mosdepth_data = coverage_files.mosdepth_summary_data  # has chr_mean_cover, genome_mean_cover, genome_bases
        self.genome_wide_median = 0.0
        self.genome_wide_mean = 0.0
        # population
        self.population_sample_ids = set()  # hold the filename of the population samples
        # work variables
        self.positions = np.zeros(0)
        self.coverage = np.zeros(0)
        self.coverage_log2 = np.zeros(0)
        self.bin_size = 0
        self.normalization_value = 1.0
        self.genome_std = 0.0
        # results by chromosome
        self.coverage_analysis = {}  # use this for CNV call NOTE: here is the coverage data under "cov_data"
        self.cnv_calls_list = {}  # use this for CNV calls
        self.cnv_calls_list_af_filtered = None  # use if not None
        self.raw_cnv_calls_list = {}  # use this for storing raw CNV call
        self.intermediate_candidates_file_location = ""  # holds output path of serialized cnv object
        self.cnv_merged = {}  # use this for CNV call
        self.existing_cnv_ids = []  # holds all already used cnv IDs
        # cnv metrics
        self.cnv_metrics = None
        # snv data
        self.snv_derived_cn_neutral = None
        self.snv_loh = None
        self.merged_candidates = None
        # output
        self.chromosome_names_out = None
        # DEV
        self.min_chr_length = 1e6
        self.dist_min_overwrite = 10000  # 10kb TODO
        self.max_std_outlier_rm = self.spectre_args.dev_max_std_outlier_rm  # 5
        self.mosdepth_cov_genome_chr_diff = self.spectre_args.dev_mosdepth_cov_genome_chr_diff  # 0.10  # 10%
        self.lower_2n_threshold = self.spectre_args.dev_lower_2n_threshold  # overwritten after data_normalization call
        self.upper_2n_threshold = self.spectre_args.dev_upper_2n_threshold  # overwritten after data_normalization call
        self.cov_diff_threshold = self.spectre_args.dev_cov_diff_threshold  # 0.80
        self.dist_proportion = self.spectre_args.dev_dist_proportion  # 0.25
        self.candidate_final_threshold = self.spectre_args.min_cnv_len  # min_cnv_len #100000  # 100kb
        self.sample_coverage_overwrite = self.spectre_args.sample_coverage_overwrite  # sample_coverage_overwrite
        # TODO
        self.output_directory = spectre_args.out_dir

    # Data normalization
    def data_normalization(self):
        """
        Normalize single chromosome bins
        """
        self.logger.debug(f'coverage file: {self.coverage_file}')
        if os.stat(self.coverage_file).st_size == 0:
            self.logger.error("Empty file")

        coverage_file_tabix = pysam.TabixFile(filename=self.coverage_file, index=f'{self.coverage_file}.csi')

        tmp_genome_wide_coverage_dict = {}
        tmp_genome_wide_position_dict = {}
        # load coverage data per chromosome
        for reference_chromosome in self.genome_info["chromosomes"]:
            # init
            tmp_positions = []
            tmp_coverage = []

            # load coverage data
            for tbx_line in coverage_file_tabix.fetch(reference_chromosome):
                [reference_chromosome, start, _, coverage] = tbx_line.split("\t")
                tmp_positions.append(int(start))
                tmp_coverage.append(float(coverage))

            self.coverage = np.array(tmp_coverage)
            self.positions = np.array(tmp_positions)

            # calculate bin size (if unknown)
            if self.bin_size == 0:
                self.bin_size = abs(self.positions[1] - self.positions[0])
                self.logger.info(f"Determined bin size from Mosdepth coverage: {self.bin_size}")
                # adjusting metadata according to bin size
                fm = FastaRef()
                self.metadata = fm.load_metadata(mdr_file_path=self.metadata_path,
                                                 blacklist_file_path=self.spectre_args.black_list,
                                                 bin_size=self.bin_size)
            # remove n regions from coverage
            self.__remove_n_region_by_chromosome(reference_chromosome)
            # add coverage to the genome-wide coverage dict
            tmp_genome_wide_coverage_dict[reference_chromosome] = self.coverage
            tmp_genome_wide_position_dict[reference_chromosome] = self.positions

        # Calculate the normalization statistics
        genome_wide_cov_values = [value for chrom_cov_values in tmp_genome_wide_coverage_dict.values() for value in
                                  chrom_cov_values]

        # Calculate the median of the values using numpy
        self.genome_wide_median = float(np.nanmedian(genome_wide_cov_values))
        self.genome_wide_mean = float(np.nanmean(genome_wide_cov_values))
        if self.sample_coverage_overwrite is not None:
            self.logger.info(f"Overwriting sample coverage with {self.sample_coverage_overwrite}")
            self.genome_wide_median = self.sample_coverage_overwrite
        self.normalization_value = self.genome_wide_median
        self.genome_std = np.nanstd(genome_wide_cov_values)

        # check if the normalization value is within the expected range the range is plus minus 1 genome_std from the
        # normalization value
        if self.snv_derived_cn_neutral is not None:
            if not (self.normalization_value - self.genome_std < self.snv_derived_cn_neutral["med"] < self.normalization_value + self.genome_std):
                self.logger.warning(f'The median coverage value between the Mosdepth data and the SNV data are different!')
                self.logger.warning(f'Coverage: Mosdepth: {self.normalization_value}X. SNV: {self.snv_derived_cn_neutral["med"]}X.')
                self.logger.warning(f'Expect trouble!')

        # Apply normalization to the coverage data
        for reference_chromosome in self.spectre_args.only_chr_list:
            self.coverage = tmp_genome_wide_coverage_dict[reference_chromosome]
            self.positions = tmp_genome_wide_position_dict[reference_chromosome]
            cov_stats, norm_stats, cov_data = self.__normalization_and_statistics(reference_chromosome)
            if not cov_stats.empty:
                self.coverage_analysis[reference_chromosome] = {"cov_data": cov_data, "statistics": cov_stats,
                                                                "norm_statistics": norm_stats}

        # Show which chromosomes did not match between the coverage file and the reference genome
        # for missing in set(coverage_file_tabix.contigs) - set(chrosome_matching_list):
        #    self.logger.warning(f"NO reference sequence found for: {missing}")
        coverage_file_tabix.close()

        # setup, calculating thresholds for deletion and duplication
        self.cnv_metrics = CNVMetrics(coverage_analysis=self.coverage_analysis,
                                      cnv_calls=self.cnv_calls_list,
                                      exclusion_zones=self.metadata,
                                      hashname=self.sample_id,
                                      ploidy=self.spectre_args.ploidy,
                                      output_dir=self.output_directory, is_cancer=self.spectre_args.is_cancer,
                                      as_dev=self.as_dev, debug_dir=self.debug_dir)
        # self.cnv_metrics.calculate_del_dup_borders(ploidy=self.spectre_args.ploidy)
        self.cnv_metrics.calculate_del_dup_borders_quantile(left_quantile=0.1, right_quantile=0.90)
        # self.cnv_metrics.calculate_del_dup_borders_sd()

        self.lower_2n_threshold = self.cnv_metrics.del_threshold
        self.upper_2n_threshold = self.cnv_metrics.dup_threshold
        self.logger.info(f'Selected threshold for DEL={np.round(self.lower_2n_threshold, 4)} and '
                         f'DUP={np.round(self.upper_2n_threshold, 4)} for a ploidy of {self.spectre_args.ploidy}')
        # clean up
        self.positions = np.zeros(0)
        self.coverage = np.zeros(0)
        self.coverage_log2 = np.zeros(0)
        pass

    def __normalization_and_statistics(self, chromosome) -> tuple:
        self.logger.info(f'Number positions to be tested on chromosome {chromosome}: {self.coverage.size}')
        # clean up
        cov_stats = CoverageStatistics()
        norm_stats = CoverageStatistics()
        cov_data = CoverageData()
        # genome_avg_cov = self.mosdepth_data.genome_mean_coverage  # TODO: do we still need that?
        if self.coverage.size > 0:
            # Data normalization by mean coverage
            # Determine if LOH (SNV) mode yes/no and normalize accordingly
            if self.snv_derived_cn_neutral is not None:
                avg, std, med = self.snv_derived_cn_neutral["avg"], self.snv_derived_cn_neutral["std"], \
                    self.snv_derived_cn_neutral["med"]
                self.normalization_value = med  # update normalization value if SNV derived CN neutral is available
            else:
                self.normalization_value = self.genome_wide_median

            # remove outliers (right tail) and re-calculate avg, std, median without outliers (right tail)
            if not self.spectre_args.disable_max_coverage:
                for idx in range(0, self.coverage.size, 1):
                    cov = self.coverage[idx]
                    if cov > self.genome_wide_mean + (self.max_std_outlier_rm * self.genome_std):
                        self.coverage[idx] = np.NaN
            else:
                self.logger.warning("Outlier max coverage is disabled, no outlier removal will be performed")
            nan_count = np.count_nonzero(np.isnan(self.coverage))
            self.logger.debug(f'nan_count {nan_count} v. {self.coverage.size} self.coverage.size')
            # show results in debug
            # compute statistics
            if chromosome in self.genome_info["chr_lengths"] and nan_count < self.coverage.size:
                [avg, std, med] = [np.nanmean(self.coverage), np.nanstd(self.coverage), np.nanmedian(self.coverage)]
                self.logger.debug(f'Chromsomal Info @ {chromosome}: \t avg: {avg}, std: {std}, med: {med}')

                # === Raw coverage data collection ===
                # data
                cov_data.positions = self.positions
                cov_data.coverage_raw = self.coverage
                cov_data.empty = False
                # stats
                cov_stats.chromosome_len = self.genome_info["chr_lengths"][chromosome]
                cov_stats.chromosome_name = chromosome
                cov_stats.average = avg
                cov_stats.std_dev = std
                cov_stats.median = med
                cov_stats.genome_wide_median = self.genome_wide_median
                cov_stats.genome_wide_mean = self.genome_wide_mean
                cov_stats.min = np.nanmin(self.coverage)
                cov_stats.max = np.nanmax(self.coverage)
                cov_stats.empty = False
                # === Normalization ===
                # normalization, based on diploid organisms, using median coverage
                normalized_candidates = self.coverage / self.normalization_value
                normalized_candidates_ploidy = normalized_candidates * self.genome_info["chr_ploidy"][chromosome]
                # statistics for normalized candidates
                [avg, std, med] = [np.nanmean(normalized_candidates_ploidy), np.nanstd(normalized_candidates_ploidy),
                                   np.nanmedian(normalized_candidates_ploidy)]
                norm_stats.chromosome_len = self.genome_info["chr_lengths"][chromosome]
                norm_stats.chromosome_name = chromosome
                norm_stats.median = med
                norm_stats.genome_wide_median = self.normalization_value  # changed depending if SNV or only coverage is used
                norm_stats.average = avg
                norm_stats.std_dev = std
                norm_stats.min = np.nanmin(normalized_candidates_ploidy)
                norm_stats.max = np.nanmax(normalized_candidates_ploidy)
                norm_stats.empty = False
                # norm data
                cov_data.normalized_cov = normalized_candidates
                cov_data.normalized_cov_ploidy = normalized_candidates_ploidy
            else:
                self.logger.warning(f"Chromosome {chromosome}: no data available!")
        return cov_stats, norm_stats, cov_data

    def __remove_n_region_by_chromosome(self, chrom) -> None:
        """
        Compare CNV candidates with "N" regions from the reference genome. If the candidate starts or ends inside a
        "N" region, the candidate will is removed by assigning a np.nan value.
        :param chrom: current chromosome
        :return: None
        """
        # we use two variables which contain all the current chromosome info
        # self.positions
        # self.coverage
        if chrom in self.metadata:
            for meta_start, meta_end in self.metadata[chrom]:
                meta_start = int(meta_start)
                meta_end = int(meta_end)
                overlap_n = np.logical_and(self.positions >= meta_start, self.positions <= meta_end)
                self.coverage[overlap_n] = np.nan

    # *** MAIN ***
    def call_cnv_coverage(self):
        # unique CNV IDs
        for each_chromosome in self.coverage_analysis.keys():
            # global
            # Section 1: Generating raw CNV candidates
            self.logger.info(f"Calculating CNVs for {self.sample_id} @ chromosome {each_chromosome}")
            cnv_caller = CNVCall(self.as_dev)
            candidates_cnv_list = cnv_caller.cnv_coverage(self.coverage_analysis[each_chromosome]["cov_data"],
                                                          self.bin_size, each_chromosome, self.sample_id,
                                                          self.cnv_metrics.del_threshold,
                                                          self.cnv_metrics.dup_threshold)

            # Generating CNV IDs
            new_ids = CNV_ID.n_id_generator(existing_ids=self.existing_cnv_ids, n=len(candidates_cnv_list))
            self.existing_cnv_ids = self.existing_cnv_ids + new_ids
            for candidate, candidate_id in zip(candidates_cnv_list, new_ids):
                candidate.set_id(candidate_id)

            # save candidates for each chromosome
            self.cnv_calls_list[each_chromosome] = candidates_cnv_list
        self.raw_cnv_calls_list = self.cnv_calls_list.copy()

    def refine_cnv_calls(self):
        self.logger.info("refining cnv calls")
        for each_chromosome in self.coverage_analysis.keys():
            candidates_cnv_list = self.cnv_calls_list[each_chromosome]
            if self.as_dev:
                self.dev_write_csv(each_chromosome)

            self.logger.debug(f'total cnv candidates in {each_chromosome}: {len(candidates_cnv_list)} before merge')
            # 1) Merging CNV candidates
            final_cnv_candidates = self.merge_candidates(candidates_cnv_list, each_chromosome)
            # 2) Scaffolding CNV candidates
            cnv_call_list = self.scaffold_candidates(final_cnv_candidates, each_chromosome) \
                if len(final_cnv_candidates) >= 2 else final_cnv_candidates
            # 3) Final CNV candidates CNV Post Processing
            # 3.1) Try extending and trimming CNVs
            cnv_call_list = self.trim_merged_candidates(cnv_call_list, 1.0)
            # 3.2) Check resulting CNVs from trim_merged_candidates and check if they are overlapping.
            cnv_call_list = self.clean_by_overlap(cnv_call_list)
            cnv_call_list = self.clean_merged_candidates(cnv_call_list, 1.0)
            self.cnv_calls_list[each_chromosome] = cnv_call_list
        pass

    def check_snfj_breakpoints(self, snfspc_file):
        bp = Breakpoints(snfspc_file, self.as_dev)
        bp.correlate_cnvs_with_breakpoints(cnv_candidates=self.cnv_calls_list, bin_size=self.bin_size)

    # Candidate CNV merge by distance and CNV type
    def merge_candidates(self, candidates_cnv_list, chromosome_name):
        merged_candidates = []
        if len(candidates_cnv_list) > 1:
            dev_candidates_string = ""
            merge_rounds = 1
            self.logger.debug(f'Current merge cycle {merge_rounds}')
            [n_merges, merged_candidates, _] = self.cnv_candidate_merge(candidates_cnv_list)
            while n_merges > 0:
                merge_rounds += 1
                self.logger.debug(f'Current merge cycle {merge_rounds}')
                # self.logger.info(f'n candidates in {chromosome_name}: {len(merged_candidates)}')
                [n_merges, merged_candidates, dev_candidates_string] = self.cnv_candidate_merge(merged_candidates)
            # self.logger.debug(dev_candidates_string)
            self.logger.debug(f'Total merge rounds: {merge_rounds}')
            self.logger.debug(f'Trying to merge over blacklisted regions for {chromosome_name}')
            merge_rounds = 1
            n_merges = 1
            while n_merges > 0:
                self.logger.debug(f'Current blacklist merge cycle {merge_rounds}')
                [n_merges, merged_candidates, _] = self.cnv_candidate_merge(cnv_candidates=merged_candidates,
                                                                            allow_over_blacklist_merging=True)
            self.logger.debug(f'Total merge rounds over blacklisted region: {merge_rounds}')

            # self.logger.debug(f'n candidates {chromosome_name}: {len(candidates_cnv_list)}')
        return merged_candidates if len(merged_candidates) > 0 else candidates_cnv_list

    def cnv_candidate_merge(self, cnv_candidates, allow_over_blacklist_merging: bool = False):
        def dev_candidates_merge(header=False):  # dev
            if header:
                return f'chr\tleft\tright\tleft_size\tright_size\tcov_left\tcov_right\tdist_OK\ttype_OK\tCN_OK'
            else:
                return (f'{cnv_candidates[idx].chromosome}\t{end0}\t{start1}\t'
                        f'{cnv_cand.size}\t{cnv_candidates[idx + 1].size}\t'
                        f'{cnv_candidates[idx].median_cov_norm}\t{cnv_candidates[idx + 1].median_cov_norm}\t'
                        f'{merge_distance_allowed}\t{same_type}\t{same_cnv_status}')

        n_candidates = len(cnv_candidates)
        # distance between end of i to start of i+1
        # first on is given
        idx = 0
        n_merges = 0
        merges_candidates = []
        cnv_cand = cnv_candidates[idx]
        cov_diff_threshold = self.cov_diff_threshold  # selected by trial and error in simulations
        # dist => distance
        dist_proportion = self.dist_proportion  # selected by trial and error, max distance proportion for merge
        dist_min_overwrite = self.dist_min_overwrite
        dev_candidates_string = [dev_candidates_merge(header=True)]
        # last candidate in the list was merged
        last_candidate_merged = False
        while idx < n_candidates - 1:
            end0 = cnv_candidates[idx].end
            start1 = cnv_candidates[idx + 1].start
            same_type = cnv_candidates[idx].type == cnv_candidates[idx + 1].type
            # copy number status difference
            same_cnv_status = abs(
                cnv_candidates[idx].median_cov_norm - cnv_candidates[idx + 1].median_cov_norm) <= cov_diff_threshold

            # distance between merging windows is a % of the average length of the windows to be merged
            # function to determine the merge distance which is primarily influenced by the size of the CNV.
            # Idea: The larger the CNV the closer the gap must be
            mean_cnv_distance = np.mean([cnv_cand.size, cnv_candidates[idx + 1].size])
            cnv_merge_threshold_percent = 40 * np.e ** (-mean_cnv_distance / 1000000) + 0.1  # dampening function
            cnv_merge_threshold = mean_cnv_distance * cnv_merge_threshold_percent / 100

            merge_distance = max(cnv_merge_threshold, dist_min_overwrite) + self.bin_size
            merge_distance_end = cnv_cand.end + merge_distance
            merge_distance_allowed = merge_distance_end >= cnv_candidates[idx + 1].start

            dev_candidates_string.append(dev_candidates_merge())
            # Evaluate if the new span is covered by any blacklisted region
            # Check if span_over_blacklisted_region = true if no overlap exists
            if cnv_cand.chromosome not in self.metadata.keys():
                # No blacklisted regions for this chromosome
                span_over_blacklisted_region = False
            else:
                # Found chromosome in blacklist and needs to check all positions
                span_over_blacklisted_region = any(
                    any(end0 <= int(val) <= start1 for val in tup) for tup in self.metadata[cnv_cand.chromosome])

            if allow_over_blacklist_merging:
                # self.logger.debug(f"MERGING OVER BLACKLISTED REGION {span_over_blacklisted_region}")
                span_over_blacklisted_region = False

            # check if
            is_outside_of_thresholds = self.is_cnv_cov_outside_of_thresholds(cnv_cand)

            if merge_distance_allowed and same_type and same_cnv_status and \
                    not span_over_blacklisted_region and is_outside_of_thresholds:
                # merge and extend
                # self.logger.debug(f"Merging {end0} and {start1} with distance {merge_distance}")
                n_merges += 1
                cnv_cand.add_candidates(
                    cnv_candidates[idx + 1].pos,
                    cnv_candidates[idx + 1].cov,
                    # cnv_candidates[idx + 1].raw_cov,
                    cnv_candidates[idx + 1].id,
                    cnv_candidates[idx + 1].merged_sample_references
                )
                merges_candidates.append(cnv_cand)
                idx += 1  # set idx to next candidate will be skipped

                # check the merged candidate was the last on in the list
                if idx >= n_candidates - 1:  # idx +1 (> not >= )is next candidate still in range of list? NO -> last candidate merged
                    last_candidate_merged = True

                # self.logger.debug(f"MERGED candidate to {cnv_cand.start}-{cnv_cand.end}")
            else:
                # save and show
                cnv_cand.median_coverage_candidates_merged()
                merges_candidates.append(cnv_cand)
                # cnv_cand.show_candidates()  # DEV

            # init next candidate.
            idx += 1
            if idx < n_candidates:
                cnv_cand = cnv_candidates[idx]

        # one last time for the remaining ones only if last merge was not successful
        if not last_candidate_merged:  # in case of a merge of the last-1 and last we have to add the extended candidate
            merges_candidates.append(cnv_cand)
        # cnv_cand.show_candidates()  # DEV
        return [n_merges, merges_candidates, "\n".join(dev_candidates_string)]

    def clean_merged_candidates(self, merged_candidates, size_multiplier: float = 1.0):
        # current threshold is 100kb
        clean_merged = []
        pre_clean_cnt = len(merged_candidates)
        clean_merged = self.clean_by_threshold(merged_candidates)
        clean_merged = self.clean_by_size(clean_merged, size_multiplier)
        # clean_merged = self.clean_by_overlap(clean_merged)
        # clean_merged = self.clean_by_coverage(clean_merged, size_multiplier)
        post_clean_cnt = len(clean_merged)
        self.logger.debug(
            f"Clean by size: {pre_clean_cnt} -> {post_clean_cnt}. Cleaned: {pre_clean_cnt - post_clean_cnt}")
        return clean_merged

    def trim_merged_candidates(self, merged_candidates, size_multiplier: float = 1.0):
        # Extend and trim candidates
        if len(merged_candidates) > 0:
            chromosome_name = merged_candidates[0].chromosome
            cpp = CNV_Post_Processing(
                ploidy_coverage=self.coverage_analysis[chromosome_name]["cov_data"].normalized_cov_ploidy,
                position=self.coverage_analysis[chromosome_name]["cov_data"].positions,
                bin_size=self.bin_size,
                as_dev=self.as_dev)
            # Try extending and trimming the candidates
            adjusted_candidate_cnv_list = cpp.cnv_fit(merged_candidates)
            return adjusted_candidate_cnv_list
        # Fallback if empty return the original
        return merged_candidates

    def clean_by_overlap(self, merged_candidates):
        """
        Merge CNVs together the first CNV ends and take all values from the second CNV where the end from the first CNV
        is located at. Thus ensuring no data duplication.
        :param merged_candidates:
        :return:
        """
        clean_merged = []
        remove_idx = []
        number_of_candidates = len(merged_candidates)
        offset = 0

        # Walk through candidates by index (idx)
        for idx in range(len(merged_candidates) - 1):
            current_cand = merged_candidates[idx]
            early_stop = False

            # Walk through all following candidates (offset)
            while idx + offset < number_of_candidates - 1 and not early_stop:
                next_cand = merged_candidates[idx + offset + 1]

                # Break if the current and next candidate are not overlapping at any point
                if not (next_cand.start < current_cand.start < next_cand.end or
                        next_cand.start < current_cand.end < next_cand.end):
                    early_stop = True
                    break

                # Check of full overlap
                current_cand_is_fully_overlapping = self.__is_candidate_fully_overlaping(current_cand, next_cand)
                next_cand_is_fully_overlapping = self.__is_candidate_fully_overlaping(next_cand, current_cand)
                if current_cand_is_fully_overlapping:
                    # remove the next candidate
                    remove_idx.append(idx + offset + 1)
                    offset += 1
                if next_cand_is_fully_overlapping:
                    # remove the current candidate
                    remove_idx.append(idx)
                    early_stop = True
                else:
                    # Partial overlap
                    # # Determine which candidate is larger and use it for the value extension as a base
                    master_cand = current_cand
                    slave_cand = next_cand
                    # remove  slave candidate by index
                    remove_candidate = idx + offset + 1
                    if current_cand.size < next_cand.size:
                        # Extend current candidate by the values of the next candidate.
                        master_cand = next_cand
                        slave_cand = current_cand
                        # remove current candidate by index
                        remove_candidate = idx

                    # # Determine where the smaller candidate overlaps (left or right side).
                    is_left_overlap = True
                    if master_cand.end < slave_cand.end:
                        is_left_overlap = False

                    # # Determine with start and end position which values should be use for the extension. (No value
                    # duplication)
                    if is_left_overlap:
                        # get values from the start of the slave candidate to the start of the master candidate
                        # determine the overlapping index
                        overlap_idx = slave_cand.pos.index(master_cand.start)  #
                        # extend the values of the master candidate by the values of the slave candidate
                        master_cand.add_candidates(cnv_cand_pos=slave_cand.pos[:overlap_idx],
                                                   cnv_cand_cov=slave_cand.cov[:overlap_idx],
                                                   cnv_id=slave_cand.id,
                                                   cnv_merged_ids=slave_cand.merged_sample_references,
                                                   push_front=True)
                    else:
                        # get values from the end of the master candidate to the end of the slave candidate
                        # determine the overlapping index
                        overlap_idx = slave_cand.pos.index(master_cand.end) + 1  # +1 to exclude the overlapping value
                        # extend the values of the master candidate by the values of the slave candidate
                        master_cand.add_candidates(cnv_cand_pos=slave_cand.pos[overlap_idx:],
                                                   cnv_cand_cov=slave_cand.cov[overlap_idx:],
                                                   cnv_id=slave_cand.id,
                                                   cnv_merged_ids=slave_cand.merged_sample_references,
                                                   push_front=False)
                    # # Remove the slave candidate from the list.
                    remove_idx.append(remove_candidate)
                    offset += 1

        # remove the candidates that overlap
        for idx in range(len(merged_candidates)):
            if idx not in remove_idx:
                clean_merged.append(merged_candidates[idx])
        return clean_merged

    def __is_candidate_fully_overlaping(self, candidate1, candidate2):
        """
        Check if one candidate is overlapping the whole other one
        :param candidate1:
        :param candidate2:
        :return:
        """
        return candidate1.start <= candidate2.start and candidate2.end <= candidate1.end

    def clean_by_size(self, merged_candidates, size_multiplier: float = 1.0):
        clean_merged = []
        pre_clean_cnt = len(merged_candidates)
        for each_cand in merged_candidates:
            if each_cand.size >= self.candidate_final_threshold * size_multiplier:
                clean_merged.append(each_cand)
        post_clean_cnt = len(clean_merged)
        self.logger.debug(
            f"Clean by size: {pre_clean_cnt} -> {post_clean_cnt}. Cleaned: {pre_clean_cnt - post_clean_cnt}")
        return clean_merged

    def clean_by_coverage(self, merged_candidates, size_multiplier: float = 1.0):
        clean_merged = []
        for each_cand in merged_candidates:
            if self.is_cnv_cov_outside_of_thresholds(each_cand):
                clean_merged.append(each_cand)
        return clean_merged

    def clean_by_threshold(self, merged_candidates):
        # check if the
        clean_merged = []
        pre_clean_cnt = len(merged_candidates)
        del_threshold, dup_threshold = self.cnv_metrics.calculate_del_dup_borders_sd(sd_clip_n=1.85)
        for each_cand in merged_candidates:
            threshold = np.mean(each_cand.cov)

            # count the number of values that are outside of del and dup thresholds
            del_outside = np.count_nonzero(np.array(each_cand.cov) < del_threshold)
            dup_outside = np.count_nonzero(np.array(each_cand.cov) > dup_threshold)
            n_outside_ci = del_outside + dup_outside
            # check if more than 50% of the values are outside the confidence interval
            has_more_values_outside_of_normal_cov_space = n_outside_ci > len(
                each_cand.cov) * 0.5  # reduce FP due to noise

            # if candidate coverage mean is outside the thresholds, keep it
            if not (del_threshold <= threshold <= dup_threshold) and has_more_values_outside_of_normal_cov_space:
                clean_merged.append(each_cand)
            else:
                pass
        post_clean_cnt = len(clean_merged)

        self.logger.debug(
            f"Clean by threshold: {pre_clean_cnt} -> {post_clean_cnt}. Cleaned: {pre_clean_cnt - post_clean_cnt}")
        return clean_merged

    def is_cnv_cov_outside_of_thresholds(self, candidate):

        result = not (self.cnv_metrics.del_threshold <= candidate.median_cov_norm <= self.cnv_metrics.dup_threshold)
        # self.logger.debug(
        #    f"Clean by coverage: {result} -> {candidate.chromosome}:{candidate.start}-{candidate.end} = c_norm_mean: {c_norm_mean}, c_norm_median: {c_norm_median}")
        return result

    # Fill in the gaps after merge
    def scaffold_candidates(self, final_candidates_list, chromosome_name):
        cnv_list = []
        n_scaffolds = 0
        cov_data_chr = self.coverage_analysis[chromosome_name]["cov_data"]
        _index = 0
        while _index < len(final_candidates_list):
            _cand = final_candidates_list[_index]
            # get the difference between each position
            position_diff = np.diff(_cand.pos)
            # get index of the position_diff where the difference is not equal to the bin size
            position_diff_index = np.where(position_diff != self.bin_size)
            if len(position_diff_index[0]) > 0:
                # fill gaps for each position_diff_index and  reversing list with [::-1] to ensure no index overwritten
                # for idx1, idx2 in zip(position_diff_index[0][::-1][:-1],position_diff_index[0][::-1][1:]):
                for idx1 in position_diff_index[0][::-1]:

                    # Indices of current  cov gap in the current candidate
                    cnv_gap_start_idx = idx1
                    cnv_gap_end_idx = idx1 + 1

                    # Get chromosomal positions of the gap in the current candidate
                    gap_in_candidate_start = int(_cand.pos[cnv_gap_start_idx])
                    gap_in_candidate_end = int(_cand.pos[cnv_gap_end_idx])

                    # get gap indices in the coverage data
                    gap_start = int(np.where(cov_data_chr.positions == gap_in_candidate_start)[0][0])
                    gap_end = int(np.where(cov_data_chr.positions == gap_in_candidate_end)[0][0])

                    # get gap data
                    gap_data = cov_data_chr.normalized_cov_ploidy[gap_start:gap_end]
                    gap_pos = cov_data_chr.positions[gap_start:gap_end]
                    scaf_cov = np.nanmedian(gap_data)

                    # get array with scaf_cov values the length of the gap
                    scaf_cov_array = np.full(len(gap_data), scaf_cov)

                    # Restricting cnv scaffolding merge over 100 following N regions (NaN values)
                    scaf_cov_nan_overflow = np.count_nonzero(np.isnan(np.array(gap_data))) > 0

                    # Apply coverage to the gap if not to many NaN values are present
                    if not scaf_cov_nan_overflow:
                        _cand.add_scaffold_candidate(cnv_gap_start_idx, cnv_gap_end_idx, scaf_cov_array.tolist(),
                                                     gap_pos.tolist())
                        n_scaffolds += 1

            _cand.median_coverage_candidates_merged()
            cnv_list.append(_cand)
            _index += 1

        return cnv_list

    def cnv_loh_candidate_merge(self):
        # autorun
        # if self.snv_loh:
        if len(self.coverage_analysis) == 0:
            self.chromosome_names_out = self.spectre_args.only_chr_list
        else:
            self.chromosome_names_out = self.coverage_analysis.keys()
        merge_cnv_loh = MergeCNVLoH(chr_list=self.chromosome_names_out, cnv=self.cnv_calls_list, loh=self.snv_loh)
        self.logger.info(f'No error occured while merging CNV and LoH candidates') if merge_cnv_loh.error else None
        if not merge_cnv_loh.error:
            self.merged_candidates = merge_cnv_loh.merged_candidates
        else:
            self.logger.warning("An error occurred during the merge, only CNV will be output")
            self.merged_candidates = self.cnv_calls_list

    # Output results
    def cnv_result_bed(self, method=""):
        if method == "":
            output_bed = self.output_bed
        else:
            output_bed = os.path.join(os.path.join(self.output_directory, f'{method}{self.sample_id}.bed.gz'))

        bed_output = outputWriter.BedOutput(output_bed)
        bed_output.make_bed(self.chromosome_names_out, self.merged_candidates)

        if self.as_dev:
            bed_output_dev = outputWriter.BedOutput(f'{self.debug_dir}/raw-cnvs{self.sample_id}.bed.gz')
            bed_output_dev.make_bed(self.chromosome_names_out, self.raw_cnv_calls_list)

    def cnv_result_vcf(self, method=""):
        # TODO add LoH if self.snv_loh is not None
        # TODO modify make_vcf
        output_vcf = os.path.join(os.path.join(self.output_directory, f'{method}{self.sample_id}.vcf.gz'))
        vcf_output = outputWriter.VCFOutput(output_vcf, self.genome_info)
        vcf_output.make_vcf(self.chromosome_names_out, self.merged_candidates, self.sample_id)

    # Plots
    def coverage_plot(self):
        for each_chromosome in self.coverage_analysis.keys():
            plot_cnv = CoveragePlot()
            plot_cnv.file_prefix = f'{self.sample_id}_coverage_plot'
            plot_cnv.output_directory = self.output_directory
            chr_cov_data = self.coverage_analysis[each_chromosome]["cov_data"]
            plot_cnv.plot_coverage(each_chromosome, {"pos": chr_cov_data.positions, "cov": chr_cov_data.raw_coverage})

    def cnv_plot(self, methode=""):
        for each_chromosome in self.coverage_analysis.keys():
            chr_cov_data = self.coverage_analysis[each_chromosome]["cov_data"]
            cov_stats = self.coverage_analysis[each_chromosome]["norm_statistics"]
            lower_bound = self.lower_2n_threshold
            upper_bound = self.upper_2n_threshold
            cov_stats.median = cov_stats.median * self.spectre_args.ploidy  # for diploid
            new_plot_device = CNVPlot()
            new_plot_device.output_directory = self.output_directory
            new_plot_device.file_prefix = methode + self.sample_id
            new_plot_device.plot_coverage_cnv(each_chromosome, cov_stats,
                                              {"pos": chr_cov_data.positions,
                                               "cov": chr_cov_data.normalized_cov_ploidy},
                                              self.cnv_calls_list[each_chromosome], [lower_bound, upper_bound])

    def get_cnv_metrics(self, refined_cnvs: bool = False):
        """
        Calculates for every CNV call a metrics like the p value with statistical tests
        :return:
        """
        self.cnv_metrics.evaluate_cnvs(refined_cnvs)

    def write_intermediate_candidates(self) -> str:
        """
        Writing intermediate object files out
        :param candidate_type:
        :return:
        """
        loh_cnv_calls_list_dict = {}
        loh_raw_cnv_calls_list_dict = {}

        intermediate_output_writer = outputWriter.IntermediateFile(self.output_directory)
        # genome_info = intermediate_output_writer.convert_candidates_to_dictionary(self.genome_info)
        cnv_calls_list_dict = intermediate_output_writer.convert_candidates_to_dictionary(self.cnv_calls_list)
        raw_cnv_calls_list_dict = intermediate_output_writer.convert_candidates_to_dictionary(self.raw_cnv_calls_list)

        # merge loh dict to  cnv calls dict
        if self.cnv_calls_list_af_filtered is not None:
            loh_raw_cnv_calls_list_dict = intermediate_output_writer.convert_candidates_to_dictionary(self.snv_loh)
            loh_cnv_calls_list_dict = intermediate_output_writer.convert_candidates_to_dictionary(
                self.cnv_calls_list_af_filtered)

        analysis_dict = {
            "metadata": {"source": "spectre", "spectre_version": "0.2.alpha", "bin_size": int(self.bin_size),
                         "min_cn_len": int(self.candidate_final_threshold),
                         "normalization_value": float(self.normalization_value), "sample_id": self.sample_id},
            "genome_info": self.genome_info,
            "raw_cnvs": raw_cnv_calls_list_dict,
            "refined_cnvs": cnv_calls_list_dict,
            "refined_loh_cnvs": loh_cnv_calls_list_dict,
            "raw_loh_cnvs": loh_raw_cnv_calls_list_dict,
            "analysis_metrics": {
                "min_chr_length": int(self.min_chr_length),
                "max_std_outlier_rm": int(self.max_std_outlier_rm),
                "mosdepth_cov_genome_chr_diff": float(self.mosdepth_cov_genome_chr_diff),
                "lower_2n_threshold": float(self.lower_2n_threshold),
                "upper_2n_threshold": float(self.upper_2n_threshold),
                "strict_del_threshold": float(self.cnv_metrics.strict_del_threshold),
                "strict_dup_threshold": float(self.cnv_metrics.strict_dup_threshold),
                "cov_diff_threshold": float(self.cov_diff_threshold),
                "dist_proportion": float(self.dist_proportion),
                "dist_min_overwrite": int(self.dist_min_overwrite),
                "candidate_final_threshold": int(self.candidate_final_threshold),
                "genome_mean": float(
                    np.mean(self.cnv_metrics.df_coverage_candidate_no_excl_zone_random_samples['coverage'])),
                "genome_sd": float(
                    np.std(self.cnv_metrics.df_coverage_candidate_no_excl_zone_random_samples['coverage'])),
                "genome_var": float(
                    np.var(self.cnv_metrics.df_coverage_candidate_no_excl_zone_random_samples['coverage']))
            }
        }
        output_path = intermediate_output_writer.write_intermediate_file(analysis_dict, f"{self.sample_id}")
        self.intermediate_candidates_file_location = output_path
        return output_path

    # ############################################
    # dev
    def dev_write_csv(self, each_chromosome):
        csv_results = pd.DataFrame(
            data={"position": self.coverage_analysis[each_chromosome]["cov_data"].positions,
                  "mosdepth_cov": self.coverage_analysis[each_chromosome]["cov_data"].coverage_raw,
                  "norm_cov": self.coverage_analysis[each_chromosome]["cov_data"].normalized_cov,
                  "ploidy_cov": self.coverage_analysis[each_chromosome]["cov_data"].normalized_cov_ploidy
                  }
        )
        csv_results["chr"] = each_chromosome
        output_file = f"{self.debug_dir}/cnv_{self.sample_id}_byCoverage_chr{each_chromosome}.csv"
        self.logger.debug(f"Writing coverage to {output_file}")
        csv_results.to_csv(output_file, index=False)
