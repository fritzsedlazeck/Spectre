import os
import gzip
import pysam
import numpy as np
import pandas as pd
import logging as logger
import util.outputWriter
import util.vcf_parser as vcf
from util.OSUtil import OSUtil as OSUt
from util.dataAnalyzer import NormaldataAnalyser as NorAn
from util.cnv_id import CNV_ID
from plots.plot import CNVPlot
from plots.plot import CoveragePlot
from analysis.coverage_stats import CoverageStatistics
from analysis.coverage_stats import CoverageData
from analysis.call_cnv_coverage import CNVCall
from analysis.call_cnv_AF import CNVCall as CNVAnalysisSNP
from analysis.cnv_metrics import CNVMetrics


class CNVAnalysis(object):
    def __init__(self, coverage_file, coverage_mosdepth_data, bin_size, output_directory, outbed, outvcf, genome_info,
                 sample_id, metadata_ref, snv_file, only_chr_list="", ploidy=2,min_cnv_len=1000000, as_dev=False,
                 dev_params=None, debug_dir=""):
        # input
        self.coverage_file = coverage_file
        self.mosdepth_data = coverage_mosdepth_data  # has chr_mean_coverage, genome_mean_coverage, genome_bases
        self.output_directory = output_directory
        self.output_bed = outbed
        self.output_vcf = outvcf
        self.sample_id = sample_id
        self.snv_file = snv_file
        self.population_sample_ids = set()  # hold the filename of the population samples
        # bin size !important
        self.bin_size = int(bin_size)
        self.genome_info = genome_info
        self.metadata = metadata_ref
        # whitelist, only this chr
        self.only_chr_list = str(only_chr_list).split(",") if only_chr_list != "" else genome_info["chromosomes"]
        self.ploidy = ploidy
        # logger
        logger.basicConfig(level=logger.DEBUG) if as_dev else logger.basicConfig(level=logger.INFO)
        self.logger = logger
        # work variables
        self.positions = np.zeros(0)
        self.coverage = np.zeros(0)
        self.coverage_log2 = np.zeros(0)
        # results by chromosome
        self.genome_analysis = {}  # use this for CNV call NOTE: here is the coverage data under "cov_data"
        self.cnv_calls_list = {}  # use this for CNV call
        self.raw_cnv_calls_list = {}  # use this for storing raw CNV call
        self.intermediate_candidates_file_location = ""  # holds output path of serialized cnv object
        self.cnv_merged = {}  # use this for CNV call
        self.existing_cnv_ids = []  # holds all already used cnv IDs
        # snv data
        self.snv_af_bed = ""
        self.snv_af_df = None
        # dev/debug + hidden params
        self.as_dev = as_dev
        self.debug_dir = debug_dir
        # cnv metrics
        self.cnv_metrics = None

        # TODO
        self.min_chr_length = 1e6
        self.max_std_outlier_rm = 5
        self.mosdepth_cov_genome_chr_diff = 0.10  # 10%
        self.lower_2n_threshold = 0.0  # are overwritten after data_normalization was called
        self.upper_2n_threshold = 0.0  # are overwritten after data_normalization was called
        self.cov_diff_threshold = 0.80
        self.dist_proportion = 0.25
        self.dist_min_overwrite = 10000  # 10kb
        self.candidate_final_threshold = min_cnv_len #100000  # 100kb

    # Data normalization
    def data_normalization(self):
        """
        Normalize single chromosome bins
        """
        self.logger.debug(self.coverage_file)
        file_size_lines = OSUt.get_lines_by_chromosome(self.coverage_file)
        if len(file_size_lines) == 0:
            self.logger.error("Empty file")

        coverage_full_path = os.path.abspath(os.path.expanduser(self.coverage_file))
        coverage_file_handler = gzip.open(coverage_full_path, "rt") if "gz" in self.coverage_file \
            else open(coverage_full_path, "r")

        list_index = 0
        previous_chromosome = ""
        for line in coverage_file_handler:
            [chromosome, start, _, coverage] = line.rstrip("\n").split("\t")
            if chromosome in self.only_chr_list:
                start = int(start)
                coverage = float(coverage)
                if previous_chromosome == "":
                    # init
                    self.positions = np.zeros(file_size_lines[chromosome])
                    self.coverage = np.zeros(file_size_lines[chromosome])
                    self.coverage_log2 = np.zeros(file_size_lines[chromosome])
                    list_index = 0
                    previous_chromosome = chromosome
                    # first elem
                    self.positions[list_index] = start
                    self.coverage[list_index] = coverage
                    self.coverage_log2[list_index] = np.log2(coverage) if coverage != 0 else np.NaN  # use np.nanFUN
                    list_index += 1
                elif previous_chromosome != chromosome:
                    # analysis here
                    self.logger.debug(previous_chromosome)
                    self.__remove_n_region_by_chromosome(previous_chromosome)
                    cov_stats, norm_stats, cov_data = self.__normalization_and_statistics(previous_chromosome)
                    self.genome_analysis[previous_chromosome] = {"cov_data": cov_data,
                                                                 "statistics": cov_stats,
                                                                 "norm_statistics": norm_stats}
                    # init new chromosome
                    self.positions = np.zeros(file_size_lines[chromosome])
                    self.coverage = np.zeros(file_size_lines[chromosome])
                    self.coverage_log2 = np.zeros(file_size_lines[chromosome])
                    list_index = 0
                    previous_chromosome = chromosome
                    # first elem
                    self.positions[list_index] = start
                    self.coverage[list_index] = coverage
                    self.coverage_log2[list_index] = np.log2(coverage) if coverage != 0 else np.NaN  # use np.nanFUN
                    list_index += 1
                else:
                    self.positions[list_index] = start
                    self.coverage[list_index] = coverage
                    self.coverage_log2[list_index] = np.log2(coverage) if coverage != 0 else np.NaN  # use np.nanFUN
                    list_index += 1
        coverage_file_handler.close()

        # compute leftover chromosome
        self.__remove_n_region_by_chromosome(previous_chromosome)
        cov_stats, norm_stats, cov_data = self.__normalization_and_statistics(previous_chromosome)
        self.genome_analysis[previous_chromosome] = {"cov_data": cov_data,
                                                     "statistics": cov_stats,
                                                     "norm_statistics": norm_stats}

        # setup, calculating border for deletion and duplications
        self.cnv_metrics = CNVMetrics(genome_analysis=self.genome_analysis,
                                      cnv_calls=self.cnv_calls_list,
                                      exclusion_zones=self.metadata,
                                      hashname=self.sample_id,
                                      ploidy=self.ploidy,
                                      output_dir=self.output_directory,
                                      as_dev=self.as_dev, debug_dir=self.debug_dir)
        # self.cnv_metrics.get_del_dup_borders()
        self.lower_2n_threshold = self.cnv_metrics.del_border
        self.upper_2n_threshold = self.cnv_metrics.dup_border

        # clean up
        self.positions = np.zeros(0)
        self.coverage = np.zeros(0)
        self.coverage_log2 = np.zeros(0)

    def __normalization_and_statistics(self, chromosome) -> tuple:
        self.logger.info(f'Number positions to be tested on chromosome {chromosome}: {self.coverage.size}')
        # clean up
        cov_stats = CoverageStatistics()
        norm_stats = CoverageStatistics()
        cov_data = CoverageData()
        cov_data.normalized_cov = np.array(0)
        if self.coverage.size > 0:
            # calculate avg, std, median and remove outliers (right tail)
            # use genome-wide average
            # NOTE: select average coverage from genome-wide or chromosome depending on how different they are
            genome_avg_cov = self.mosdepth_data.genome_mean_coverage
            chromosome_avg_coverage = np.nanmean(self.coverage)
            diff_genome_chr_coverage = abs(
                genome_avg_cov - chromosome_avg_coverage) < genome_avg_cov * self.mosdepth_cov_genome_chr_diff
            use_this_avg_cov = genome_avg_cov if diff_genome_chr_coverage else chromosome_avg_coverage
            [avg, std] = [float(use_this_avg_cov), np.nanstd(self.coverage)]
            for idx in range(0, self.coverage.size, 1):
                cov = self.coverage[idx]
                if cov > avg + (self.max_std_outlier_rm * std):
                    self.coverage[idx] = np.NaN

            # re-calculate avg, std, median without outliers (right tail)
            chromosome_avg_coverage = np.nanmean(self.coverage)
            use_this_avg_cov = genome_avg_cov  # if diff_genome_chr_coverage else chromosome_avg_coverage
            [avg, std, med] = [float(use_this_avg_cov), np.nanstd(self.coverage), np.nanmedian(self.coverage)]
            self.logger.debug([f'avg: {genome_avg_cov}|{chromosome_avg_coverage}, std: {std}, med: {med}'])

            if chromosome in self.genome_info["chr_lengths_by_name"]:
                # data
                cov_data.positions = self.positions
                cov_data.coverage_raw = self.coverage
                cov_data.coverage_log2 = self.coverage_log2
                # stats
                cov_stats.chromosome_len = self.genome_info["chr_lengths_by_name"][chromosome]
                cov_stats.chromosome_name = chromosome
                cov_stats.average = avg
                cov_stats.std_dev = std
                cov_stats.median = med
                cov_stats.min = np.nanmin(self.coverage)
                cov_stats.max = np.nanmax(self.coverage)

                # normalization, based on diploid organisms
                normalize_by = med  # med:median | avg:average
                normalized_candidates = NorAn.normalize_candidates(self.coverage, normalize_by)
                normalized_candidates_ploidy = normalized_candidates * self.ploidy
                if len(normalized_candidates) < 1:
                    normalized_candidates_ploidy = normalized_candidates

                # statistics for normalized bin
                avg, std, min_val, max_val, med = NorAn.get_candidate_statistics(normalized_candidates)
                norm_stats.chromosome_len = self.genome_info["chr_lengths_by_name"][chromosome]
                norm_stats.chromosome_name = chromosome
                norm_stats.median = med
                norm_stats.average = avg
                norm_stats.std_dev = std
                norm_stats.min = min_val
                norm_stats.max = max_val

                # norm data
                cov_data.normalized_cov = normalized_candidates
                cov_data.normalized_cov_ploidy = normalized_candidates_ploidy
            else:
                self.logger.warning(f"Chromosome:{chromosome} not found in gene info!")
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
    def call_cnv_coverage(self, cnv_id_list=None, write_csv=False):
        if cnv_id_list is None:
            cnv_id_list = []
        # unique CNV IDs
        for each_chromosome in self.genome_analysis.keys():
            # global
            # Section 1: Generating raw CNV candidates
            self.logger.info(f"Calculating CNVs for {self.sample_id} @ chromosome {each_chromosome}")
            cnv_caller = CNVCall(self.as_dev)
            candidates_cnv_list = cnv_caller.cnv_coverage(self.genome_analysis[each_chromosome]["cov_data"],
                                                          self.bin_size, each_chromosome, self.sample_id,
                                                          self.cnv_metrics.del_border, self.cnv_metrics.dup_border)

            # Generating CNV IDs
            new_ids = CNV_ID.n_id_generator(len(candidates_cnv_list), cnv_id_list)
            self.existing_cnv_ids = self.existing_cnv_ids + new_ids
            for candidate, id in zip(candidates_cnv_list, new_ids):
                candidate.set_id(id)

            # save candidates for each chromosome
            self.cnv_calls_list[each_chromosome] = candidates_cnv_list
        self.raw_cnv_calls_list = self.cnv_calls_list.copy()

    def refine_cnv_calls(self, write_csv=False):
        self.logger.info("refining cnv calls")
        for each_chromosome in self.genome_analysis.keys():
            candidates_cnv_list = self.cnv_calls_list[each_chromosome]
            if write_csv:
                self.dev_write_csv(each_chromosome)

            self.logger.debug(f'total cnv candidates in {each_chromosome}: {len(candidates_cnv_list)} before merge')
            final_cnv_candidates = self.merge_candidates(candidates_cnv_list, each_chromosome)
            cnv_call_list = self.scaffold_candidates(final_cnv_candidates, each_chromosome) \
                if len(final_cnv_candidates) >= 2 else final_cnv_candidates
            self.cnv_calls_list[each_chromosome] = cnv_call_list
        pass

    # Candidate CNV merge by distance and CNV type
    def merge_candidates(self, candidates_cnv_list, chromosome_name):
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
            candidates_cnv_list = self.clean_merged_candidates(merged_candidates)
            # self.logger.debug(f'n candidates {chromosome_name}: {len(candidates_cnv_list)}')
        return candidates_cnv_list

    def cnv_candidate_merge(self, cnv_candidates):
        def dev_candidates_merge(header=False):  # dev
            if header:
                return f'chr\tleft\tright\tleft_size\tright_size\tcov_left\tcov_right\tdist_OK\ttype_OK\tCN_OK'
            else:
                return (f'{cnv_candidates[idx].chromosome}\t{end0}\t{start1}\t'
                        f'{cnv_cand.size}\t{cnv_candidates[idx + 1].size}\t'
                        f'{cnv_candidates[idx].median_cov_norm}\t{cnv_candidates[idx + 1].median_cov_norm}\t'
                        f'{(start1 - end0) < max_merge_distance}\t{same_type}\t{same_cnv_status}')

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
        while idx < n_candidates - 1:
            end0 = cnv_candidates[idx].end
            start1 = cnv_candidates[idx + 1].start
            same_type = cnv_candidates[idx].type == cnv_candidates[idx + 1].type
            # copy number status difference
            same_cnv_status = abs(
                cnv_candidates[idx].median_cov_norm - cnv_candidates[idx + 1].median_cov_norm) <= cov_diff_threshold
            # distance between merging windows is a % of the average length of the windows to be merged
            max_merge_distance = max(np.mean([cnv_cand.size, cnv_candidates[idx + 1].size]) * dist_proportion,
                                     dist_min_overwrite)
            dev_candidates_string.append(dev_candidates_merge())
            # Evaluate if the new span is covered by any blacklisted region
            # Check if span_over_blacklisted_region = true if no overlap exists
            span_over_blacklisted_region = any(
                any(end0 <= int(val) <= start1 for val in tup) for tup in self.metadata[cnv_cand.chromosome])

            if (start1 - end0) < max_merge_distance and same_type and same_cnv_status and \
                    not span_over_blacklisted_region:
                # merge and extend
                n_merges += 1
                cnv_cand.add_candidates(
                    cnv_candidates[idx + 1].pos,
                    cnv_candidates[idx + 1].cov,
                    cnv_candidates[idx + 1].id,
                    cnv_candidates[idx + 1].merged_sample_references
                )
            else:
                # save and show
                cnv_cand.median_coverage_candidates_merged()
                merges_candidates.append(cnv_cand)
                # cnv_cand.show_candidates()  # DEV
                # init
                cnv_cand = cnv_candidates[idx + 1]
            idx += 1
        # one last time for the remaining ones
        merges_candidates.append(cnv_cand)
        # cnv_cand.show_candidates()  # DEV
        return [n_merges, merges_candidates, "\n".join(dev_candidates_string)]

    def clean_merged_candidates(self, merged_candidates):
        # current threshold is 100kb
        clean_merged = []
        for each_cand in merged_candidates:
            if each_cand.size >= self.candidate_final_threshold:
                clean_merged.append(each_cand)
        return clean_merged

    # Fill in the gaps after merge
    def scaffold_candidates(self, final_candidates_list, chromosome_name):
        cnv_list = []
        n_scaffolds = 0
        cov_data_chr = self.genome_analysis[chromosome_name]["cov_data"]
        _index = 0
        _cand = final_candidates_list[_index]
        while _index < len(final_candidates_list) - 1:
            _cand_next = final_candidates_list[_index + 1]
            if _cand.type == _cand_next.type:
                gap_start = np.where(cov_data_chr.positions == _cand.end)
                gap_end = np.where(cov_data_chr.positions == _cand_next.start)
                gap_data = cov_data_chr.normalized_cov_ploidy[gap_start[0][0]:gap_end[0][0]]
                scaf_cov = np.nanmedian(gap_data)

                # Restricting cnv scaffolding merge over 100 following N regions (NaN values)
                scaf_cov_nan_overflow = np.count_nonzero(np.isnan(np.array(gap_data))) > 0
                if abs(
                        scaf_cov - np.nanmedian(list(_cand.cov) + list(_cand_next.cov))
                ) <= self.cov_diff_threshold and not scaf_cov_nan_overflow:
                    _cand.add_candidates(_cand_next.pos, _cand_next.cov, _cand_next.id,
                                         _cand_next.merged_sample_references)
                    n_scaffolds += 1
                else:
                    _cand.median_coverage_candidates_merged()
                    cnv_list.append(_cand)
                    _cand = _cand_next
            else:
                _cand.median_coverage_candidates_merged()
                cnv_list.append(_cand)
                # init with next
                _cand = _cand_next
            _index += 1
        # one final
        _cand.median_coverage_candidates_merged()
        cnv_list.append(_cand)
        return cnv_list

    # Output results
    def cnv_result_bed(self, method=""):
        if method == "":
            output_bed = self.output_bed
        else:
            output_bed = os.path.join(os.path.join(self.output_directory, f'{method}{self.sample_id}.bed'))

        bed_output = util.outputWriter.BedOutput(output_bed)
        bed_output.make_bed(self.genome_analysis.keys(), self.cnv_calls_list)

    def cnv_result_vcf(self, method=""):
        output_vcf = os.path.join(os.path.join(self.output_directory, f'{method}{self.sample_id}.vcf'))
        vcf_output = util.outputWriter.VCFOutput(output_vcf, self.genome_info)
        vcf_output.make_vcf(self.genome_analysis.keys(), self.cnv_calls_list, self.sample_id)

    # Plots
    def coverage_plot(self):
        for each_chromosome in self.genome_analysis.keys():
            plot_cnv = CoveragePlot()
            plot_cnv.file_prefix = f'plot1-coverage-{self.sample_id}'
            plot_cnv.output_directory = self.output_directory
            chr_cov_data = self.genome_analysis[each_chromosome]["cov_data"]
            plot_cnv.plot_coverage(each_chromosome, {"pos": chr_cov_data.positions, "cov": chr_cov_data.raw_coverage})

    def cnv_plot(self, methode=""):
        for each_chromosome in self.genome_analysis.keys():
            chr_cov_data = self.genome_analysis[each_chromosome]["cov_data"]
            cov_stats = self.genome_analysis[each_chromosome]["norm_statistics"]
            lower_bound = self.lower_2n_threshold
            upper_bound = self.upper_2n_threshold
            cov_stats.median = cov_stats.median * self.ploidy  # for diploid
            new_plot_device = CNVPlot()
            new_plot_device.output_directory = self.output_directory
            new_plot_device.file_prefix = methode + self.sample_id
            new_plot_device.plot_coverage_cnv(each_chromosome, cov_stats,
                                              {"pos": chr_cov_data.positions,
                                               "cov": chr_cov_data.normalized_cov_ploidy},
                                              self.cnv_calls_list[each_chromosome], [lower_bound, upper_bound])

    def convert_vcf_to_tabular(self, snv_af_bed_output):
        self.logger.info("Parsing VCF to BED file")
        vcf_parser = vcf.VCFSNVParser(self.min_chr_length, self.as_dev)
        self.logger.debug("Parsing: VCF -> DataFrame")
        vcf_df = vcf_parser.vcf_to_dataframe(self.snv_file)
        self.logger.debug("Parsing: DataFrame -> BED")
        self.snv_af_df = vcf_parser.dataframe_to_tabular_file(vcf_df, self.coverage_file, snv_af_bed_output)
        self.snv_af_bed = snv_af_bed_output

    def call_cnv_af_region(self):
        # self.cnv_calls_list[each_chromosome]
        cnv_by_af = CNVAnalysisSNP(genome_info=self.genome_info,
                                   output_directory=self.output_directory,
                                   sample_id=self.sample_id,
                                   as_dev=self.as_dev)
        self.logger.info("Calculating CNV events based on SNV data")
        self.cnv_calls_list = cnv_by_af.af_cnv_call_region(self.cnv_calls_list, self.snv_af_bed)
        cnv_by_af.call_cnv_af(self.snv_af_df, write_csv=True)

    def get_cnv_metrics(self, refined_cnvs: bool = False):
        """
        Calculates for every CNV call a metrics like the p value with statistical tests
        :return:
        """
        self.cnv_metrics.evaluate_cnvs(self.cnv_calls_list, refined_cnvs)

    def write_intermediate_candidates(self, candidate_type: str = "") -> str:
        """
        Writing intermediate object files out
        :param candidate_type:
        :return:
        """

        intermediate_output_writer = util.outputWriter.IntermediateFile(self.output_directory)
        #genome_info = intermediate_output_writer.convert_candidates_to_dictionary(self.genome_info)
        cnv_calls_list_dict = intermediate_output_writer.convert_candidates_to_dictionary(self.cnv_calls_list)
        raw_cnv_calls_list_dict = intermediate_output_writer.convert_candidates_to_dictionary(self.raw_cnv_calls_list)

        analysis_dict = {
                "metadata": {"source":"spectre","spectre_version": "0.1"},
                "genome_info": self.genome_info,
                "raw_cnvs": raw_cnv_calls_list_dict,
                "refined_cnvs": cnv_calls_list_dict,
                "analysis_metrics": {
                    "min_chr_length": self.min_chr_length,
                    "max_std_outlier_rm": self.max_std_outlier_rm,
                    "mosdepth_cov_genome_chr_diff": self.mosdepth_cov_genome_chr_diff,
                    "lower_2n_threshold": self.lower_2n_threshold,
                    "upper_2n_threshold": self.upper_2n_threshold,
                    "cov_diff_threshold": self.cov_diff_threshold,
                    "dist_proportion": self.dist_proportion,
                    "dist_min_overwrite": self.dist_min_overwrite,
                    "candidate_final_threshold": self.candidate_final_threshold,
                    "genome_mean": np.mean(self.cnv_metrics.df_coverage_candidate_no_excl_zone_random_samples['coverage']),
                    "genome_sd": np.std(self.cnv_metrics.df_coverage_candidate_no_excl_zone_random_samples['coverage']),
                    "genome_var": np.var(self.cnv_metrics.df_coverage_candidate_no_excl_zone_random_samples['coverage'])
                }
            }
        output_path = intermediate_output_writer.write_intermediate_file(analysis_dict, f"{self.sample_id}")
        self.intermediate_candidates_file_location = output_path
        return output_path

    # ############################################
    # dev
    def dev_write_csv(self, each_chromosome):
        csv_results = pd.DataFrame(
            data={"position": self.genome_analysis[each_chromosome]["cov_data"].positions,
                  "mosdepth_cov": self.genome_analysis[each_chromosome]["cov_data"].coverage_raw,
                  "norm_cov": self.genome_analysis[each_chromosome]["cov_data"].normalized_cov,
                  "ploidy_cov": self.genome_analysis[each_chromosome]["cov_data"].normalized_cov_ploidy
                  }
        )
        csv_results["chr"] = each_chromosome
        output_file = f"{self.debug_dir}/cnv_{self.sample_id}_byCoverage_chr{each_chromosome}.csv"
        self.logger.debug(f"Writing coverage to {output_file}")
        csv_results.to_csv(output_file, index=False)
