import gzip
import os
import traceback
import json

from spectre.util import logger
from spectre.classes.loh_candidate import LoHCandidate
from spectre.util.outputWriter import VCFOutput
from spectre.classes.cnv_candidate import CNVCandidate


class SpectrePopulation(object):
    def __init__(self, sample_id, output_dir, genome_info, reciprocal_overlap=0.8, discard_quality_control=False,
                 as_dev=False):

        self.logger = logger.setup_log(__name__, as_dev)
        self.sample_id = sample_id
        self.output_dir = output_dir
        self.final_candidates = {}  # Holds all CNVs calls that have been analyzed and filtered by Spectre
        self.raw_candidates = {}  # Holds all raw CNVs calls before any filtering was applied by Spectre
        self.cnv_call_list = {}  # Holds resulting cNV calls.
        self.genome_info = genome_info
        self.bin_size = 0  # binsize of all samples
        self.reciprocal_overlap = reciprocal_overlap  # the minimum position overlap to be considered a support CNV call
        self.discard_quality_control = discard_quality_control  # if True, all CNVs that do not pass the quality
        # control will be considered for the lookup
        self.as_dev = as_dev  # TODO get outside propagation

    def merge_genome_info(self, new_genome_info: dict) -> None:
        """
        Merging genome information incase different genomic information is stored in the .spc files.
         This ensures the correct metadata information in the resulting population output VCF.
        :param new_genome_info: dict of genome information
        :return: None
        """
        if not self.genome_info:
            self.genome_info = new_genome_info.copy()
        else:
            if self.genome_info.keys() != new_genome_info.keys():
                logger.warning("Genomics information of provided samples are not matching!")
            for key1, values1 in self.genome_info.items():
                for key2, values2 in new_genome_info.items():
                    if key1 == key2:
                        for value2 in values2:
                            if value2 not in values1:
                                values1.append(value2)

    def load_files(self, files: list) -> None:
        """
        Determination of loading process for the candidate generation.
        :param files: list of paths to the .spc
        :return: None
        """
        file_counter = 1
        for file in files:
            try:
                if os.path.exists(file):
                    if str(file).__contains__(".spc") or str(file).__contains__(".spc.gz"):
                        self.load_candidates_from_intermediate_file_spc(file, file_counter)
                else:
                    raise
            except:
                self.logger.error(f"File does not exist! Provided : {file}")
            file_counter += 1
        pass

    def load_candidates_from_intermediate_file_spc(self, file, file_counter) -> None:
        """
        Load candidates from a json format file and converting the content to candidates.
        :param file: path to pickled file
        :return: None
        """
        spc_dict = dict()
        try:
            if '.gz' in file:
                with gzip.open(file, "rt", encoding="UTF-8") as input_file:
                    self.logger.info(f"Loading file: {file}")
                    spc_dict = json.load(input_file)
            else:
                with open(file, "rt") as input_file:
                    spc_dict = json.load(input_file)
            # Convert dictionary to candidates
            filename = os.path.basename(file).split(".spc")[0]
            sample_name = spc_dict["metadata"]["sample_id"]
            # initialize bin size
            if self.bin_size < 1:
                self.bin_size = int(spc_dict["metadata"]["bin_size"])
            # check if sample bin size is matching the other samples
            if self.bin_size != int(spc_dict["metadata"]["bin_size"]):
                self.logger.warning(
                    f"Bin size of {filename} ({sample_name}) does not match the bin size of the other samples! ")
                self.logger.warning(f" Skipping: {filename} ({sample_name})")
            else:
                self.convert_spc_to_candidate_list(spc_dict, file_counter)
        except:
            self.logger.error(traceback.print_exc())
            self.logger.error(f"Check if file meets the JSON standard. Error in file {file}")
        pass

    def convert_dict_to_candidate_list(self, filename, candidates_dict, candidate_class: str = "CNVCandidate") -> dict:
        result = dict([(chrom, []) for chrom in candidates_dict.keys()])
        for chrom, candidates in candidates_dict.items():
            for candidate in candidates:
                if candidate_class == "LoHCandidate":
                    new_candidate = LoHCandidate()
                else:
                    new_candidate = CNVCandidate()

                for candidate_key, candidate_value in candidate.items():
                    if candidate_key in new_candidate.__dict__.keys():
                        new_candidate.__setattr__(candidate_key, candidate_value)
                        new_candidate.__setattr__('sample_origin', filename)  # update source to filename
                        pass
                result[chrom].append(new_candidate)
        return result

    def convert_spc_to_candidate_list(self, candidate_dict: dict, file_counter: int):
        filename = candidate_dict["metadata"]["sample_id"]
        if filename in self.final_candidates.keys():
            filename = f'{filename}_{file_counter}'
        if "spectre" not in candidate_dict["metadata"]["source"]:
            self.logger.warning("Provided .spc file does not originate from Spectre.")
            self.logger.warning("Trying to convert the provided file")

        # Get genome information
        if "genome_info" in candidate_dict.keys():
            self.merge_genome_info(candidate_dict['genome_info'])
        # Refined CNVs
        # Get final/refined cnvs
        if 'refined_cnvs' in candidate_dict.keys():
            # Init candidate dictionary
            if filename not in self.final_candidates.keys():
                self.final_candidates[filename] = dict()
            self.final_candidates[filename] = self.convert_dict_to_candidate_list(filename,
                                                                                  candidate_dict['refined_cnvs'])
        if 'refined_loh_cnvs' in candidate_dict.keys():
            # Init candidate dictionary
            if filename not in self.final_candidates.keys():
                self.final_candidates[filename] = dict()
            refined_loh_candidates = self.convert_dict_to_candidate_list(filename, candidate_dict['refined_loh_cnvs'],
                                                                         candidate_class="LoHCandidate")
            self.final_candidates[filename] = self.merge_candidate_dicts(self.final_candidates[filename],
                                                                         refined_loh_candidates)

        # Raw CNVs
        if self.discard_quality_control:
            # Get raw cnvs
            if 'raw_cnvs' in candidate_dict.keys():
                if filename not in self.raw_candidates.keys():
                    self.raw_candidates[filename] = dict()
                self.raw_candidates[filename] = self.convert_dict_to_candidate_list(filename,
                                                                                    candidate_dict['raw_cnvs'])

            else:
                self.logger.warning("No raw CNVs found in .spc file.")
            if 'raw_loh_cnvs' in candidate_dict.keys():
                if filename not in self.raw_candidates.keys():
                    self.raw_candidates[filename] = dict()
                raw_loh_cnv_candidates = self.convert_dict_to_candidate_list(filename, candidate_dict['raw_loh_cnvs'],
                                                                             candidate_class="LoHCandidate")

                self.raw_candidates[filename] = self.merge_candidate_dicts(self.raw_candidates[filename],
                                                                           raw_loh_cnv_candidates)
        pass

    def merge_candidate_dicts(self, dict1, dict2):
        """
        Merging two dictionaries
        :param dict1: dictionary 1
        :param dict2: dictionary 2
        :return: merged dictionary
        """

        result = dict1.copy()
        for key2, values2 in dict2.items():
            if key2 not in dict1.keys():
                result[key2] = values2
            else:
                result[key2] += values2
        return result



    def cnv_call_population(self) -> None:
        """
        Starts CNV population calling
        :return: None
        """
        self.logger.info("Starting population CNV calls")
        self.call_cnv()

    @staticmethod
    def candidate_overlapping(cnv1: CNVCandidate, cnv2: CNVCandidate) -> bool:
        """
        Checking if two CNV candidates are overlapping
        :param cnv1: CNV candidate 1
        :param cnv2: CNV candidate 2
        :return: True if candidates are overlapping.
        """
        return cnv1.start <= cnv2.start <= cnv1.end or cnv1.start <= cnv2.end <= cnv1.end

    @staticmethod
    def candidate_overlapping_reciprocal(cnv1: CNVCandidate, cnv2: CNVCandidate, reciprocal_overlap: float) -> bool:
        """
        Checking if two CNV candidates are overlapping reciprocal. Finding the maximum starting position, minimum end
        position and calculating the overlap. The overlap is divided by the minimum size of the two CNVs.
        :param cnv1: CNV candidate 1
        :param cnv2: CNV candidate 2
        :param reciprocal_overlap: percentage of reciprocal overlap
        :return: True if candidates are reciprocal overlapping with a certain overlap.
        """
        overlap = min(cnv1.end, cnv2.end) - max(cnv1.start, cnv2.start)
        # if the overlap is negative the CNVs do not overlap
        if overlap <= 0:
            return False
        cnv_reciprocal_overlap = overlap / min(cnv1.size, cnv2.size)
        return cnv_reciprocal_overlap >= reciprocal_overlap

    @staticmethod
    def candidate_same_cn(cnv1: CNVCandidate, cnv2: CNVCandidate) -> bool:
        """
        Checking if two CNV candidates have the same copy number status
        :param cnv1: CNV candidate 1
        :param cnv2: CNV candidate 2
        :return: True if the copy number status matches.
        """
        return cnv1.cn_status == cnv2.cn_status

    @staticmethod
    def candidate_same_cn_type(cnv1: CNVCandidate, cnv2: CNVCandidate) -> bool:
        """
        Checking if two CNV candidates have the same copy number type
        :param cnv1: CNV candidate 1
        :param cnv2: CNV candidate 2
        :return: True if the copy number type matches.
        """
        return cnv1.type == cnv2.type

    @staticmethod
    def candidate_similar_size(cnv1: CNVCandidate, cnv2: CNVCandidate, differ_ratio: float) -> bool:
        """
        Checking if two CNV candidates have a similar size. The size of the CNVs is divided by the maximum size of the
        two CNVs.
        :param cnv1: CNV candidate 1
        :param cnv2: CNV candidate 2
        :param similarity_ratio: percentage of similarity
        :return: True if the size of the CNVs is similar.
        """
        return abs(cnv1.size - cnv2.size) / max(cnv1.size, cnv2.size) <= differ_ratio

    def call_cnv_final_candidates(self) -> None:
        """
        Creating a structure which holds all overlapping CNVs of the final CNV candidates.
        :return: None
        """
        self.logger.info("Searching for overlapping CNVs in final candidates")
        sample1_cnt = 0
        sample2_cnt = 0
        final_chr_pos_cn = set()
        # everything against everything
        for samples_key1, samples_values1 in self.final_candidates.items():
            self.logger.info(
                f"Checking final CNVs {samples_key1}: {sample1_cnt + 1}/{len(self.final_candidates.keys())}")
            sample2_cnt = 0
            for samples_key2, samples_values2 in self.final_candidates.items():
                for sample_chrom1, sample_values1 in samples_values1.items():
                    for sample_chrom2, sample_values2 in samples_values2.items():
                        # self.logger.info(f"Checking {sample1_cnt}:{sample_chrom1}:{samples_key1} against {sample2_cnt}:{sample_chrom2}:{samples_key2}")
                        # check for same chromosome
                        if sample_chrom1 != sample_chrom2:
                            continue
                        # all individual samples against each other
                        for sample1 in sample_values1:
                            for sample2 in sample_values2:
                                # check if chr, pos, cn are already in the list
                                if f"{sample1.chromosome}{sample1.start}{sample1.cn_status}" in final_chr_pos_cn:
                                    continue
                                else:
                                    final_chr_pos_cn.add(f"{sample1.chromosome}{sample1.start}{sample1.cn_status}")
                                differ_ratio = 1 - self.reciprocal_overlap
                                is_similar_size = self.candidate_similar_size(cnv1=sample1, cnv2=sample2,
                                                                              differ_ratio=differ_ratio)
                                is_overlapping = self.candidate_overlapping_reciprocal(sample1, sample2,
                                                                                       self.reciprocal_overlap)
                                is_same_cn = self.candidate_same_cn(sample1, sample2)
                                is_same_cn_type = self.candidate_same_cn_type(sample1, sample2)
                                if is_overlapping and is_same_cn and is_same_cn_type and is_similar_size:

                                    # check if all keys are available in sample1
                                    if sample1.support_cnv_calls.keys() != self.final_candidates.keys():
                                        for key in self.final_candidates.keys():
                                            if key not in sample1.support_cnv_calls.keys():
                                                sample1.support_cnv_calls[key] = set()
                                        # sort dictionary by key
                                        sample1.support_cnv_calls = dict(sorted(sample1.support_cnv_calls.items()))

                                    # add sample2 to the support cnvs in sample1
                                    sample1.support_cnv_calls[samples_key2].add(sample2)

                                    # check if chromosome exists in cnv_call_list
                                    if sample_chrom1 not in self.cnv_call_list.keys():
                                        self.cnv_call_list[sample_chrom1] = []  # create list for

                                    # check if sample1 is in the list
                                    if sample1 not in self.cnv_call_list[sample_chrom1]:
                                        self.cnv_call_list[sample_chrom1].append(sample1)

                    sample2_cnt += 1
            sample1_cnt += 1

    def cnv_lookup_in_raw_candidates(self) -> None:
        """
        Check if any of the final CNVs are supported by any non-refined CNVs. If a CNV is found it will be
        added to the candidate.
        :return: None
        """
        self.logger.info("Searching for overlapping CNVs in raw CNVs")
        variant_cnt = 0
        variant_len = len(self.cnv_call_list.keys())
        raw_candidate_len = len(self.raw_candidates.keys())
        candidate_cnt = 0
        missed_cnt = 0
        for variant_key, variants in self.cnv_call_list.items():
            # self.logger.info(f"Checking {variant_key}:{variant_cnt + 1}/{variant_len}")
            for variant in variants:
                candidate_cnt = 0
                for raw_sample_key, value in self.raw_candidates.items():
                    for chrom_key, candidates in value.items():
                        if chrom_key != variant.chromosome:
                            continue
                        self.logger.info(
                            f"Final variant:{variant_key}:{variant_cnt + 1}/{variant_len} vs. raw:{candidate_cnt + 1}/{raw_candidate_len}")
                        for candidate in candidates:
                            if variant.sample_origin != candidate.sample_origin:
                                # qualification checks
                                is_overlapping = self.candidate_overlapping_reciprocal(variant, candidate,
                                                                                       self.reciprocal_overlap)
                                is_same_cn = self.candidate_same_cn(variant, candidate)
                                is_same_cn_type = self.candidate_same_cn_type(variant, candidate)
                                if is_overlapping and is_same_cn and is_same_cn_type:
                                    variant.support_cnv_calls[candidate.sample_origin].add(candidate)
                                    missed_cnt += 1
                    candidate_cnt += 1
            variant_cnt += 1

    def call_cnv(self) -> None:
        """
        Starts CNV calling with CNVs from multiple samples.
        :return: None
        """
        self.logger.info(f"Starting population mode with samples: {', '.join(list(self.final_candidates.keys()))}")
        # generating union table of final overlaps
        self.call_cnv_final_candidates()
        if self.discard_quality_control:
            self.logger.info("Quality control is disabled. Searching also in raw CNVs for supporting CNVs.")
            # look up if all missing fields are covered by any raw cnv call
            self.cnv_lookup_in_raw_candidates()
        # Writing results to disk
        output_file = f"{self.output_dir}/population_mode_{self.sample_id}.vcf.gz"
        self.logger.info(f"Writing population VCF @: {output_file}")
        vcf_output = VCFOutput(output_file=output_file, genome_info=self.genome_info)
        vcf_output.make_vcf(
            chromosome_list=self.cnv_call_list.keys(), cnv_candidate_list=self.cnv_call_list,
            sample_id=self.final_candidates.keys(), population_sample_ids=list(self.final_candidates.keys()))
