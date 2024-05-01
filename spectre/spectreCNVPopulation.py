import gzip
import os
import traceback
import json

import pandas as pd

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
        self.cnv_merged_call_list = {}  # Holds merged CNV calls
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
        cnv_reciprocal_overlap_1 = overlap / cnv1.size
        cnv_reciprocal_overlap_2 = overlap / cnv2.size
        return cnv_reciprocal_overlap_1 >= reciprocal_overlap and cnv_reciprocal_overlap_2 >= reciprocal_overlap

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
        :param differ_ratio: percentage of similarity
        :return: True if the size of the CNVs is similar.
        """
        return abs(cnv1.size - cnv2.size) / max(cnv1.size, cnv2.size) <= differ_ratio

    def search_for_supporting_candidates(self, candidates: dict) -> dict:
        """
        Search for supporting CNVs in all samples. Qualifying CNVs must meet the requirements of the reciprocal overlap,
         same CN status and same CN type.
        :param candidates: dictionary of candidates
        :return: dictionary of CNVs which are supported by other CNVs
        """
        # 1) get names of all  samples from the candidates and generate a sorted dictionary
        sample_dict = dict(sorted(candidates.items()))
        total_sample_len = len(sample_dict.keys())
        # 1.2) create a cnv_call_list and initialize it with all the chromosomes as keys and an empty set as value
        cnv_call_list = dict([(chrom, set()) for chrom in self.genome_info["chromosomes"]])
        # 2) walk through all samples
        for chrom in self.genome_info["chromosomes"]:
            self.logger.info(f"Searching for supporting CNVs in chromosome: {chrom}")
            # 3) walk through all chromosomes
            sample1_cnt = 0
            for sample_origin_name, sample_candidates in sample_dict.items():
                # 4) walk through all candidates
                if chrom not in sample_candidates.keys():
                    continue
                for candidate1 in sample_candidates[chrom]:
                    # 5) walk through all samples
                    sample2_cnt = 0
                    for sample_origin2_name, sample_candidates2 in sample_dict.items():
                        # self.logger.info(
                        #     f"Checking {sample_origin_name}:{sample1_cnt + 1}/{total_sample_len} vs. {sample_origin2_name}:{sample2_cnt + 1}/{total_sample_len}")
                        # 6) walk through all candidates
                        for candidate2 in sample_candidates2[chrom]:
                            # 6.1) check if the same sample id
                            if candidate1.id == candidate2.id:
                                continue
                            # 6.2) check if the same CN status
                            if candidate1.cn_status != candidate2.cn_status:
                                continue
                            # 6.3) check if the same CN type
                            if candidate1.type != candidate2.type:
                                continue
                            # 6.4) check if the candidates are overlapping
                            if self.candidate_overlapping_reciprocal(candidate1, candidate2, self.reciprocal_overlap):
                                # 7) Add to candidate1 at sample_origin2_name the candidate2
                                candidate1.support_cnv_calls[sample_origin2_name].add(candidate2)
                                # 8) Add candidate1 to the cnv_call_list
                                cnv_call_list[chrom].add(candidate1)
                        sample2_cnt += 1
                sample1_cnt += 1
        pass
        return cnv_call_list

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

    def convert_candidates_dict_to_dataframe(self, candidate_dict: dict):
        df = pd.DataFrame(
            columns=["sample_origin", "id", "chromosome", "start", "end", "size", "cn_status", "type", "chr_idx"])
        for sample_origin, sample_candidates in candidate_dict.items():
            df_sample_origin_tmp = pd.DataFrame(columns=df.columns)
            for candidates in sample_candidates.values():
                tmp_candidate_add_list = []
                for candidate_idx, candidate in enumerate(candidates):
                    tmp_candidate_add_list.append([candidate.sample_origin, candidate.id, candidate.chromosome,
                                                   candidate.start, candidate.end, candidate.size, candidate.cn_status,
                                                   candidate.type, candidate_idx])
                # create a dataframe and concat it to df_sample_origin_tmp
                df_sample_origin_tmp = pd.concat(
                    [df_sample_origin_tmp, pd.DataFrame(tmp_candidate_add_list, columns=df.columns)])
            # concat df_sample_origin_tmp to df
            df = pd.concat([df, df_sample_origin_tmp])
        df.reset_index(drop=True, inplace=True)
        return df

    def apply_similarity_merge(self):
        n_merges = 2
        n_merge_rounds = 1
        while n_merges > 1:
            self.logger.info(f"Starting merge round: {n_merge_rounds}")
            n_merges = self.merge_similar_candidates()
            n_merge_rounds += 1

    def merge_similar_candidates(self):
        """
        Merging similar CNVs together. Similar CNVs are CNVs which have the same size, start and end position and the
        same type. The CNVs are added to the base (first) CNV population_merge_candidate_list and removed from the final
        CNV candidates list. Thus, the merged CNVs will no longer be used in the population mode.
        :return:
        """
        # define wobble of start and end position which should be considered as the same starting point
        wobble_position = 20 * self.bin_size
        wobble_size = 4 * wobble_position
        merge_idx = {}  # dictionary with the base index as key and the merge candidates as values
        new_merge_idx = {}  # Key superseding index, value list of indices which are superseded,
        # added to the index and removed from the final_candidates
        # lookup dataframe with sample_origin, sample_id, chromosome, start, end, size, cn_status, type
        # convert dictionary to dataframe
        df_final_candidates = self.convert_candidates_dict_to_dataframe(self.final_candidates)
        # 1) for each variant type -> avoid interacting of variant types: [DEL, DUP] vs LOH
        for cn_type in df_final_candidates.type.unique():
            df = df_final_candidates[df_final_candidates.type == cn_type]
            # 2) for each chromosome
            for chrom in df.chromosome.unique():
                df_chrom = df[df.chromosome == chrom]
                # 3) for each CN
                for cn_status in df_chrom.cn_status.unique():
                    df_chrom_cn = df_chrom[df_chrom.cn_status == cn_status]
                    if len(df_chrom_cn) > 1:

                        df_chrom_cn = df_chrom_cn.sort_values(by=["start"])
                        df_chrom_cn.reset_index(drop=False, inplace=True)

                        df_chrom_cn["start_base_df_idx"] = df_chrom_cn.apply(lambda x: set([x["index"]]), axis=1)
                        df_chrom_cn["end_base_df_idx"] = df_chrom_cn.apply(lambda x: set([x["index"]]), axis=1)

                        # 4) Check wobble states loop through all start_wobble, end_wobble and index values
                        for row in df_chrom_cn.itertuples():
                            # check in all df_chrom_cn rows if any of the start values is within the start value of the
                            # row +/- wobble_position.
                            start = row.start
                            end = row.end
                            base_df_idx = row.index
                            # get all rows where the start is within the wobble position and add the base_df_idx to the
                            # array stored in  start_base_df_idx.
                            mask = (df_chrom_cn["start"].between(start - wobble_position, start + wobble_position))
                            # append the base_df_idx to the array stored in  start_base_df_idx
                            df_chrom_cn.loc[mask, "start_base_df_idx"] = df_chrom_cn.loc[
                                mask, "start_base_df_idx"].apply(lambda x: x.union(set([base_df_idx])))
                            # get all rows where the end is within the wobble position and add the base_df_idx to the
                            # array stored in  end_base_df_idx
                            mask = (df_chrom_cn["end"].between(end - wobble_position, end + wobble_position))
                            # append the base_df_idx to the array stored in  end_base_df_idx
                            df_chrom_cn.loc[mask, "end_base_df_idx"] = df_chrom_cn.loc[mask, "end_base_df_idx"].apply(
                                lambda x: x.union(set([base_df_idx])))

                        # 5) Determine is the start and end position are the same by intersecting the
                        # start_base_df_idx and end_base_df_idx
                        df_chrom_cn["same_pos"] = False
                        df_idx_set = {}
                        for row in df_chrom_cn.itertuples():
                            start_idx_set = row.start_base_df_idx
                            end_idx_set = row.end_base_df_idx
                            overlapping_idx = tuple(start_idx_set.intersection(end_idx_set))
                            if len(overlapping_idx) > 1:
                                df_chrom_cn.at[row.Index, "same_pos"] = True
                                # add overlapping_idx to df_idx_set the smallest index is the base index
                                # get smallest index
                                base_idx = min(overlapping_idx)
                                # get all other indices
                                superseded_idx = list(overlapping_idx)
                                superseded_idx.remove(base_idx)
                                df_idx_set[base_idx] = superseded_idx
                        new_merge_idx |= df_idx_set  # merging the dictionaries

        # 6) merge the CNVs
        remove_candidate_list = set()
        for base_idx, superseded_idx_list in new_merge_idx.items():
            df_base_candidate = df_final_candidates.loc[base_idx]
            base_candidate = self.final_candidates[df_base_candidate["sample_origin"]][df_base_candidate["chromosome"]][
                df_base_candidate["chr_idx"]]
            self.logger.debug(f"Collapsing CNVs to: {base_candidate.id}")
            quality_list = [base_candidate.quality]
            for idx in superseded_idx_list:
                # merge the CNV
                df_merge_candidate = df_final_candidates.loc[idx]
                merge_candidate = \
                    self.final_candidates[df_merge_candidate["sample_origin"]][df_merge_candidate["chromosome"]][
                        df_merge_candidate["chr_idx"]]
                merge_candidate.merged_sample = True  # Mark the variant as merged
                self.logger.debug(f"Removing CNV: {merge_candidate.id}")
                # add to the support list
                # base_candidate.population_merge_candidate_list.append(merge_candidate)
                base_candidate.support_cnv_calls[merge_candidate.sample_origin].add(merge_candidate)
                # add candidates which should be removed
                quality_list.append(merge_candidate.quality)
                remove_candidate_list.add(merge_candidate)
            # 6.1) calculate new quality by choosing the best quality
            base_candidate.quality = max(quality_list)
            # 6.2) add the base_candidate to the self.cnv_merged_call_list
            if base_candidate.chromosome not in self.cnv_merged_call_list.keys():
                self.cnv_merged_call_list[base_candidate.chromosome] = []
            self.cnv_merged_call_list[base_candidate.chromosome].append(base_candidate)
            pass

        # 7) Remove all candidates which are in the remove_candidate_list
        for candidate in remove_candidate_list:
            self.final_candidates[candidate.sample_origin][candidate.chromosome].remove(candidate)
        # self.logger.info(f"Removed {len(merge_idx)} similar CNVs")
        self.logger.info(f"Merged {len(remove_candidate_list)} similar CNVs")
        return len(new_merge_idx.keys())

    def initialize_support_cnv_calls(self, candidates_dict: dict) -> dict:
        """
        Initialize the support CNV calls in the candidates dictionary.
        :param candidates_dict: dictionary of candidates
        :return: dictionary of candidates with initialized support CNV calls
        """
        sample_dict = dict(sorted(candidates_dict.items()))
        for sample_origin, sample_candidates in candidates_dict.items():
            for chrom, candidates in sample_candidates.items():
                for candidate in candidates:
                    candidate.support_cnv_calls = dict([(key, set()) for key in sample_dict.keys()])
        return candidates_dict

    def add_merged_candidates_to_support_list(self) -> None:
        """
        Add merged candidates to the support list of the candidate in the self.cnv_call_list. A merged candidate is only
        added, if the candidate is not already in the list. If the candidate is already in the list, the supporting CNVs
        from the merged candidate are added to the support list of the candidate in the self.cnv_call_list.
        :return:
        """
        for chrom, candidates in self.cnv_merged_call_list.items():
            for candidate in candidates:
                # if the chromosome is in the cnv_call_list
                if candidate.chromosome not in self.cnv_call_list.keys():
                    self.cnv_call_list[candidate.chromosome] = []
                # Candidate was not found in the list and will be added
                if candidate not in self.cnv_call_list[candidate.chromosome]:
                    self.cnv_call_list[candidate.chromosome].add(candidate)
                    continue
                # Found the same candidate in the list
                # add the supporting CNVs to the candidate
                for supporting_cnv in self.cnv_call_list[candidate.chromosome]:
                    if supporting_cnv.id == candidate.id:
                        for sample_origin, supporting_candidates in candidate.support_cnv_calls.items():
                            supporting_cnv.support_cnv_calls[sample_origin] |= supporting_candidates
                            pass



    def call_cnv(self) -> None:
        """
        Starts CNV calling with CNVs from multiple samples.
        :return: None
        """
        # Initialize support CNV calls
        self.final_candidates = self.initialize_support_cnv_calls(self.final_candidates)
        self.logger.info("Merge similar variants together.")
        self.apply_similarity_merge()
        self.logger.info(f"Starting population mode with samples: {', '.join(list(self.final_candidates.keys()))}")
        # generating union table of final overlaps
        # self.call_cnv_final_candidates()
        self.cnv_call_list = self.search_for_supporting_candidates(self.final_candidates)
        # Add merged candidates to the support list
        self.add_merged_candidates_to_support_list()

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
