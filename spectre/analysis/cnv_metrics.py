import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import hashlib
import random
from scipy.stats import ks_2samp, norm
from spectre.util import logger


class CNVMetrics(object):
    def __init__(self, coverage_analysis, exclusion_zones, cnv_calls=None,
                 hashname: str = "mock_dataname_of_coverage.csv",ploidy:float=2.0,
                 output_dir: str = "", is_cancer:bool =False, as_dev: bool = False, debug_dir=""):
        self.as_dev = as_dev
        self.logger = logger.setup_log(__name__, self.as_dev)
        self.coverage_analysis = coverage_analysis
        self.cnv_calls = cnv_calls
        self.hashname = hashname
        self.ploidy = ploidy
        self.is_cancer = is_cancer
        self.output_dir = output_dir
        # conversion
        self.df_mosdepth_coverage = None
        self.__convert_coverage_analysis_to_coverage_dataframe()
        self.blacklist_exclusions = None
        self.__convert_blacklist_exclusion_zone(exclusion_zones)
        self.__prepare_cnv_evaluation()
        self.cnv_exclusion = None
        # calculate threholds for deletion and duplication
        self.del_threshold, self.dup_threshold = None, None
        self.strict_del_threshold, self.strict_dup_threshold = None, None
        self.debug_dir = debug_dir

    def __convert_coverage_analysis_to_coverage_dataframe(self):
        """
        Convert the standard Spectre coverage_analysis into a pandas dataframe, using the minimal required fields
        to evaluate the CNV calls.
        :param coverage_analysis: standard Spectre coverage_analysis dictionary
        :return: pandas dataframe holding at least column values ["chr", "position", "coverage"]
        """
        self.logger.debug("Convert Genome analysis to coverage dataframe")
        list_modsepth_coverage = [[key, self.coverage_analysis[key]["cov_data"].positions,
                                   self.coverage_analysis[key]["cov_data"].normalized_cov_ploidy] for key in
                                   self.coverage_analysis.keys()]
        self.df_mosdepth_coverage = pd.DataFrame(list_modsepth_coverage)
        self.df_mosdepth_coverage.columns = ["chr", "position", "coverage"]
        self.df_mosdepth_coverage = self.df_mosdepth_coverage.explode(["position", "coverage"])
        self.df_mosdepth_coverage.reset_index(inplace=True)

    def __convert_cnv_calls_to_exclusions(self):
        """
        Converts the Spectre into a specific form of dictionary, which will be used to perform the statistical
        analysis. Furthermore, the cnv calls play a role in randomly sampling coverage values which will be used
        in the statistical analysis.
        """
        self.logger.debug("Converting cnv calls to exclusion zone")
        result = {}
        for key in self.cnv_calls.keys():
            if key not in result.keys():
                result[key] = []
            for candidate in self.cnv_calls[key]:
                result[key].append({"start": candidate.start, "end": candidate.end, "info": "cnv"})
        self.cnv_exclusion = result

    def __convert_blacklist_exclusion_zone(self, blacklist: dict):
        """
        Converts the given blacklist which contains a dict with lists of tupels into a dict with lists
        and items as dicts.
        :param blacklist: standard Spectre blacklist
        :return:
        """
        self.logger.debug("Converting blacklist to exclusion zone")
        result = {}
        for key in blacklist.keys():
            if key not in result.keys():
                result[key] = []
            # convert tupel into dictionary
            for start, end in blacklist[key]:
                result[key].append({"start": start, "end": end, "info": "blacklist"})
        self.blacklist_exclusions = result

    def __merge_dict_with_lists_as_value(self, dict_1: dict, dict_2: dict) -> dict:
        """
        Merges two dictionaries which hold a list as value and creating a completely new dictionary
        in which the lists are appended to the other one.
        :param dict_1: primary dictionary. (values will be first in the merged list)
        :param dict_2: secondary dictionary. (values will be appended to the dict_1 values)
        :return: merged and newly created dictionary
        """
        self.logger.debug("Merging exclusion zones")
        dict_result = dict_1.copy()
        for key in dict_2:
            if key not in dict_1.keys():
                dict_result[key] = dict_2[key].copy()
            else:
                dict_result[key] = [*dict_result[key], *dict_2[key].copy()]
        return dict_result

    def get_exclusion_zones_indices_in_coverage_data(self, df_mosdepth_coverage, dict_excl_zone):
        """
        Gets all mosdepth coverage dataframe indices, which are within the exclusion zone.
        :param df_mosdepth_coverage: dataframe with mosdepth coverage
        :param dict_excl_zone: list of indices
        :return:
        """
        self.logger.debug("Getting exclusion zone indices in mosdepth coverage data")
        excl_indices = []
        for excl_key in dict_excl_zone.keys():
            df_chr = df_mosdepth_coverage.loc[df_mosdepth_coverage.chr == excl_key]
            for excl_zone in dict_excl_zone[excl_key]:
                in_excl_zone = df_chr.loc[
                    df_chr.position.between((int(excl_zone["start"]) - 1), int(excl_zone["end"]))].index  # ["index"]

                if len(in_excl_zone) > 0:
                    excl_indices = excl_indices + list(in_excl_zone)
        return excl_indices

    def random_sample_selection(self, df_source, amount: float = 0.1):
        """
        Selecting random samples outside any exclusion zone (blacklist, cnv, ...)
        :param df_source:
        :param amount:
        :return:
        """
        self.logger.debug("Generating random coverage sample indices")
        result_hash_of_filename = hashlib.md5(self.hashname.encode()).hexdigest()  # used for the random samples
        random.seed(result_hash_of_filename)

        fraction = amount
        last_index = len(df_source.index) - 1
        cov_window_amount = int(last_index * fraction)

        samples_indices = random.sample(range(0, last_index), cov_window_amount)
        sorted_samples_indices = sorted(samples_indices)
        return sorted_samples_indices

    def get_ks_test(self, cnv_coverage, df_random_sample_coverage):
        """
        Calculating the ks-Score for given cnv coverage values and comparing the mean and sd to the distribution of
        randomly samples coverage datasets
        :param cnv_coverage: list of coverage values
        :param df_random_sample_coverage: dataframe holding randomly selected coverage data
        :return: dictionary with: cnv_call_mean, score (z-score), pvalue, ...
        """
        # self.logger.debug("Performing the KS-test")
        result = {}

        cnv_coverage = np.array(cnv_coverage)
        cnv_coverage = cnv_coverage[~np.isnan(cnv_coverage)]
        random_sample_coverage = df_random_sample_coverage["coverage"]
        random_sample_coverage = random_sample_coverage.dropna()
        statistic, pvalue = ks_2samp(cnv_coverage, random_sample_coverage)

        result["cnv_call_mean"] = round(np.nanmean(cnv_coverage), 3)
        result["pvalue"] = pvalue
        result["sample_score"] = abs(pvalue * 255)
        result["score"] = np.nan
        result["statistics"] = statistic
        result["test_type"] = "ks-score"
        return result

    def get_z_score(self, cnv_coverage, df_random_sample_coverage):
        """
        Calculating the Z-Score for given cnv coverage values and comparing the mean and sd to the distribution of
        randomly samples coverage datasets
        :param cnv_coverage: list of coverage values
        :param df_random_sample_coverage: dataframe holding randomly selected coverage data
        :return: dictionary with: cnv_call_mean, score (z-score), pvalue, ...
        """
        # self.logger.debug("Performing the Z-Test")
        result = {}
        cnv_coverage = np.array(cnv_coverage)
        cnv_coverage = cnv_coverage[~np.isnan(cnv_coverage)]  # exclude NaN values
        cnv_coverage = cnv_coverage[~np.isinf(cnv_coverage)]  # exclude -inf/inf values

        random_sample_coverage = df_random_sample_coverage["coverage"]
        random_sample_coverage.replace([np.inf, -np.inf], np.nan, inplace=True)  # replace -inf/inf values with nan
        random_sample_coverage = random_sample_coverage.dropna()  # exclude all nan values

        mean = np.mean(cnv_coverage)
        sd = np.std(cnv_coverage)

        # Z score
        z = (np.mean(random_sample_coverage) - mean) / (sd / np.sqrt(len(random_sample_coverage)))
        # pvalue
        z_capped = min(z, 10) if z > 1 else max(-10, z)
        # pvalue = 2 * (1 - norm.cdf(abs(z_capped)))  # Two-tailed test
        pvalue = norm.sf(abs(z_capped)) * 2
        # pvalue = 2 * (1 - norm.cdf(abs(z)))  # Two-tailed test
        gq = -np.log10(pvalue) * 10
        gq_capped = min(60, int(gq))

        result["cnv_call_mean"] = np.nanmean(cnv_coverage)
        result["score"] = z_capped
        result["pvalue"] = pvalue
        result["sample_score"] = gq_capped
        result["statistics"] = None
        result["test_type"] = "z-score"
        return result

    def __prepare_cnv_evaluation(self) -> None:
        """
        Preparing all the necessary variables, to performe the CNV evaluation on the cnv candidates
        :return: None
        """
        self.logger.debug("Preparing all parameters for the cnv evaluation")
        df_coverage_candidate = self.df_mosdepth_coverage.copy()
        df_coverage_candidate["blacklist"] = False  # holds only the values if row is in blacklist
        # df_coverage_candidate["cnv"] = False  # holds only values that are in cnv_events
        # df_coverage_candidate["excl_zone"] = False  # holds the combination of blacklist and cnv
        df_coverage_candidate.replace([np.inf, -np.inf], np.nan, inplace=True)
        df_coverage_candidate_no_Nans = df_coverage_candidate.dropna(subset=["coverage"])  # contains no nan values

        # get exclusion indices
        excl_zone_blacklist_indices = self.get_exclusion_zones_indices_in_coverage_data(df_coverage_candidate_no_Nans,
                                                                                        self.blacklist_exclusions)

        # Add exclusion indices
        df_coverage_candidate_no_Nans.iloc[
            excl_zone_blacklist_indices, df_coverage_candidate_no_Nans.columns.get_loc("blacklist")] = True

        # load curated data
        df_coverage_candidate_no_blacklist = df_coverage_candidate_no_Nans.loc[df_coverage_candidate.blacklist == False]

        # get random samples
        no_blacklist_indices = self.random_sample_selection(df_coverage_candidate_no_blacklist, 0.1)
        self.df_coverage_candidate_no_blacklist_random_samples = df_coverage_candidate_no_blacklist.iloc[
            no_blacklist_indices]

        # add for later use
        if self.cnv_calls:
            self.__recalculate_exlustion_zone()
        else:
            self.df_coverage_candidate_no_excl_zone_random_samples = self.df_coverage_candidate_no_blacklist_random_samples.copy()

        if self.as_dev:
            # write random samples to file
            self.logger.debug("CNV-Metrics: Writing random samples to file")
            self.df_coverage_candidate_no_blacklist_random_samples.to_csv(
                f"{self.output_dir}/debug/cnv_metrics_random_samples.csv.gz", sep="\t", index=False, compression="gzip")

    def __recalculate_exlustion_zone(self) -> None:
        """
        Recalculates the random sample coverage, based on new cnv_calls. This ensures that no coverage sample from
        the blacklist or an existing CNV call will be used to give the CNV calls a score.
        :param cnv_calls:
        :return: None
        """
        # setup
        # fill in self.cnv_exclusion
        self.__convert_cnv_calls_to_exclusions()
        self.blacklist_cnv_exclusion = self.__merge_dict_with_lists_as_value(self.blacklist_exclusions,
                                                                             self.cnv_exclusion)

        df_coverage_candidate = self.df_mosdepth_coverage.copy()
        # df_coverage_candidate["blacklist"] = False  # holds only the values if row is in blacklist
        # df_coverage_candidate["cnv"] = False  # holds only values that are in cnv_events
        df_coverage_candidate["excl_zone"] = False  # holds the combination of blacklist and cnv
        df_coverage_candidate.replace([np.inf, -np.inf], np.nan, inplace=True)
        df_coverage_candidate_no_Nans = df_coverage_candidate.dropna(subset=["coverage"])  # contains no nan values

        # get inclustio
        # excl_zone_cnv_indices = self.get_exclusion_zones_indices_in_coverage_data(
        #    df_coverage_candidate_no_Nans,self.cnv_exclusion)
        excl_zone_backlist_cnv_indices = self.get_exclusion_zones_indices_in_coverage_data(
            df_coverage_candidate_no_Nans, self.blacklist_cnv_exclusion)

        # Add exclusion indices
        # df_coverage_candidate_no_Nans.loc[excl_zone_cnv_indices, ["cnv"]] = True

        if len(excl_zone_backlist_cnv_indices) > 0:
            df_coverage_candidate_no_Nans.loc[excl_zone_backlist_cnv_indices, "excl_zone"] = True

        # df_coverage_candidate_no_cnv = df_coverage_candidate_no_Nans.loc[df_coverage_candidate.cnv == False]
        df_coverage_candidate_no_excl_zone = df_coverage_candidate_no_Nans.loc[df_coverage_candidate.excl_zone == False]

        # Get random samples for samples with the same random
        # no_cnv_indices = self.random_sample_selection(df_coverage_candidate_no_cnv, 0.1)
        # self.df_coverage_candidate_no_cnv_random_samples = df_coverage_candidate_no_cnv.iloc[no_cnv_indices]

        no_excl_zone_indices = self.random_sample_selection(df_coverage_candidate_no_excl_zone, 0.1)
        self.df_coverage_candidate_no_excl_zone_random_samples = df_coverage_candidate_no_excl_zone.iloc[
            no_excl_zone_indices]
        pass

    @staticmethod
    def get_lower_border(mean: float, sd: float, sd_multiplier: float, epsilon: float = 0.01,
                         epsilon_multiplier: float = 0):
        return mean - (sd * sd_multiplier + (epsilon * epsilon_multiplier))

    @staticmethod
    def get_upper_border(mean: float, sd: float, sd_multiplier: float, epsilon: float = 0.01,
                         epsilon_multiplier: float = 0):
        return mean + (sd * sd_multiplier + (epsilon * epsilon_multiplier))

    def calculate_del_dup_borders(self, del_border_multiplier: float = 1.0, dup_border_multiplier: float = 1.0,
                                  ploidy: float = 2.0) -> tuple:
        """
        Calculate dynamic borders for deletions (DEL) and duplications (DUP), based on the randomly selected coverage
        samples.
        :param del_border_multiplier: Multiplier for the tightness of the lower (deletion) border
        :param dup_border_multiplier: Multiplier for the tightness of the upper (duplication) border
        :return: Tupel with deletion border and duplication border (DEL,DUP)
        """
        # Find optimal border
        ploidy_offset_deletion_threshold = ploidy - 1
        ploidy_offset_duplication_threshold = ploidy + 1
        epsilon = 0.01
        step_cnt = 1
        max_steps = 1000
        found_new_deletion_border = False
        found_new_duplication_border = False
        ratio = 0.2  # default 0.4
        sd_clip_n = 3 # SD multiplier for clipping the random coverage samples

        mean = float(np.mean(self.df_coverage_candidate_no_excl_zone_random_samples['coverage']))
        # determine the SD of the random samples without considering the outliers
        sd_estimate = float(np.std(self.df_coverage_candidate_no_excl_zone_random_samples['coverage']))
        # select only coverage values which are less than the mean + (SD * sd_clip_n), to avoid extreme outliers
        df_random_clipped = self.df_coverage_candidate_no_excl_zone_random_samples['coverage'][
            self.df_coverage_candidate_no_excl_zone_random_samples['coverage'] < mean + (sd_clip_n * sd_estimate)]
        sd = float(np.std(df_random_clipped))
        var = np.var(self.df_coverage_candidate_no_excl_zone_random_samples['coverage'])
        self.logger.info(f"Random sample mean: {mean} \t SD: {sd} \t VAR: {var}")

        # Rough border estimation
        lower_threshold = mean - (sd * del_border_multiplier)
        upper_threshold = mean + (sd * dup_border_multiplier)
        pass

        # lower and upper threshold under or over ploidy -/+ 1
        if lower_threshold < ploidy_offset_deletion_threshold:
            self.logger.warning(f"Estimated DEL threshold {lower_threshold} is under the allowed minimum threshold of {ploidy_offset_deletion_threshold}. ")
            lower_threshold = (mean + ploidy_offset_deletion_threshold)/2
            self.logger.warning(f"Recalculating DEL threshold using fallback to threshold entrypoint set to {lower_threshold}")
        if upper_threshold > ploidy_offset_duplication_threshold:
            self.logger.warning(f"Estimated DUP threshold {upper_threshold} is over the allowed maximum threshold of {ploidy_offset_duplication_threshold}. ")
            upper_threshold = (mean + ploidy_offset_duplication_threshold)/2
            self.logger.warning(f"Recalculating DUP threshold using fallback to threshold entrypoint set to {upper_threshold}")

        while step_cnt < max_steps:
            # Calculate DEL border
            if not found_new_deletion_border:
                del_threshold_to_ploidy_offset = abs(lower_threshold - ploidy_offset_deletion_threshold) * (1 - ratio)
                del_threshold_to_ploidy = abs(lower_threshold - ploidy) * ratio
                self.logger.debug(f"DEL: {lower_threshold} \t DEL to ploidy offset: {del_threshold_to_ploidy_offset} \t DEL to ploidy: {del_threshold_to_ploidy}")

                if lower_threshold > ploidy_offset_deletion_threshold:
                    if del_threshold_to_ploidy_offset > del_threshold_to_ploidy:
                        lower_threshold = self.get_lower_border(mean, sd, del_border_multiplier, epsilon, step_cnt)
                    else:
                        self.logger.debug(f"Stopping with deletion threshold calculation. Found threshold: {lower_threshold}")
                        found_new_deletion_border = True
                else:
                    self.logger.debug(
                        f"Stopping with deletion threshold calculation. Under ploidy Found threshold: {lower_threshold}, at step: {step_cnt}")
                    lower_threshold = ploidy_offset_deletion_threshold
                    found_new_deletion_border = True

            # Calculate DUP border
            if not found_new_duplication_border:
                dup_thresholds_to_ploidy_offset = abs(upper_threshold - ploidy_offset_duplication_threshold) * (1 - ratio)
                dup_thresholds_to_ploidy = abs(upper_threshold - ploidy) * ratio
                self.logger.debug(f"DUP: {upper_threshold} \t DUP to ploidy offset: {dup_thresholds_to_ploidy_offset} \t DUP to ploidy: {dup_thresholds_to_ploidy}")

                if upper_threshold < ploidy_offset_duplication_threshold:
                    if dup_thresholds_to_ploidy_offset > dup_thresholds_to_ploidy:
                        upper_threshold = self.get_upper_border(mean, sd, dup_border_multiplier, epsilon, step_cnt)
                    else:
                        self.logger.debug(f"Stopping with duplication threshold calculation. Found threshold: {upper_threshold}")
                        found_new_duplication_border = True
                else:
                    self.logger.debug(
                        f"Stopping with duplication threshold calculation. Over ploidy Found threshold: {upper_threshold}, at step: {step_cnt}")
                    upper_threshold = ploidy_offset_duplication_threshold
                    found_new_duplication_border = True

            step_cnt += 1
        self.logger.debug(f"ploidy_offset_deletion_threshold: {ploidy_offset_deletion_threshold}")
        self.logger.debug(f"ploidy: {ploidy}")
        self.logger.debug(f"ploidy_offset_duplication_threshold: {ploidy_offset_duplication_threshold}")
        self.logger.debug(f"NEW Calculated thresholds for DEL: {lower_threshold} and DUP: {upper_threshold}")
        # save values
        self.del_threshold = lower_threshold
        self.dup_threshold = upper_threshold

    def calculate_del_dup_borders_sd(self, sd_clip_n: float = 2.0):
        """
        Determine the deletion and duplication thresholds based on the SD of the randomly selected coverage
        :param sd_clip_n: SD multiplier for clipping the random coverage samples
        :return:
        """
        random_sample_mean = float(np.mean(self.df_coverage_candidate_no_excl_zone_random_samples['coverage']))
        # determine the SD of the random samples without considering the outliers
        random_sample_sd = float(np.std(self.df_coverage_candidate_no_excl_zone_random_samples['coverage']))
        # select only coverage values which are less than the mean + (SD * sd_clip_n), to avoid extreme outliers

        del_threshold = random_sample_mean - (random_sample_sd * sd_clip_n)
        dup_threshold = random_sample_mean + (random_sample_sd * sd_clip_n)

        del_threshold, dup_threshold = self.check_del_dup_borders(del_threshold, dup_threshold)

        self.strict_del_threshold = del_threshold
        self.strict_dup_threshold = dup_threshold

        return del_threshold, dup_threshold
        # self.del_threshold = del_threshold
        # self.dup_threshold = dup_threshold


    def calculate_del_dup_borders_quantile(self, left_quantile: float = 0.1, right_quantile: float = 0.9):
        """
        Determine the deletion and duplication thresholds based on the quantiles of the randomly selected coverage
        :param left_quantile: DEL threshold quantile 10th quantile
        :param right_quantile: DUP threshold quantile 90th quantile
        :return:
        """
        self.logger.debug(f"Calculationg DEL/DUP thresholds based on quantiles: {left_quantile} and {right_quantile}")

        mean = float(np.mean(self.df_coverage_candidate_no_excl_zone_random_samples['coverage']))
        # determine the SD of the random samples without considering the outliers
        sd_estimate = float(np.std(self.df_coverage_candidate_no_excl_zone_random_samples['coverage']))
        # select only coverage values which are less than the mean + (SD * sd_clip_n), to avoid extreme outliers

        sd_clip_n = 2  # SD multiplier for clipping the random coverage samples
        df_random_clipped = self.df_coverage_candidate_no_excl_zone_random_samples['coverage'][
            self.df_coverage_candidate_no_excl_zone_random_samples['coverage'] < mean + (sd_clip_n * sd_estimate)]

        del_threshold = df_random_clipped.quantile(left_quantile)
        dup_threshold = df_random_clipped.quantile(right_quantile)

        del_threshold, dup_threshold = self.check_del_dup_borders(del_threshold, dup_threshold)
        self.logger.debug(f"Calculated thresholds for DEL: {del_threshold} and DUP: {dup_threshold}")
        self.del_threshold = del_threshold
        self.dup_threshold = dup_threshold

    def check_del_dup_borders(self, del_threshold: float, dup_threshold: float) -> tuple:
        self.logger.debug(f"CNV-Metrics: Calculated thresholds for DEL: {del_threshold} and DUP: {dup_threshold}")
        if not self.is_cancer:
            # If cancer mode is present no check of thresholds between ploidy -/+ 1 is performed
            if del_threshold < self.ploidy - 1:
                self.logger.warning(
                    f"Estimated DEL threshold {del_threshold} is under the allowed minimum threshold of {self.ploidy - 1}. ")
                del_threshold = self.ploidy - 1
                self.logger.warning(f"Using fallback to threshold of: {del_threshold}")
            if dup_threshold > self.ploidy + 1:
                self.logger.warning(
                    f"Estimated DUP threshold {dup_threshold} is over the allowed maximum threshold of {self.ploidy + 1}. ")
                dup_threshold = self.ploidy + 1
                self.logger.warning(f"Using fallback to threshold of: {dup_threshold}")
        return del_threshold, dup_threshold
    def evaluate_cnvs(self, cnv_calls=None, refined_cnvs: bool = False) -> dict:
        """
        Evaluating all submitted CNV calls, by applying statistical tests. Optionally, plotting the location of the CNVs
        on the global coverage samples per chr.
        :return: Modified CNVCandidates
        """
        self.logger.debug("Evaluating all CNVs")

        if self.cnv_calls is not None:
            self.logger.debug("Recalculating")
            self.__recalculate_exlustion_zone()

        self.logger.info(f"DEL border:{self.del_threshold} \t DUP border:{self.dup_threshold}")
        border_del = self.del_threshold
        border_dup = self.dup_threshold

        for cnv_chr in self.cnv_calls.keys():
            for cnv in self.cnv_calls[cnv_chr]:
                # Z-Test 2 - DIY
                if refined_cnvs:
                    cnv_metrics_Z = self.get_z_score(cnv.cov, self.df_coverage_candidate_no_excl_zone_random_samples)
                else:
                    cnv_metrics_Z = {}
                    cnv_metrics_Z["statistics"] = None
                    cnv_metrics_Z["score"] = 0
                    cnv_metrics_Z["pvalue"] = 0
                    cnv_metrics_Z["sample_score"] = 0
                    cnv_metrics_Z["cnv_call_mean"] = np.nanmean(cnv.cov)

                cnv.statistics["z-score"] = {"statistics": cnv_metrics_Z["statistics"],
                                             "score": cnv_metrics_Z["score"],
                                             "pvalue": cnv_metrics_Z["pvalue"],
                                             "sample_score": cnv_metrics_Z["sample_score"]}

                if self.as_dev and False:  # BUG line 385
                    if not np.isnan(cnv_metrics_Z["cnv_call_mean"]):

                        # pre check --> avoide div by 0
                        if cnv_metrics_Z["cnv_call_mean"] == 0:
                            x_index = 0
                        else:
                            if not np.isinf(cnv_metrics_Z["cnv_call_mean"]):
                                x_index = int(cnv_metrics_Z["cnv_call_mean"] / width)
                            else:
                                x_index = int(x[-1])
                        # check if out of bounds
                        if x_index >= len(x):
                            cnv_y = x[-1]
                        else:
                            cnv_y = int(x[x_index])
                        cnv_y = int(cnv_y)
                        annotation_endpoint_offset = 5 * (100 - x_index) * (
                                1 - (cnv_y / 150)) + (max(x) / 10)  # np.random.randint(100,300)

                        # Add
                        ax.annotate(f"{cnv.id}", xy=(cnv_metrics_Z["cnv_call_mean"], cnv_y),
                                    xytext=(cnv_metrics_Z["cnv_call_mean"], cnv_y + annotation_endpoint_offset),
                                    rotation=60,
                                    fontsize=12, arrowprops=dict(facecolor='green', shrink=0.05))
        if self.as_dev:
            self.logger.debug(f"Creating CNV distribution plot")
            fig, ax = plt.subplots()
            x, bins, p = ax.hist(self.df_coverage_candidate_no_excl_zone_random_samples['coverage'], bins=100,
                                 range=[0.0, 5.0])
            # FIXME: BUG -> AttributeError: 'silent_list' object has no attribute 'patches'
            # width = p.patches[0].get_width()

            ax.axvline(border_del, ymax=1, color="red", linestyle=":")
            ax.axvline(border_dup, ymax=1, color="green", linestyle=":")
            ax.axhline()
            ax.text(2.7 - 0.65, max(x), f"DUP: {round(self.dup_threshold, 4)}")
            ax.text(1.2 - 0.65, max(x), f"DEL: {round(self.del_threshold, 4)}")
            ax.set_title(f"Coverage distribution on random global bin samples")
            ax.set_xlabel('Coverage')
            ax.set_ylabel('Count')
            plot_path = f"{self.debug_dir}/genome_coverage_distribution-{self.hashname}.png"
            self.logger.debug(f'Plot path: {plot_path}')
            fig.savefig(plot_path)
