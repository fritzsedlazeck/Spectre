import pandas as pd
import numpy as np
import logging as logger


class CNV_Post_Processing(object):
    def __init__(self, ploidy_coverage, position, bin_size: int = 1000, as_dev: bool = False):
        self.df_coverage = pd.DataFrame({"position": position, "ploidy_cov": ploidy_coverage})
        self.bin_size = bin_size
        self.as_dev = as_dev
        self.logger = logger
        pass

    @staticmethod
    def get_coverage_slope(coverage, value_offset):
        """
        Calculate the slope of a given coverage array for n + value_offset positions.
        :param coverage:
        :param value_offset:
        :return: numpy array with the slope values
        """
        d = np.zeros(len(coverage))
        n = value_offset
        slope = 0
        # for index,y_0 in enumerate(coverage,n):
        for index in range(0, len(coverage), n):
            y_0 = coverage[index]
            y_1 = coverage[index:index + n + 1]
            if len(y_1) > 0:
                slope = (y_0 - np.mean(y_1)) / 1
            d[index:index + n] = slope
        return d

    def cnv_fit(self, final_cnvs):
        new_final_cnvs = []
        for final_cnv in final_cnvs:
            df_extended_cnv, new_start, new_end = self.__cnv_grow(final_cnv)
            new_start_trim, new_end_trim = self.__cnv_trim(df_extended_cnv, final_cnv, new_start, new_end)
            new_final_cnv = self.__cnv_overwrite(df_extended_cnv, final_cnv, new_start_trim, new_end_trim)
            new_final_cnvs.append(new_final_cnv)
        return new_final_cnvs

    def __cnv_grow(self, final_cnv):
        """
        Extends the coverage values of a CNV to the left and right up to bigger value of either 100*bin_size or
        25% of the CNV size.

        :param final_cnv:
        :return:
        """
        binsize = self.bin_size
        n = max([100 * binsize, final_cnv.size * 0.25])

        # get the coverage for the CNV
        df_cnv = self.df_coverage.loc[
                 (self.df_coverage["position"] >= final_cnv.start) & (self.df_coverage["position"] <= final_cnv.end), :]
        # get the coverage for the 100 extra bases before and after the CNV
        df_before = self.df_coverage.loc[
                    (self.df_coverage["position"] >= final_cnv.start - n) & (
                            self.df_coverage["position"] < final_cnv.start), :]
        df_after = self.df_coverage.loc[
                   (self.df_coverage["position"] > final_cnv.end) & (self.df_coverage["position"] <= final_cnv.end + n),
                   :]
        df_extended_cnv_coverage = pd.concat((df_before, df_cnv, df_after))

        # Calculate the slope
        df_extended_cnv_coverage["slope_1"] = self.get_coverage_slope(df_extended_cnv_coverage["ploidy_cov"].values, 1)
        df_extended_cnv_coverage["slope_5"] = self.get_coverage_slope(df_extended_cnv_coverage["ploidy_cov"].values, 5)
        df_extended_cnv_coverage["slope_10"] = self.get_coverage_slope(df_extended_cnv_coverage["ploidy_cov"].values,
                                                                       10)

        cnv_sd = np.std(df_cnv['ploidy_cov'].values) * 3
        cnv_mean = np.mean(df_cnv['ploidy_cov'].values)

        # determine if outside of cnv sd
        df_extended_cnv_coverage["is_in_sd"] = df_extended_cnv_coverage["ploidy_cov"].apply(
            lambda x: 1 if cnv_mean - cnv_sd < x < cnv_mean + cnv_sd else 0)

        upper_slope10_limit = df_extended_cnv_coverage["slope_10"].mean() + df_extended_cnv_coverage["slope_10"].std()
        lower_slope10_limit = df_extended_cnv_coverage["slope_10"].mean() - df_extended_cnv_coverage["slope_10"].std()

        # get new borders for CNV where is_in_sd 0 outside the position of final_cnv.start and final_cnv.end and
        # where the slope starts to be outside of the upper and lower limits get the first position where is_in_sd 0
        # after the final_cnv.end position
        new_end_list = df_extended_cnv_coverage.loc[
            (df_extended_cnv_coverage["is_in_sd"] == 0) & (df_extended_cnv_coverage["position"] > final_cnv.end) & (
                    lower_slope10_limit < df_extended_cnv_coverage["slope_10"]) & (
                    df_extended_cnv_coverage["slope_10"] < upper_slope10_limit), "position"]
        new_end = self.__get_end(new_end_list, final_cnv.end)

        # get the fist position where is_in_sd 0 before the final_cnv.start position
        new_start_list = df_extended_cnv_coverage.loc[
            (df_extended_cnv_coverage["is_in_sd"] == 0) & (df_extended_cnv_coverage["position"] < final_cnv.start) & (
                    lower_slope10_limit < df_extended_cnv_coverage["slope_10"]) & (
                    df_extended_cnv_coverage["slope_10"] < upper_slope10_limit), "position"]

        new_start = self.__get_start(new_start_list, final_cnv.start)

        df_extended = df_extended_cnv_coverage.loc[(df_extended_cnv_coverage["position"] >= new_start) & (
                df_extended_cnv_coverage["position"] <= new_end), :]
        return df_extended, new_start, new_end

    @staticmethod
    def __get_start(position_list, current_position):
        if len(position_list) > 0:
            return position_list.values[-1]
        return current_position

    @staticmethod
    def __get_end(position_list, current_position):
        if len(position_list) > 0:
            return position_list.values[0]
        return current_position

    def __cnv_trim(self, df_extended, final_cnv, new_start, new_end):
        """
        Trims the extended CNV to the left and right up to the first position where the slope is the biggest
        in the regions between the new borders and the final_cnv.start and final_cnv.end.
        :param df_extended:
        :param final_cnv:
        :param new_start:
        :param new_end:
        :return:
        """
        df_extended["slope_diff"] = abs(df_extended["slope_5"].diff())

        # get position of first biggest value of slope_diff after final_cnv.end and new_end
        end_max_positions = df_extended.loc[(df_extended["position"] > final_cnv.end) &
                                            (df_extended["position"] < new_end)]
        new_end_trim_list = end_max_positions.loc[
            end_max_positions["slope_diff"] == end_max_positions["slope_diff"].max(), "position"]

        new_end_trim = self.__get_end(new_end_trim_list, new_end)

        start_max_positions = df_extended.loc[(df_extended["position"] < final_cnv.start) &
                                              (df_extended["position"] > new_start)]

        new_start_trim_list = start_max_positions.loc[
            start_max_positions["slope_diff"] == start_max_positions["slope_diff"].max(), "position"]
        new_start_trim = self.__get_start(new_start_trim_list, new_start)

        return new_start_trim, new_end_trim

    def __cnv_overwrite(self, df_extended, final_cnv, new_start_trim, new_end_trim):
        """
        Overwrites the final_cnv with the new start and end positions and the new coverage.
        :param df_extended:
        :param final_cnv:
        :param new_start_trim:
        :param new_end_trim:
        :return:
        """
        # get coverage for the new CNV with the new start and end positions
        df_cnv = \
        df_extended.loc[(df_extended["position"] >= new_start_trim) & (df_extended["position"] <= new_end_trim), :][
            ["position", "ploidy_cov"]]

        ploidy_cov = df_cnv["ploidy_cov"].tolist()
        position = df_cnv["position"].tolist()
        start = int(position[0])
        end = int(position[-1])
        size = abs(end - start)
        final_cnv.start = start
        final_cnv.end = end
        final_cnv.size = size
        final_cnv.cov = ploidy_cov
        final_cnv.pos= position
        return final_cnv
