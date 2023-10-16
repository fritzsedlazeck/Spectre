import numpy as np


class NormaldataAnalyser:

    @classmethod
    def get_candidate_statistics(cls, normalized_candidates) -> tuple:
        """
        :param normalized_candidates:
        :return:
        """
        min_val = np.nanmin(normalized_candidates)
        max_val = np.nanmax(normalized_candidates)
        avg = np.nanmean(normalized_candidates)
        std = np.nanstd(normalized_candidates)
        med = np.nanmedian(normalized_candidates)
        return avg, std, min_val, max_val, med

    @classmethod
    def normalize_candidates(cls, candidates: np.array(list), use_value: float) -> np.array(list):
        """
        Normalizes an i
        :param candidates:
        :param use_value: can be median or mean
        :return:
        """
        # return candidates/use_value
        return candidates/use_value  # the median will be CN2, we derive the rest from them



    @classmethod
    def get_slope(cls, n:int, normalized_candidates: np.array(list)) -> np.array(list):
        x = np.array(range(0, n))
        y = np.array(list)
        slope_array = np.zeros(normalized_candidates.size)
        l = 0
        slope = 0
        for i in range(0,normalized_candidates.size,1):
            if float(normalized_candidates[i]) != 0:
                l = i - n
                if i > n:
                    y = normalized_candidates[l:i].astype(np.float)
                    slope = cls.get_slope_from_values(x,y)
                else:
                    slope = 0
                slope_array[i] = float(slope)
        return slope_array


    @classmethod
    def get_slope_from_values(cls, x, y) -> float:
        """
        Calculates the slope of the regression of the given x and y values.
        :param x: Numpy array with walues in the X-axis
        :param y: Numpy array with values in the Y-axis
        :return: float of the slope value.
        """
        return ((x*y).mean(axis=0) - x.mean()*y.mean(axis=0)) / ((x**2).mean() - (x.mean())**2)
