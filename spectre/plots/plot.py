import logging as logger
import matplotlib.pyplot as plot_engine
from matplotlib import gridspec
import numpy as np


class CoveragePlot:
    def __init__(self, as_dev=False):
        logger.basicConfig(level=logger.DEBUG) if as_dev else logger.basicConfig(level=logger.INFO)
        self.logger = logger
        # the plot
        self.figure = plot_engine.figure(figsize=(8, 4))
        gs = gridspec.GridSpec(1, 1)
        self.main_plot = plot_engine.subplot(gs[0])        # colors
        self.coverage_color = "#67a9cf"
        # legends
        self.main = ""
        self.x_axis = ""
        self.y_axis = ""
        self.file_prefix = "test"
        self.output_directory = "./"

    def plot_coverage(self, current_chromosome="", coverage=None):
        self.logger.debug("plotting coverage")
        # main plot
        self.main_plot.plot(np.array(coverage["pos"]), np.array(coverage["cov"]), color=self.coverage_color,
                            linewidth='0.5')
        # save and close
        self.figure.savefig(f'{self.output_directory}/{self.file_prefix}-{current_chromosome}.png', dpi=300)
        self.logger.info(f'Plot saved: {self.file_prefix}-{current_chromosome}.png')
        self.figure.clf()


class CNVPlot:
    def __init__(self, as_dev=False):
        logger.basicConfig(level=logger.DEBUG) if as_dev else logger.basicConfig(level=logger.INFO)
        self.logger = logger
        # the plot
        self.figure = plot_engine.figure(figsize=(8, 4))
        gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
        self.main_plot = plot_engine.subplot(gs[0])
        self.candidates_plot = plot_engine.subplot(gs[1])
        self.candidates_plot.axes.get_yaxis().set_visible(False)
        # colors
        self.coverage_color = "#67a9cf"
        self.cnv_color = {"DUP": "#d73027", "DEL": "#1a9850"}
        # legends
        self.main = ""
        self.x_axis = ""
        self.y_axis = ""
        self.axis_ylim = {"bottom": 0, "top": 6}  # not showing over 6x coverage, min can not be lower than 0
        self.file_prefix = "test"
        self.output_directory = "./"

    def plot_coverage_cnv(self, current_chromosome="", stats=None, coverage=None, cnv_cand_list=None, bounds=None):
        if stats is None or coverage is None:
            self.logger.error("bot parameters are needed")
        self.logger.debug("plotting coverage + CNV")
        # main plot
        self.main_plot.plot(np.array(coverage["pos"]), np.array(coverage["cov"]), color=self.coverage_color,
                            linewidth='0.5')
        self.main_plot.axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
        start = coverage["pos"][0]
        end = coverage["pos"][-1]
        self.candidates_plot.plot(np.array([start, end]), np.array([0, 0]), linewidth='0', color="#ffffff")
        # add CNV candidates
        self.logger.info(f"CNVs in chromosome: {current_chromosome}")
        if cnv_cand_list is not None:
            for cnv in cnv_cand_list:
                start = cnv.start
                end = cnv.end
                cnv_color = self.cnv_color[cnv.type]
                self.candidates_plot.plot(np.array([start, end]), np.array([0, 0]), linewidth='5', color=cnv_color)
        # save and close
        self.main_plot.plot(np.array([1, stats.chromosome_len]), np.array([stats.median, stats.median]),
                            linewidth='1', color="#000000")
        if bounds is not None:
            [upperb, lowerb] = bounds if len(bounds) == 2 else [np.NaN, np.NaN]
            self.main_plot.plot(np.array([1, stats.chromosome_len]), np.array([lowerb, lowerb]),
                                linewidth='1', color="#dd3497")
            self.main_plot.plot(np.array([1, stats.chromosome_len]), np.array([upperb, upperb]),
                                linewidth='1', color="#dd3497")
        self.figure.suptitle(f'{self.file_prefix} chromosome: {current_chromosome}')
        self.figure.savefig(f'{self.output_directory}/plot-coverage-{self.file_prefix}-{current_chromosome}.png', dpi=300)
        self.logger.info(f'Plot saved: plot-coverage-{self.file_prefix}-{current_chromosome}.png')
        self.figure.clf()
