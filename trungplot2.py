import os
import pickle
import math

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

from spectre.plots.plot import CNVPlot

# https://stackoverflow.com/questions/46683735/python-subplots-with-shared-axis-loop


def load_all_res() -> dict:
    all_res = dict()

    for filename in file_list:
        fullfilename = os.path.join(plot_data_dir, filename)
        # print(fullfilename)
        with open(fullfilename, 'rb') as f:
            all_res[filename.replace(".pkl", "")] = pickle.load(f)
    return all_res

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


if __name__ == "__main__":
    cnvplot = CNVPlot()
    # dpi = 300
    dpi = None
    scl = 8

    sample_id = "22-726"
    base_dir = f"/Users/kylematoba/Spectre/{sample_id}"

    plot_data_dir = os.path.join(base_dir, "plot_data")

    raw_filelist = sorted(set(os.listdir(plot_data_dir)) - {'chrX.pkl', 'chrY.pkl'})
    file_list = sorted(raw_filelist, key=lambda _: int(_.replace("chr", "").replace(".pkl", "")))

    # file_list = ['chr1.pkl']
    all_res = load_all_res()
    all_y_cat = np.concat([v["coverage_y"] for k, v in all_res.items()])
    num_bins = [v["coverage_x"].size for k, v in all_res.items()]

    breakpoints = np.cumsum(num_bins)

    grand_mean = np.nanmean(all_y_cat)
    num_files = len(file_list)

    xticks = []
    labels = [_.replace(".pkl", "") for _ in file_list]
    figsize = (num_files * scl, 1 * scl / 2)
    # figsize = (3 * scl, 1 * scl / 2)
    fig, ax = plt.subplots(1, 1, figsize=figsize, squeeze=True)
    # ax.plot(all_y_cat)
    to_plot_y = all_y_cat / grand_mean
    ax.semilogy(to_plot_y, ".", markersize=.4, base=2)

    breakpoints0 = np.array([0] + breakpoints.tolist())
    for idx, (k, v) in enumerate(all_res.items()):

        avg = np.nanmean(v["coverage_y"])
        xmin = breakpoints0[idx] / breakpoints0[-1]
        xmax = breakpoints0[idx + 1] / breakpoints0[-1]
        print(idx, k, avg, (xmin, xmax))
        ax.axhline(avg, xmin=xmin, xmax=xmax, color="r")

    to_plot_y_smoothed = pd.Series(to_plot_y).rolling(5000, min_periods=1000).mean().to_numpy()
    # ax.semilogy(to_plot_y_smoothed, base=2, color="r")
    # ax.plot(all_y_cat / grand_mean, ".", markersize=.4)

    plt.xticks(breakpoints, labels=labels)
    # plt.xticks(breakpoints, labels=labels)
    # for breakpoint in breakpoints:
    #     ax.vline(breakpoint)
    plt.xticks(rotation=-45)
    ax.set_xlim(0, len(all_y_cat))
    ax.grid()

    # for file_idx, filename in enumerate(file_list):
    #     print(file_idx, filename)

    # for row_idx in range(num_rows):
    #     for col_idx in range(num_cols):
    #         flat_idx = row_idx * num_cols + col_idx
    #         print(flat_idx)
    #         if flat_idx < num_files:
    #             data = all_res[file_list[flat_idx]]
    #
    #             # gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
    #             # gss = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[flat_idx],  height_ratios=[5, 1],
    #             #                                        hspace=.2)
    #             gss = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[flat_idx],  height_ratios=[5, 1])
    #             ax0 = fig.add_subplot(gss[0])
    #             # ax1 = fig.add_subplot(gss[1], sharex=ax0)
    #             ax1 = fig.add_subplot(gss[1])
    #
    #             ax0.set_title('{} {}'.format(sample_id, data["current_chromosome"]))
    #             # ax0.set_ylim(0, 6)
    #             ax0.set_ylim(0 - grand_mean, 6 - grand_mean)
    #
    #             # main_plot = plt.subplot(gs[0])
    #             # candidates_plot = plt.subplot(gs[1])
    #             # candidates_plot.axes.get_yaxis().set_visible(False)
    #
    #             # axs[row_idx, col_idx].plot(data["coverage_x"], data["coverage_y"])
    #             ax0.plot(data["coverage_x"], data["coverage_y"] - grand_mean, linewidth=.4, color="#67a9cf")
    #             ax0.axhline(0.0, color="k", linewidth=.8)
    #             for idx, cnv in enumerate(data["cand_tuples"]):
    #                 cnv_color = cnvplot.cnv_color[cnv[2]]
    #                 ax1.plot(np.array([cnv[0], cnv[1]]),
    #                                           np.array([0, 0]),
    #                                           linewidth='5',
    #                                           color=cnv_color)
    #             ax1.get_yaxis().set_visible(False)
    ax.axhline(1.0, color="k", linewidth=.8)
    fullfilename_out = os.path.join(base_dir, "trungplot2.png")
    fig.tight_layout()
    fig.savefig(fullfilename_out, dpi=dpi)
    # fig.savefig(fullfilename_out)
    print("Done")
