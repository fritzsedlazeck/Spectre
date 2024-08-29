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
            all_res[filename] = pickle.load(f)
    return all_res


if __name__ == "__main__":
    cnvplot = CNVPlot()
    dpi = 300
    # num_cols = 1
    num_cols = 23
    scl = 8

    sample_id = "22-726"
    base_dir = f"/Users/kylematoba/Spectre/{sample_id}"

    plot_data_dir = os.path.join(base_dir, "plot_data")
    file_list = sorted(os.listdir(plot_data_dir))

    # file_list = ['chr1.pkl']
    all_res = load_all_res()
    all_y = [v["coverage_y"] for k, v in all_res.items()]
    grand_mean = np.nanmean(np.concat(all_y))

    num_files = len(file_list)

    num_rows = math.ceil(num_files / num_cols)
    # fig, axs = plt.subplots(num_rows,
    #                         num_cols,
    #                         figsize=(num_cols * scl, num_rows * scl / 2),
    #                         squeeze=False)

    fig = plt.figure(figsize=(num_cols * scl, num_rows * scl / 2))
    # gs = gridspec.GridSpec(num_rows, num_cols, hspace=0.6, wspace=0.3)
    gs = gridspec.GridSpec(num_rows, num_cols)

    # fullfilename1 = "/Users/kylematoba/Spectre/22-726/cnv_22-726_byCoverage_chrchr1.csv"
    # # fullfilename2 =
    #
    # file1 = pd.read_csv(fullfilename1)

    for row_idx in range(num_rows):
        for col_idx in range(num_cols):
            flat_idx = row_idx * num_cols + col_idx
            print(flat_idx)
            if flat_idx < num_files:
                data = all_res[file_list[flat_idx]]

                # gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
                # gss = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[flat_idx],  height_ratios=[5, 1],
                #                                        hspace=.2)
                gss = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[flat_idx],  height_ratios=[5, 1])
                ax0 = fig.add_subplot(gss[0])
                # ax1 = fig.add_subplot(gss[1], sharex=ax0)
                ax1 = fig.add_subplot(gss[1])

                ax0.set_title('{} {}'.format(sample_id, data["current_chromosome"]))
                # ax0.set_ylim(0, 6)
                ax0.set_ylim(0 - grand_mean, 6 - grand_mean)

                # main_plot = plt.subplot(gs[0])
                # candidates_plot = plt.subplot(gs[1])
                # candidates_plot.axes.get_yaxis().set_visible(False)

                # axs[row_idx, col_idx].plot(data["coverage_x"], data["coverage_y"])
                ax0.plot(data["coverage_x"], data["coverage_y"] - grand_mean, linewidth=.4, color="#67a9cf")
                ax0.axhline(0.0, color="k", linewidth=.8)

                for idx, cnv in enumerate(data["cand_tuples"]):
                    cnv_color = cnvplot.cnv_color[cnv[2]]
                    ax1.plot(np.array([cnv[0], cnv[1]]),
                                              np.array([0, 0]),
                                              linewidth='5',
                                              color=cnv_color)
                ax1.get_yaxis().set_visible(False)

    fullfilename_out = os.path.join(base_dir, "trungplot.png")
    fig.tight_layout()
    fig.savefig(fullfilename_out, dpi=dpi)
    # fig.savefig(fullfilename_out)
    print("Done")
