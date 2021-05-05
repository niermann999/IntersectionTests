#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ===============================================================================
# author: joana.niermann@cern.ch
# Plot a crude speedup overview comparing float vs double and avx vs sse
#
# Input needs to be given manually
# 
# ===============================================================================

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path


if __name__ == "__main__":
    plot_file_name = "./aggregate_plots/speedup_comparision.png"
    file_path  = Path("./aggregate_plots/")
    if not file_path.exists():
        os.makedirs("./aggregate_plots/")

    # medians
    eigen_median_avx_f  = 1.
    eigen_median_avx_d  = 1.
    eigen_median_sse_f  = 1.
    eigen_median_sse_d  = 1.

    vert_median_avx_f   = 0.90244178
    vert_median_avx_d   = 1.04607321
    vert_median_sse_f   = 1.02844139
    vert_median_sse_d   = 1.03390465

    hybrid_median_avx_f = 1.49827263
    hybrid_median_avx_d = 1.3589161
    hybrid_median_sse_f = 0.47665337
    hybrid_median_sse_d = 0.37315663

    horiz_median_avx_f  = 2.62233525
    horiz_median_avx_d  = 1.67191108
    horiz_median_sse_f  = 0.66469471
    horiz_median_sse_d  = 0.39744582

    # std deviations
    eigen_stddev_avx_f  = 0
    eigen_stddev_avx_d  = 0
    eigen_stddev_sse_f  = 0
    eigen_stddev_sse_d  = 0

    vert_stddev_avx_f   = 0
    vert_stddev_avx_d   = 0
    vert_stddev_sse_f   = 0
    vert_stddev_sse_d   = 0

    hybrid_stddev_avx_f = 0
    hybrid_stddev_avx_d = 0
    hybrid_stddev_sse_f = 0
    hybrid_stddev_sse_d = 0

    horiz_stddev_avx_f  = 0
    horiz_stddev_avx_d  = 0
    horiz_stddev_sse_f  = 0
    horiz_stddev_sse_d  = 0

    labels = ['sse double', 'sse float', 'axv double', 'avx float']
    eigen = np.array([eigen_median_sse_d,eigen_median_sse_f,eigen_median_avx_d,eigen_median_avx_f])
    eigen = np.round(eigen/eigen_median_sse_d, 1)
    eigen_err = np.round(np.array([eigen_stddev_sse_d,eigen_stddev_sse_f,eigen_stddev_avx_d,eigen_stddev_avx_f]), 1)

    AoS = np.array([vert_median_sse_d,vert_median_sse_f,vert_median_avx_d,vert_median_avx_f])
    AoS = np.round(AoS/eigen_median_sse_d, 1)
    AoS_err = np.round(np.array([vert_stddev_sse_d,vert_stddev_sse_f,vert_stddev_avx_d,vert_stddev_avx_f]), 1)

    hybrid = np.array([hybrid_median_sse_d,hybrid_median_sse_f,hybrid_median_avx_d,hybrid_median_avx_f])
    hybrid = np.round(hybrid/eigen_median_sse_d, 1)
    hybrid_err = np.round(np.array([hybrid_stddev_sse_d,hybrid_stddev_sse_f,hybrid_stddev_avx_d,hybrid_stddev_avx_f]), 1)

    horiz = np.array([horiz_median_sse_d,horiz_median_sse_f,horiz_median_avx_d,horiz_median_avx_f])
    horiz = np.round(horiz/eigen_median_sse_d, 1)
    horiz_err = np.round(np.array([horiz_stddev_sse_d,horiz_stddev_sse_f,horiz_stddev_avx_d,horiz_stddev_avx_f]), 1)


    x = np.arange(len(labels))  # the label locations
    width = 0.15  # the width of the bars

    fig, ax = plt.subplots()
    bars_eigen  = ax.bar(x - 1.5*width, eigen, width, label='eigen', yerr=eigen_err)
    bars_AoS    = ax.bar(x - 0.5*width, AoS, width, label='AoS', yerr=AoS_err)
    bars_hybrid = ax.bar(x + 0.5*width, hybrid, width, label='hybrid', yerr=hybrid_err)
    bars_horiz  = ax.bar(x + 1.5*width, horiz, width, label='horiz', yerr=horiz_err)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('speedup')
    ax.set_title('speed up relative to Eigen sse(double) - preliminary')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    #autolabel(bars_eigen, ax)
    #autolabel(bars_AoS, ax)
    #autolabel(bars_hybrid, ax)
    #autolabel(bars_horiz, ax)

    fig.tight_layout()

    plot_file_name = "./aggregate_plots/speedup_comparision.png"
    fig.savefig(plot_file_name, dpi=100)
    fig.clf()
    plt.close(fig)

# ==============================================================================
