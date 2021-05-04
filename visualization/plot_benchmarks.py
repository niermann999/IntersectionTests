#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ===============================================================================
# author: joana.niermann@cern.ch
# Visualization of the vectorization benchmarks
#
# 0) Warning: Can produce A LOT of files, depending on the number of data 
#    that was collected!
#
# 1) Create json output by running the benchmark executable:
#   ./src/BenchmarkIntersectionBench 
#        --benchmark_out=../visualization/Eigen_vs_Vc_bench.json        
#        --benchmark_out_format=json
#
# 2) Run the script in the directory where data file lies:
#   'python plot_benchmarks.py'
# 
# Each intersector type is benchmarked over a range of number of surfaces to 
# intersect. These are then repeated a number of times for better static certanty
# and to monitor fluctuations accross identical benchmarks. A mean, median and
# the standard deviation for each of those benchmarks are automatically 
# calculated by google benchmark, too.
#
# In this file the realtime of the benchmarks and their speedup compared to the 
# Eigen4F implementation are plotted.
# 
# ===============================================================================


import json
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
import pandas as pd


if __name__ == "__main__":

    if not os.path.exists('./benchmark_plots/'):
        os.makedirs('./benchmark_plots/')
    if not os.path.exists('./aggregate_plots/'):
        os.makedirs('./aggregate_plots/')

    matplotlib.rcParams.update({'font.size': 8})

    #
    # Get input data
    #
    print("Preparing input data...")

    input_file = json.load(open('./Eigen_vs_Vc_bench.json'))
    data = pd.DataFrame(input_file["benchmarks"])

    #prepare data for plots

    # Find max exec time to scale plots
    #max_time = data.loc[data['real_time'].idxmax()]['real_time']
    max_time = data.real_time.max()

    # Extract the number of intersected surfaces from benchmark label
    def extract_nSurfaces(label):
        sub_str = str.split(label, '/')
        return sub_str[1]

    data['nSurfaces'] = data['name'].apply(extract_nSurfaces)

    # Extract the intersector method
    def extract_intr_method(label):
        sub_str = str.split(label, '/')
        return sub_str[0] + '-' + sub_str[1]

    data['name'] = data['name'].apply(extract_intr_method)


    #
    # Outlier correction (TODO)
    #

    #
    # Data visualization plots
    #

    # Build dictionary of plots that need to be made
    plot_dict = {}
    for name in data['name']:
        name, n_surf = str.split(name, "-")
        if name not in plot_dict:
            plot_dict[name] = set()
        plot_dict[name].add(n_surf)

    # Plot function for the progression of exec time with repetition and the
    # histogram of exec times
    def create_data_plots(data_fr, method, n_surf):
        width = 0.6
        fig, axes = plt.subplots(nrows=2, ncols=1)
        fig.subplots_adjust(hspace=.3)

        # Timing over repetitions
        data_fr[(data_fr.run_type == 'iteration')].plot.bar(ax = axes[0], x='repetition_index', y='real_time', label="Exec. time",
                      width=width)
        # Settings plot #1
        axes[0].set_xlabel('#repetitions')
        axes[0].set_ylabel('t [ns]')
        xtick_labels = np.array(data_fr.repetition_index[(data_fr.run_type == 'iteration')])
        axes[0].set_xticklabels(labels=xtick_labels.astype(int), fontsize=7, rotation=90)
        axes[0].set_title(method+' - '+n_surf+" surface intersections")
        for label in axes[0].get_xticklabels()[::2]:
            label.set_visible(False)
        axes[0].legend(bbox_to_anchor=(1.05, 1.15))

        # Distribution of timings
        data_fr[(data_fr.run_type == 'iteration')].real_time.plot.hist(ax = axes[1], label="Time hist")

        # Settings plot #2
        axes[1].set_xlabel('t [ns]')

        # histogram legend
        handles, labels = axes[1].get_legend_handles_labels()
        # add figures of merit
        legend_str = ('mean: ' +  
                      str(np.around(data_fr.loc[data_fr['aggregate_name'] == 'mean']['real_time'].values[0], 2)) + '\n' + 
                      'median: ' + 
                      str(np.around(data_fr.loc[data_fr['aggregate_name'] == 'median']['real_time'].values[0], 2)) + '\n' + 
                      'stddev: ' + 
                      str(np.around(data_fr.loc[data_fr['aggregate_name'] == 'stddev']['real_time'].values[0], 2)))
        patch = mpatches.Patch(alpha = 0.0, color='grey', label=legend_str)
        handles.append(patch)
        axes[1].legend(handles=handles)

        plot_file_name = "./benchmark_plots/BenchmarkData_"+name+".png"
        fig.savefig(plot_file_name, dpi=200)
        plt.close()

    # Go through each method and number of surfaces and plot the timing for the repetitions, as well as the timing distribution
    print("Plotting benchmarks...")
    for method, surfs in plot_dict.items():
        for n_surf in surfs:
            name = method + "-" + n_surf
            plot_data = data[(data.name == name) & (data.nSurfaces == n_surf)]
            create_data_plots(plot_data, method, n_surf)


    #
    # Plot figures of merit and speedup
    #

    # Attach a text label above each bar in *rects*, displaying its height.
    """def autolabel(rects, axis):
        for rect in rects:
            height = rect.get_height()
            axis.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')"""

    print("Plotting aggregate data...")

    # Keep the mean, median, stddev
    data_dict = {}
    for method, surfs in plot_dict.items():
        data_dict[method + "_mean"] = []
        data_dict[method + "_median"] = []
        data_dict[method + "_stddev"] = []
        data_dict[method + "_surfs"] = []
        for n_surf in surfs:
            name = method + "-" + n_surf
            data_fr = data[(data.name == name) & (data.nSurfaces == n_surf)]

            mean = np.around(data_fr.loc[data_fr['aggregate_name'] == 'mean']['real_time'].values[0], 2)

            median = np.around(data_fr.loc[data_fr['aggregate_name'] == 'median']['real_time'].values[0], 2)

            stddev = np.around(data_fr.loc[data_fr['aggregate_name'] == 'stddev']['real_time'].values[0], 2)
            
            data_dict[method + "_mean"].append(mean)
            data_dict[method + "_median"].append(median)
            data_dict[method + "_stddev"].append(stddev)
            data_dict[method + "_surfs"].append(n_surf)
            
    data_frame = pd.DataFrame.from_dict(data_dict)
    data_frame = data_frame.assign(Eigen4D_surfs=data_frame.Eigen4D_surfs.astype(int)).sort_values(by=['Eigen4D_surfs'], ascending=True, inplace=False)

    #width = 0.6
    fig, axes = plt.subplots(nrows=2, ncols=1)
    fig.subplots_adjust(hspace=.3)
    data_frame.plot(ax=axes[0], x='Eigen4D_surfs', y=['Eigen4D_median','VcVert_median','VcHybrid_median','VcHoriz_median'], label=['Eigen4D','VcVert','VcHybrid','VcHoriz'], sharex=True)
    axes[0].set_ylabel('t [ns]')
    axes[0].set_title("Median of exec. Time vs. no. Surfaces")
    axes[0].legend()

    data_frame.plot(ax=axes[1], x='Eigen4D_surfs', y=['Eigen4D_wres_median','VcVert_wres_median','VcHybrid_wres_median','VcHoriz_wres_median'], label=['Eigen4D','VcVert','VcHybrid','VcHoriz'], sharex=True)
    axes[1].set_xlabel('no. surfaces')
    axes[1].set_ylabel('t [ns]')
    axes[1].set_title("Median of exec. Time vs. no. Surfaces (with vector container)")
    axes[1].legend()

    plot_file_name = "./aggregate_plots/Benchmark_scaling.png"
    fig.savefig(plot_file_name, dpi=200)
    plt.close()


    # Plot speedup
    data_frame[['Eigen4D_median','VcVert_median','VcHybrid_median','VcHoriz_median']] = data_frame[['Eigen4D_median','VcVert_median','VcHybrid_median','VcHoriz_median']].div(data_frame['Eigen4D_median'].values,axis=0)

    data_frame[['Eigen4D_wres_median','VcVert_wres_median','VcHybrid_wres_median','VcHoriz_wres_median']] = data_frame[['Eigen4D_wres_median','VcVert_wres_median','VcHybrid_wres_median','VcHoriz_wres_median']].div(data_frame['Eigen4D_wres_median'].values,axis=0)

    data_frame['VcVert_median'] = np.reciprocal(data_frame['VcVert_median'])
    data_frame['VcHybrid_median'] = np.reciprocal(data_frame['VcHybrid_median'])
    data_frame['VcHoriz_median'] = np.reciprocal(data_frame['VcHoriz_median'])
    data_frame['VcVert_wres_median'] = np.reciprocal(data_frame['VcVert_wres_median'])
    data_frame['VcHybrid_wres_median'] = np.reciprocal(data_frame['VcHybrid_wres_median'])
    data_frame['VcHoriz_wres_median'] = np.reciprocal(data_frame['VcHoriz_wres_median'])

    fig, axes = plt.subplots(nrows=2, ncols=1)
    fig.subplots_adjust(hspace=.3)
    data_frame.plot(ax=axes[0], x='Eigen4D_surfs', y=['Eigen4D_median','VcVert_median','VcHybrid_median','VcHoriz_median'], label=['Eigen4D','VcVert','VcHybrid','VcHoriz'], sharex=True)
    axes[0].set_ylabel('t [ns]')
    axes[0].set_title("Speedup vs. no. Surfaces")
    axes[0].legend()

    data_frame.plot(ax=axes[1], x='Eigen4D_surfs', y=['Eigen4D_wres_median','VcVert_wres_median','VcHybrid_wres_median','VcHoriz_wres_median'], label=['Eigen4D','VcVert','VcHybrid','VcHoriz'], sharex=True)
    axes[1].set_xlabel('no. surfaces')
    axes[1].set_ylabel('t [ns]')
    axes[1].set_title("Speedup vs. no. Surfaces (with vector container)")
    axes[1].legend()

    plot_file_name = "./aggregate_plots/Benchmark_speedup.png"
    fig.savefig(plot_file_name, dpi=200)
    plt.close()


# ==============================================================================
