# !/usr/bin/env python3

# Author: burcakotlu

# Contact: burcakotlu@eng.ucsd.edu


"""

This python file plots a grid where rows are annotated regions, and columns are signatures.
Cells contain circles in which radius reflects fold changes and color reflects the -log10 qvalue.

"""


import multiprocessing
import math
import os
import pandas as pd
import numpy as np

from scipy import stats
import statsmodels.stats.multitest
from functools import reduce
from statsmodels.stats.weightstats import ztest

from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import colors

from SigProfilerTopography.source.annotatedregion.AnnotatedRegionAnalysis import read_lncRNA_GRCh37_annotation_file
from SigProfilerTopography.source.annotatedregion.AnnotatedRegionAnalysis import create_intervaltree

from SigProfilerTopography.source.commons.TopographyCommons import natural_key
from SigProfilerTopography.source.commons.TopographyCommons import GRCh37
from SigProfilerTopography.source.commons.TopographyCommons import GRCh38
from SigProfilerTopography.source.commons.TopographyCommons import Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import get_genome_length

def get_best_annotated_regions(combined_df,
                               signatures,
                               num_of_best_annotated_regions,
                               column_for_sorting,
                               column_for_sorting_ascending,
                               column_for_significance_level,
                               significance_level):

    print('combined_df.columns.values:', combined_df.columns.values)

    all_annotated_regions = set()
    selected_signatures = []

    for signature in signatures:
        sorting_column_name = '%s_%s' %(signature, column_for_sorting)
        significance_column_name = '%s_%s' % (signature, column_for_significance_level)
        if significance_column_name in combined_df.columns.values:
            df = combined_df[combined_df[significance_column_name] < significance_level].copy()
            df.sort_values(by=[sorting_column_name], ascending=column_for_sorting_ascending, inplace=True) # TODO check A value is trying to be set on a copy of a slice from a DataFrame
            annotated_regions = df['name'][0:num_of_best_annotated_regions].values
            qvalues = df[significance_column_name][0:num_of_best_annotated_regions].values
            all_annotated_regions = all_annotated_regions.union(annotated_regions)
            print('DEBUG', signature, 'annotated_regions:', annotated_regions)
            print('DEBUG', signature, 'qvalues:', qvalues)
            if len(annotated_regions) > 0:
                selected_signatures.append(signature)

    print('DEBUG len(all_annotated_regions):', len(all_annotated_regions))
    print('DEBUG all_annotated_regions:', all_annotated_regions)
    print('DEBUG selected_signatures:', selected_signatures)

    return selected_signatures, all_annotated_regions

def apply_multiple_testing(combined_df, signatures, pvalue_of_interest):
    print('combined_df.columns.values:', combined_df.columns.values)

    for signature in signatures:
        pvalue_column_name = '%s_%s' % (signature, pvalue_of_interest)
        qvalue_column_name = '%s_qvalue' % (signature)

        if pvalue_column_name in combined_df.columns.values:
            pvalue_array = combined_df[pvalue_column_name].values
            pvalue_array = pvalue_array.astype('float64')

            nan_mask_array = np.isnan(pvalue_array)
            pvalues_reduced_array = pvalue_array[~nan_mask_array]

            qvalues_array = np.zeros(pvalue_array.shape)
            qvalues_array[nan_mask_array] = np.nan

            if len(pvalues_reduced_array) > 0:
                rejected, qvalues_reduced_array, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(
                    pvalues_reduced_array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

                qvalues_array[~nan_mask_array] = qvalues_reduced_array

            combined_df[qvalue_column_name] = qvalues_array


def fill_arrays(combined_df,
                annotated_regions,
                signatures,
                column_for_circle_radius,
                column_for_circle_color,
                significance_level):

    # sort the annotated regions in ascending order
    annotated_regions_sorted_list = list(annotated_regions)
    annotated_regions_sorted_list.sort()

    # sort the signatures in ascending order
    signatures = sorted(signatures, key=natural_key)

    circle_radius_array = np.zeros((len(annotated_regions_sorted_list), len(signatures)), dtype=float)
    circle_color_array = np.zeros((len(annotated_regions_sorted_list), len(signatures)), dtype=float)

    for i, annotated_region in enumerate(annotated_regions_sorted_list):
        for j, signature in enumerate(signatures):
            circle_radius_column_name = '%s_%s' % (signature, column_for_circle_radius)
            circle_color_column_name = '%s_%s' % (signature, column_for_circle_color)

            if (circle_radius_column_name in combined_df.columns.values) and (circle_color_column_name in combined_df.columns.values):
                radius_value = combined_df[combined_df['name'] == annotated_region][circle_radius_column_name].values[0]
                color_value = combined_df[combined_df['name'] == annotated_region][circle_color_column_name].values[0]
                if color_value < significance_level:
                    circle_radius_array[i][j] = radius_value
                    circle_color_array[i][j] = color_value
                else:
                    circle_radius_array[i][j] = np.nan
                    circle_color_array[i][j] = np.nan
            else:
                circle_radius_array[i][j] = np.nan
                circle_color_array[i][j] = np.nan

    print('DEBUG len(annotated_regions_sorted_list):', len(annotated_regions_sorted_list))
    print('DEBUG annotated_regions_sorted_list:', annotated_regions_sorted_list)
    print('DEBUG len(signatures):', len(signatures))
    print('DEBUG signatures:', signatures)
    print('DEBUG circle_radius_array.shape:', circle_radius_array.shape)
    print('DEBUG circle_radius_array:', circle_radius_array)

    # which signatures are all nans?
    circle_radius_array_columns_all_nans = np.all(np.isnan(circle_radius_array), axis=0)

    # remove the columns all nans
    circle_radius_array = circle_radius_array[:, ~circle_radius_array_columns_all_nans]

    # remove the qvalue_array columns all nans
    circle_color_array = circle_color_array[:, ~circle_radius_array_columns_all_nans]

    print('DEBUG signatures:', signatures)
    print('DEBUG signatures_all_nans:', circle_radius_array_columns_all_nans)

    # Remove the signatures all with nans
    signatures_array = np.array(signatures)[~circle_radius_array_columns_all_nans]

    # replace nans with zero
    circle_radius_array[np.isnan(circle_radius_array)] = 0
    circle_color_array[np.isnan(circle_color_array)] = 0

    print('BEFORE RETURN')
    print('annotated_regions_sorted_list:', annotated_regions_sorted_list)
    print('signatures_array.shape:', signatures_array.shape)
    print('circle_radius_array.shape:', circle_radius_array.shape)
    print('circle_color_array.shape:', circle_color_array.shape)

    print('np.sum(np.isnan(circle_radius_array)):', np.sum(np.isnan(circle_radius_array)))
    print('np.sum(np.isnan(circle_color_array)):', np.sum(np.isnan(circle_color_array)))

    return annotated_regions_sorted_list, \
           signatures_array, \
           circle_radius_array, \
           circle_color_array

def get_radius(row_index, column_index, column_sorted_array, circle_radius_max, radius_adjustment_value):
    # row_max = np.max(column_sorted_array[row_index,:])
    # radius = (column_sorted_array[row_index][column_index] / (row_max + radius_adjustment_value)) / 2
    radius = (column_sorted_array[row_index][column_index] / (circle_radius_max + radius_adjustment_value)) / 2
    return radius


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def calculate_radius_add_patch_updated(row_index,
                                       column_index,
                                       circle_radius_array,
                                       circle_color_array,
                                       significance_level,
                                       radius_adjustment_value,
                                       cmap,
                                       colorbar_v_min,
                                       colorbar_v_max,
                                       circle_radius_max,
                                       ax):


    v_min = colorbar_v_min
    v_max = colorbar_v_max
    # Very important: You have to normalize
    norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)

    radius = get_radius(row_index,
                        column_index,
                        circle_radius_array,
                        circle_radius_max,
                        radius_adjustment_value)

    qvalue = circle_color_array[row_index][column_index]

    try:
        minus_log10_qvalue = -math.log10(qvalue)
    except ValueError as e:
        minus_log10_qvalue = 1000

    print('calculate_radius_add_patch_updated:', 'row_index:', row_index, 'column_index:', column_index,
          'qvalue:', qvalue, 'radius:', radius, 'minus_log10_qvalue:', minus_log10_qvalue)

    if ((radius > 0) and (qvalue < significance_level)):
        ax.add_patch(plt.Circle((column_index + 0.5, row_index + 0.5), radius, color=cmap(norm(minus_log10_qvalue)), fill=True, zorder=5))

def plot_radiusbar_vertical(output_path, circle_radius_min, circle_radius_max, annotated_regions_filename):

    diameter_labels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    if circle_radius_min < 1.1:
        circle_radius_min = 1.1

    circle_radius_middle = (circle_radius_min + circle_radius_max)/2
    diameter_ticklabels = [np.round(circle_radius_min,1), '', '', '', np.round(circle_radius_middle,1), '', '', '', '', np.round(circle_radius_max,1)]

    row_labels = ['circle']

    # Make a figure and axes with dimensions as desired.
    fig, ax = plt.subplots(figsize=(5, 8))

    ax.grid(which="major", color="black", linestyle='-') #  linewidth=2

    # Make the borders black
    plt.setp(ax.spines.values(), color='black')

    # make aspect ratio square
    ax.set_aspect(1.0)

    ax.set_facecolor('white')

    for row_index, row_label in enumerate(row_labels):
        for diameter_index, diameter_label in enumerate(diameter_labels):
            circle = plt.Circle((row_index + 0.5, diameter_index + 0.5),
                                radius=(diameter_label / (2 * 1.08)), color='gray', fill=True)
            ax.add_artist(circle)

    # CODE GOES HERE TO CENTER X-AXIS LABELS...
    ax.set_xlim([0, len(row_labels)])
    ax.set_xticklabels([])

    ax.tick_params(axis='x', which='minor', length=0, labelsize=12)
    # major ticks
    ax.set_xticks(np.arange(0, len(row_labels), 1))
    # minor ticks
    ax.set_xticks(np.arange(0, len(row_labels), 1) + 0.5, minor=True)

    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='major',  # both major and minor ticks are affected
        bottom=False)  # labels along the bottom edge are off


    # CODE GOES HERE TO CENTER Y-AXIS LABELS...
    ax.set_ylim([0, len(diameter_labels)])
    ax.set_yticklabels([])

    ax.tick_params(axis='y', which='minor', length=0, labelsize=25)
    # major ticks
    ax.set_yticks(np.arange(0, len(diameter_labels), 1))
    # minor ticks
    ax.set_yticks(np.arange(0, len(diameter_labels), 1) + 0.5, minor=True)
    ax.set_yticklabels(diameter_ticklabels, minor=True)

    ax.yaxis.set_ticks_position('right')

    plt.tick_params(
        axis='y',  # changes apply to the x-axis
        which='major',  # both major and minor ticks are affected
        left=False,  # ticks along the bottom edge are off
        right=False)  # labels along the bottom edge are off

    ax.yaxis.set_label_position("right")
    ax.set_ylabel('Fold\nchange', fontsize=30, labelpad=50, rotation=0)

    filename = '%s_radius_bar.png' %(annotated_regions_filename)
    figureFile = os.path.join(output_path, filename)

    fig.savefig(figureFile)
    plt.close()

def plot_colorbar_vertical(output_path, cmap, v_min, v_max, annotated_regions_filename):
    fig = plt.figure(figsize=(4, 10))
    ax = fig.add_axes([0.05, 0.05, 0.1, 0.9])

    # If a ListedColormap is used, the length of the bounds array must be
    # one greater than the length of the color list.  The bounds must be
    # monotonically increasing.

    bounds = np.arange(v_min, v_max + 1, 2)
    norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)

    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=bounds, spacing='proportional', orientation='vertical')
    cb.ax.tick_params(labelsize=30)
    cb.set_label("-log10\n(q-value)", verticalalignment='center', rotation=0, labelpad=100, fontsize=40)

    filename = '%s_color_bar_vertical.png' %(annotated_regions_filename)
    figureFile = os.path.join(output_path, filename)

    fig.savefig(figureFile)
    plt.close()


def plot_circles_on_grid(output_path,
                         rows_list,
                         columns_array,
                         circle_radius_array,
                         circle_color_array,
                         significance_level,
                         radius_adjustment_value,
                         cmap,
                         colorbar_v_min,
                         colorbar_v_max,
                         circle_radius_max,
                         annotated_regions_filename):

    # annotated_regions_sorted_list,
    # signatures_array,
    # circle_radius_array,
    # circle_color_array,
    # radius_adjustment_value,
    # colorbar_v_min,
    # colorbar_v_max

    print('plot_circles_on_grid')
    print('rows_list', rows_list)
    print('columns_array', columns_array)
    print('circle_radius_array.shape:', circle_radius_array.shape)
    print('circle_radius_array:', circle_radius_array)
    print('circle_color_array.shape:', circle_color_array.shape)
    print('circle_color_array:', circle_color_array)

    figwidth = None
    figheight = None
    font_size = None
    title_font_size = None
    label_font_size = None
    cell_font_size = None

    if len(rows_list) > 0 or len(columns_array) > 0:
        dpi = 100
        squaresize = 200  # pixels

        if figwidth is None:
            if len(columns_array) <= 3:
                figwidth = 2 * len(columns_array) * squaresize / float(dpi)
            else:
                figwidth = len(columns_array) * squaresize / float(dpi)
        if figheight is None:
            figheight = len(rows_list) * squaresize / float(dpi)

        # Set font size w.r.t. number of rows
        if font_size is None:
            font_size = 10 * np.sqrt(len(circle_radius_array))
        if title_font_size is None:
            title_font_size = font_size
        if label_font_size is None:
            label_font_size = font_size

        # Set cell font size w.r.t. number of columns
        if cell_font_size is None:
            cell_font_size = font_size * 2 / 3

        fig, ax = plt.subplots(1, figsize=(figwidth, figheight), dpi=dpi)
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

        ax.set_facecolor('white')

        # make aspect ratio square
        ax.set_aspect(1.0)

        # reverse such that at index 0 we have the former last element
        rows_list.reverse()

        # CODE GOES HERE TO CENTER X-AXIS LABELS...
        ax.set_xlim([0, len(columns_array)])
        ax.set_xticklabels([])
        ax.tick_params(axis='x', which='minor', labelsize=font_size,
                       bottom=True, top=True,
                       labelbottom=True, labeltop=True)  # labelsize=35 width=1 length=6

        ax.tick_params(axis='x', which='major', labelsize=font_size,
                       bottom=False, top=False,
                       labelbottom=False, labeltop=False)  # labelsize=35 width=1 length=6

        # major ticks
        ax.set_xticks(np.arange(0, len(columns_array), 1))
        # minor ticks
        ax.set_xticks(np.arange(0, len(columns_array), 1) + 0.5, minor=True)
        ax.set_xticklabels(columns_array, minor=True, fontweight='bold', fontname='Arial', rotation=90)

        # CODE GOES HERE TO CENTER Y-AXIS LABELS...
        ax.set_ylim([0, len(rows_list)])
        ax.set_yticklabels([])
        ax.tick_params(axis='y', which='minor', labelsize=font_size,
                       left=True, right=True,
                       labelleft=True, labelright=True)  # labelsize=35 width=1 length=6,

        ax.tick_params(axis='y', which='major', labelsize=font_size,
                       left=False, right=False,
                       labelleft=False, labelright=False)  # labelsize=35 width=1 length=6

        # major ticks
        ax.set_yticks(np.arange(0, len(rows_list), 1))
        # minor ticks
        ax.set_yticks(np.arange(0, len(rows_list), 1) + 0.5, minor=True)
        ax.set_yticklabels(rows_list, minor=True, fontname="Arial", fontweight='bold')  # fontsize

        # Gridlines based on major ticks
        ax.grid(which='major', color='black', zorder=1)
        ax.grid(which='minor', color='gray', linestyle='--', zorder=0)

        # reverse the rows of the numpy array such that at index 0 we have the former last element
        circle_radius_array = np.flip(circle_radius_array, axis=0)
        circle_color_array = np.flip(circle_color_array, axis=0)

        for row_index, row in enumerate(rows_list):
            for column_index, column in enumerate(columns_array):
                print('plot_circles_on_grid:', 'row_index:', row_index, 'column_index:', column_index,
                      'row:', row, 'column:', column)
                calculate_radius_add_patch_updated(row_index,
                                                   column_index,
                                                   circle_radius_array,
                                                   circle_color_array,
                                                   significance_level,
                                                   radius_adjustment_value,
                                                   cmap,
                                                   colorbar_v_min,
                                                   colorbar_v_max,
                                                   circle_radius_max,
                                                   ax)

        filename = '%s.png' %(annotated_regions_filename)
        figureFile = os.path.join(output_path, filename)
        fig.savefig(figureFile, dpi=dpi, bbox_inches="tight")
        plt.close()


def calculate_binomtest_fisherexacttest_ztest_pvalue(signature,
                                                     annotated_region_name,
                                                     real_num_of_hits,
                                                     sims_max_hits_at_annotated_region,
                                                     expected_num_of_hits,
                                                     sims_num_of_hits,
                                                     signature_num_of_mutations,
                                                     p,
                                                     contingency_table):

    # Binomial Test (with replacement)
    result = stats.binomtest(real_num_of_hits, signature_num_of_mutations, p, alternative='greater')
    binomtest_pvalue = result.pvalue

    # Fisher Exact Test
    oddsratio, fisherexacttest_pvalue = stats.fisher_exact(contingency_table)

    # ztest
    zstat1, ztest_pvalue = ztest(sims_num_of_hits, [real_num_of_hits], alternative='two-sided')

    # qvalue
    qvalue = np.nan

    # foldchange
    foldchange = np.nan

    if sims_max_hits_at_annotated_region != 0:
        foldchange = real_num_of_hits /sims_max_hits_at_annotated_region

    return (signature,
            annotated_region_name,
            real_num_of_hits,
            sims_max_hits_at_annotated_region,
            expected_num_of_hits,
            sims_num_of_hits,
            binomtest_pvalue,
            fisherexacttest_pvalue,
            ztest_pvalue,
            qvalue,
            foldchange)

def get_np_array_rows_runs_columns_ordered_annotated_regions(signature, signatures_df):
    arr_stacked = None
    signature_based_num_of_mutations = None

    if np.any(signatures_df[signatures_df['signature'] == signature]['number_of_mutations'].values):
        signature_based_num_of_mutations = \
        signatures_df[signatures_df['signature'] == signature]['number_of_mutations'].values[0]
        print(signature, 'signature_based_num_of_mutations:', signature_based_num_of_mutations)

        if signature == 'aggregatedsubstitutions':
            output_path = os.path.join(
                '/restricted/alexandrov-ddn/users/burcak/SigProfilerTopographyRuns/Mutographs_ESCC_552/topography/all_escc_552_samples/data/lncRNA/aggregatedsubstitutions')
        elif signature == 'aggregatedindels':
            output_path = os.path.join(
                '/restricted/alexandrov-ddn/users/burcak/SigProfilerTopographyRuns/Mutographs_ESCC_552/topography/all_escc_552_samples/data/lncRNA/aggregatedindels')
        else:
            output_path = os.path.join(
                '/restricted/alexandrov-ddn/users/burcak/SigProfilerTopographyRuns/Mutographs_ESCC_552/topography/all_escc_552_samples/data/lncRNA/signaturebased')

        arr_list = []
        for sim_num in range(101):
            if sim_num == 0:
                if signature != 'aggregatedsubstitutions' and signature != 'aggregatedindels':
                    file_name = '%s_AccumulatedSignalArray.txt' % (signature)
                else:
                    file_name = 'all_escc_552_samples_AccumulatedNumofHits.txt'

            else:
                if signature != 'aggregatedsubstitutions' and signature != 'aggregatedindels':
                    file_name = '%s_sim%s_AccumulatedSignalArray.txt' % (signature, sim_num)
                else:
                    file_name = 'all_escc_552_samples_sim%d_AccumulatedNumofHits.txt' % (sim_num)

            if os.path.exists(os.path.join(output_path, file_name)):
                arr = np.loadtxt(os.path.join(output_path, file_name), dtype=float, delimiter='\t')
                arr_list.append(arr)
                print('sim_num:', sim_num, 'arr.shape:', arr.shape, 'arr.size:', arr.size, 'arr.dtype:', arr.dtype,
                      'arr:', arr)

        if len(arr_list) > 0:
            arr_stacked = np.stack(arr_list, axis=0)
            print('arr_stacked.shape:', arr_stacked.shape)  # (101, 19974)

    return signature_based_num_of_mutations, arr_stacked


def convert_dict_2_df(signature, signature_num_of_mutations, signature_2_annotated_region_2_result_dict):
    df = None

    if (signature_num_of_mutations is not None) and (signature in signature_2_annotated_region_2_result_dict):
        # result -> list
        # [0] real_num_of_hits,
        # [1] sims_max_hits_at_annotated_region
        # [2] expected_num_of_hits,
        # [3] np.mean(sims_num_of_hits),
        # [4] np.std(sims_num_of_hits),
        # [5] sims_num_of_hits,
        # [6] binomtest_pvalue,
        # [7] fisherexacttest_pvalue,
        # [8] ztest_pvalue,
        # [9] qvalue
        # [10] foldchange

        L = sorted([(signature,
                     annotated_region_name,
                     result[0], # real_num_of_hits
                     result[1], # sims_max_hits_at_annotated_region
                     result[2], # expected_num_of_hits
                     result[3], # np.mean(sims_num_of_hits)
                     result[4], # np.std(sims_num_of_hits)
                     result[5], # sims_num_of_hits
                     result[6], # binomtest_pvalue
                     result[7], # fisherexacttest_pvalue
                     result[8], # ztest_pvalue
                     result[9], # qvalue
                     result[10] # foldchange
                     )
                    for signature, annotated_region_2_result_dict in signature_2_annotated_region_2_result_dict.items()
                    for annotated_region_name, result in annotated_region_2_result_dict.items()])

        df = pd.DataFrame(L, columns=['signature',
                                      'name',
                                      'real_num_of_hits',
                                      'sims_max_hits_at_annotated_region',
                                      'expected_num_of_hits',
                                      'sims_mean',
                                      'sims_std',
                                      'sims_num_of_hits',
                                      'binomtest_pvalue',
                                      'fisherexacttest_pvalue',
                                      'ztest_pvalue',
                                      'qvalue',
                                      'foldchange'])

    return df

def plot_highly_mutated_annotated_regions_figure(annotated_regions_filename,
                  genome_length,
                  signatures,
                  signature_type,
                  signatures_df,
                  ordered_annotated_regions_array,
                  annotated_region_2_length_dict,
                  run_parallel,
                  num_of_best_annotated_regions_per_signature,
                  pvalue_column_of_interest,
                  column_for_sorting,
                  column_for_sorting_ascending,
                  column_for_circle_radius,
                  column_for_circle_color,
                  column_for_significance_level,
                  significance_level,
                  radius_adjustment_value,
                  cmap,
                  colorbar_v_min,
                  colorbar_v_max,
                  output_path):

    df_list = []

    for signature in signatures:

        signature_2_annotated_region_2_result_dict = {}
        signature_num_of_mutations, signature_rows_runs_columns_regions_arr = get_np_array_rows_runs_columns_ordered_annotated_regions(signature,
                                                                                                                                       signatures_df)
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)

        def write_result(result_tuple):

            signature = result_tuple[0]
            annotated_region_name = result_tuple[1] # e.g.: LINC002237
            real_num_of_hits = result_tuple[2]
            sims_max_hits_at_annotated_region = result_tuple[3]
            expected_num_of_hits = result_tuple[4]
            sims_num_of_hits = result_tuple[5]
            binomtest_pvalue = result_tuple[6]
            fisherexacttest_pvalue = result_tuple[7]
            ztest_pvalue = result_tuple[8]
            qvalue = result_tuple[9]
            foldchange = result_tuple[10]

            # result -> list
            # [0] real_num_of_hits,
            # [1] expected_num_of_hits,
            # [2] np.mean(sims_num_of_hits),
            # [3] np.std(sims_num_of_hits),
            # [4] sims_num_of_hits,
            # [5] binomtest_pvalue,
            # [6] fisherexacttest_pvalue,
            # [7] ztest_pvalue,
            # [8] qvalue
            # [9] foldchange

            if (np.all(real_num_of_hits > sims_num_of_hits) and (foldchange is not None) and (foldchange > 1) and
                    (real_num_of_hits > expected_num_of_hits)):

                if signature in signature_2_annotated_region_2_result_dict:
                    signature_2_annotated_region_2_result_dict[signature][annotated_region_name] = [real_num_of_hits,
                                                                           sims_max_hits_at_annotated_region,
                                                                           expected_num_of_hits,
                                                                           np.mean(sims_num_of_hits),
                                                                           np.std(sims_num_of_hits),
                                                                           sims_num_of_hits,
                                                                           binomtest_pvalue,
                                                                           fisherexacttest_pvalue,
                                                                           ztest_pvalue,
                                                                           qvalue,
                                                                           foldchange]
                else:
                    signature_2_annotated_region_2_result_dict[signature] = {}
                    signature_2_annotated_region_2_result_dict[signature][annotated_region_name] = [real_num_of_hits,
                                                                           sims_max_hits_at_annotated_region,
                                                                           expected_num_of_hits,
                                                                           np.mean(sims_num_of_hits),
                                                                           np.std(sims_num_of_hits),
                                                                           sims_num_of_hits,
                                                                           binomtest_pvalue,
                                                                           fisherexacttest_pvalue,
                                                                           ztest_pvalue,
                                                                           qvalue,
                                                                           foldchange]

        # signature_rows_runs_columns_regions_arr:     rows: real and simulations topography runs
        # signature_rows_runs_columns_regions_arr:  columns: ordered annotated regions
        if (signature_num_of_mutations is not None) and (signature_rows_runs_columns_regions_arr is not None):
            # j refers to an annotated region
            for j in range(signature_rows_runs_columns_regions_arr.shape[1]):
                annotated_region_name = ordered_annotated_regions_array[j]
                num_of_real_sims_hits = signature_rows_runs_columns_regions_arr[:,j]
                real_num_of_hits = int(num_of_real_sims_hits[0]) # e.g. 4.0 -> 4
                sims_num_of_hits = num_of_real_sims_hits[1:]
                expected_num_of_hits = (annotated_region_2_length_dict[annotated_region_name] * signature_num_of_mutations) / genome_length

                # For Binom Test
                annotated_region_length = annotated_region_2_length_dict[annotated_region_name]
                p = annotated_region_length / genome_length  # TODO get the length of non coding regions

                # For Fisher Exact Test
                real_hits_at_annotated_region = real_num_of_hits
                sims_max_hits_at_annotated_region = math.ceil(max(p * signature_num_of_mutations, np.mean(sims_num_of_hits)))
                real_not_at_annotated_region = signature_num_of_mutations - real_hits_at_annotated_region
                sims_not_at_annotated_region = signature_num_of_mutations - sims_max_hits_at_annotated_region
                contingency_table = [[real_hits_at_annotated_region, sims_max_hits_at_annotated_region],
                                     [real_not_at_annotated_region, sims_not_at_annotated_region]]

                if run_parallel:
                    jobs = []
                    jobs.append(pool.apply_async(calculate_binomtest_fisherexacttest_ztest_pvalue,
                                                 args=(signature,
                                                       annotated_region_name,
                                                        real_num_of_hits,
                                                        sims_max_hits_at_annotated_region,
                                                        expected_num_of_hits,
                                                        sims_num_of_hits,
                                                        signature_num_of_mutations,
                                                        p,
                                                        contingency_table,),
                                                 callback=write_result))

                    # wait for all jobs to finish
                    for _ in jobs:
                        print(signature, annotated_region_name)

                else:
                    result_tuple = calculate_binomtest_fisherexacttest_ztest_pvalue(signature,
                                                                                    annotated_region_name,
                                                                                    real_num_of_hits,
                                                                                    sims_max_hits_at_annotated_region,
                                                                                    expected_num_of_hits,
                                                                                    sims_num_of_hits,
                                                                                    signature_num_of_mutations,
                                                                                    p,
                                                                                    contingency_table)

                    write_result(result_tuple)

        pool.close()
        pool.join()

        # convert dict to a dataframe
        df = convert_dict_2_df(signature, signature_num_of_mutations, signature_2_annotated_region_2_result_dict)

        if df is not None:
            df.columns = [column if column == 'signature' or column == 'name' else signature + '_' + column for column in df.columns]
            print(df.shape)

            # sort w.r.t column_to_be_sorted
            df.sort_values(by=['%s_%s' %(signature, column_for_sorting)], ascending=column_for_sorting_ascending, inplace=True)
            filename = '%s_%s.txt' %(signature, annotated_regions_filename)
            df.to_csv(os.path.join(output_path, filename), sep="\t", index=False, header=True)
            df.drop(['signature'], inplace=True, axis=1)
            df_list.append(df)

    combined_df = reduce(lambda df1, df2: df1.merge(df2, how='inner', on='name'), df_list)
    apply_multiple_testing(combined_df, signatures, pvalue_column_of_interest)
    filename = '%s_common_%s.txt' %(signature_type,annotated_regions_filename)
    combined_df.to_csv(os.path.join(output_path, filename), sep="\t",
              index=False, header=True)

    combined_df = reduce(lambda df1, df2: df1.merge(df2, how='outer', on='name'), df_list)
    apply_multiple_testing(combined_df, signatures, pvalue_column_of_interest)
    filename = '%s_all_%s.txt' %(signature_type,annotated_regions_filename)
    combined_df.to_csv(os.path.join(output_path, filename), sep="\t",
              index=False, header=True)

    # TODO DEFINE PROMOTER REGION
    # TODO what about copy number alterations? how to integrate them?
    # sort and select
    selected_signatures, annotated_regions = get_best_annotated_regions(combined_df,
                                                                        signatures,
                                                                        num_of_best_annotated_regions_per_signature,
                                                                        column_for_sorting,
                                                                        column_for_sorting_ascending,
                                                                        column_for_significance_level,
                                                                        significance_level)

    annotated_regions_sorted_list, signatures_array, circle_radius_array, circle_color_array = fill_arrays(combined_df,
                                                                                                            annotated_regions,
                                                                                                            signatures, # selected_signatures,
                                                                                                            column_for_circle_radius,
                                                                                                            column_for_circle_color,
                                                                                                           significance_level)


    plot_colorbar_vertical(output_path,
                           cmap,
                           colorbar_v_min,
                           colorbar_v_max,
                           annotated_regions_filename)

    circle_radius_min = np.min(circle_radius_array[np.nonzero(circle_radius_array)])
    circle_radius_max = np.max(circle_radius_array[np.nonzero(circle_radius_array)])

    plot_radiusbar_vertical(output_path,
                            circle_radius_min,
                            circle_radius_max,
                            annotated_regions_filename)

    plot_circles_on_grid(output_path,
                         annotated_regions_sorted_list,
                         signatures_array,
                         circle_radius_array,
                         circle_color_array,
                         significance_level,
                         radius_adjustment_value,
                         cmap,
                         colorbar_v_min,
                         colorbar_v_max,
                         circle_radius_max,
                         annotated_regions_filename)


# left here
# send parameters from SigProfilerTopography
# TODO modify to integrate with SigProfilerTopography
# TODO run from start to end with SigProfilerTopography
# TODO apply for protein coding genes + miRNAs
def annotated_regions_figures(genome,
                              output_path,
                              jobname,
                              numberofSimulations,
                              annotated_regions_filename,
                              log_file,
                              verbose):

    run_parallel = True
    # run_parallel = False

    genome_length = get_genome_length(genome)

    # read
    sbs_signatures_file_path = os.path.join(
        '/restricted/alexandrov-ddn/users/burcak/SigProfilerTopographyRuns/Mutographs_ESCC_552/topography/all_escc_552_samples/data/Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability.txt')
    signatures_df = pd.read_csv(os.path.join(sbs_signatures_file_path),
                                sep='\t', header=0,
                                dtype={'cutoff': np.float32,
                                       'signature': str,
                                       'number_of_mutations': np.int32,
                                       'average_probability': np.float32})

    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(
            os.path.join(outputDir, jobname, DATA,
                         Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename), sep='\t', header=0,
            dtype={'cutoff': np.float32, 'signature': str, 'number_of_mutations': np.int32,
                   'average_probability': np.float32})
    if (DBS in sigprofiler_extractor_mutation_types_contexts):
        dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(
            os.path.join(outputDir, jobname, DATA,
                         Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename), sep='\t', header=0,
            dtype={'cutoff': np.float32, 'signature': str, 'number_of_mutations': np.int32,
                   'average_probability': np.float32})
    if (ID in sigprofiler_extractor_mutation_types_contexts):
        indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(
            os.path.join(outputDir, jobname, DATA,
                         Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability_Filename), sep='\t', header=0,
            dtype={'cutoff': np.float32, 'signature': str, 'number_of_mutations': np.int32,
                   'average_probability': np.float32})

    # read
    lncRNA_file_path = os.path.join('/restricted/alexandrov-ddn/users/burcak/data/GENCODE/GRCh37/lncRNA',
                                    'gencode.v43lift37.long_noncoding_RNAs.gff3')
    annotated_region_df = read_lncRNA_GRCh37_annotation_file(lncRNA_file_path)
    annotated_region_df['length'] = annotated_region_df['end'] - annotated_region_df['start'] + 1
    annotated_region_2_length_dict = {k: v for k, v in zip(annotated_region_df['gene_name'], annotated_region_df['length'])}
    print("annotated_region_2_length_dict['CASC6']", annotated_region_2_length_dict['CASC6'])

    # real and simulated number of hits are stored in the numpy arrays in the order of ordered_annotated_regions_array
    num_of_annotated_regions, ordered_annotated_regions_array, chrom_2_tree_dict = create_intervaltree(annotated_region_df)

    print('num_of_annotated_regions:', num_of_annotated_regions)
    print('len(ordered_annotated_regions_array):', len(ordered_annotated_regions_array))

    # signatures of interest
    # signatures = ['SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS13', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS40', 'aggregatedsubstitutions']
    # signatures = ['SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS13', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS40']
    signatures = ['SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS5', 'SBS8', 'SBS10a', 'SBS10b', 'SBS12', 'SBS13', 'SBS14',
                  'SBS15', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS20', 'SBS22', 'SBS23', 'SBS25', 'SBS28',
                  'SBS30', 'SBS33', 'SBS34', 'SBS39', 'SBS40', 'SBS44', 'SBS288P']
    # signatures = ['ID1', 'ID2', 'ID3', 'ID4', 'ID6', 'ID8', 'ID9', 'ID11', 'ID14', 'ID17', 'aggregatedindels']
    # 'aggregatedsubstitutions'
    # 'aggregateddinucs'
    # 'aggregatedindels'

    # signature_type = 'indels'
    signature_type = 'subs'

    # hyperparameters
    num_of_best_annotated_regions_per_signature = 15

    # pvalue of interest for multiple testing correction
    # pvalue_column_of_interest = 'binomtest_pvalue'
    pvalue_column_of_interest = 'fisherexacttest_pvalue'

    # sort w.r.t. which column?
    # circle size w.r.t. which column?
    # circle color w.r.t. which column?
    column_for_sorting = 'real_num_of_hits' # can be foldchange or real_num_of_hits
    column_for_sorting_ascending = False
    column_for_circle_radius = 'foldchange'
    column_for_circle_color = 'qvalue'
    column_for_significance_level = 'qvalue'
    significance_level = 0.05
    radius_adjustment_value = 2 # must be modified wrt range of column_for_sorting
    colorbar_v_min = 0
    colorbar_v_max = 20

    cmap = cm.get_cmap('Reds')
    cmap = truncate_colormap(cmap, 0.2, 1)
    output_path = os.path.join('/oasis/tscc/scratch/burcak/temp/Mutographs_ESCC_552_lncRNA')

    plot_highly_mutated_annotated_regions_figure(annotated_regions_filename,
        genome_length,
        signatures,
        signature_type,
        signatures_df,
        ordered_annotated_regions_array,
        annotated_region_2_length_dict,
        run_parallel,
        num_of_best_annotated_regions_per_signature,
        pvalue_column_of_interest,
        column_for_sorting,
        column_for_sorting_ascending,
        column_for_circle_radius,
        column_for_circle_color,
        column_for_significance_level,
        significance_level,
        radius_adjustment_value,
        cmap,
        colorbar_v_min,
        colorbar_v_max,
        output_path)