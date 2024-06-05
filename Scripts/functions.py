# -*- coding: utf-8 -*-


import os
import datetime
import sys
import gzip
import pickle
import numpy as np
import pandas as pd
import subprocess
import collections as coll
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
from scipy import stats


__version__ = '0.6.2'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def nest_counter():
    """
    Create nested counters
    """
    return coll.Counter()


def nest_ddict():
    """
    Returns a defaultdict to allow creation of nested dictionaries
    """
    return coll.defaultdict(list)


def opener(file_path, open_mode):
    """
    :param file_path: path to file to be opened
    :param open_mode: mode by which to open the file (e.g. w/r/a)
    :return: the appropriate file opening command (open or gzip.open)
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, open_mode + 't')
    else:
        return open(file_path, open_mode)


def list_to_df(list_to_df, headers, rename):
    """
    Convert a list to a (long) dataframe. Note that first entry becomes the index if chosen
    :param list_to_df: List of list entries (with each position in each list corresponding to a column)
    :param headers: List of column headers. First column should be unique, becoming the rownames, if rename = True
    :param rename: Option to rename row IDs by first colum
    :return: sorted pandas dataframe
    """
    df = pd.DataFrame(list_to_df)
    df = df.rename(index=str, columns=dict(list(zip(list(range(len(headers))), headers))))
    df = df.sort_values(by=[headers[0]])
    if rename == True:
        df = df.set_index(headers[0], drop=True)
    return df


def open_pickle(file_name):
    """
    :param file_name: the file name of the pickle archive to be opened
    :return: the nested peptide dictionary it contained
    """
    pkl_file = opener(file_name, 'rb')
    pkl_dat = pickle.load(pkl_file)
    pkl_file.close()
    return pkl_dat


def save_pickle(file_name, dict):
    """
    :param file_name: the file name to save the pkl as
    :param dict: the dictionary to be saved
    :return:
    """
    if not file_name.endswith('.pkl'):
        file_name = file_name + '.pkl'
    output = opener(file_name, 'wb')
    pickle.dump(dict, output)
    output.close()


def check_scripts_dir():
    """
    Check we're in the right directory (Scripts)
    :return:
    """
    if not os.getcwd().endswith('/Scripts'):
        if 'Scripts' in os.listdir(os.getcwd()):
            os.chdir('Scripts')
        else:
            print("Check your current working directory - this is designed to be run from root or Scripts folders")
            sys.exit()


def get_timepoint(donor, date, meta):
    """
    Given a particular sample date, return which of the three time-points it corresponds to for that patient
    :param donor: patient ID
    :param date: date of particular sample
    :param meta: dataframe of patient metadata
    :return: 1/2/3 (/error), depending on which of the three time points this sample corresponds to, using the metadata
    """
    if meta.loc[donor, 'TP1-date'] == date:
        return 1
    elif meta.loc[donor, 'TP2-date'] == date:
        return 2
    elif meta.loc[donor, 'TP3-date'] == date:
        return 3
    else:
        print("Error: unable to match donor to sample date:\t", donor + '\t' + date)
        print("Maybe one of these:\t\t\t\t\t" + meta.loc[donor, 'TP1-date'] + \
              ' ' + meta.loc[donor, 'TP3-date'] + ' ' + meta.loc[donor, 'TP3-date'])


def get_limits(df):
    mins = {}
    maxs = {}
    for x in df:
        try:
            if min(df[x]) > 0:
                min_len = len(str(min(df[x])).replace('-', '').split('.')[0])
                mins[x] = round(min(df[x]), -min_len + 1)
            else:
                mins[x] = round(min(df[x]), 1)
            max_len = len(str(max(df[x])).replace('-', '').split('.')[0])
            # maxs[x] = round(max(df[x]), max_len - 1)
            maxs[x] = round_up(max(df[x]))
        except Exception:
            pass
    return mins, maxs


def convert_tcr_data():
    # TODO docstring

    # First check whether the files are there or not first
    conv_tsvs = [x for x in os.listdir(conv_data_dir) if '.tsv' in x]
    if len(conv_tsvs) == 60:
        return

    # If not, continue to generate the required AIRR-seq format TSVs
    print("No AIRR-seq format TSV files detected: converting from raw immunoSEQ data")
    all_in_files = [x for x in os.listdir(raw_data_dir) if '.tsv' in x]
    all_in_files.sort()

    make_check_dir(conv_data_dir)

    for f in all_in_files:
        base = f.split('.')[0]
        cmd = 'python immunoseq2airr.py -i ' + raw_data_dir + f + ' -o ' + conv_data_dir + base + ' -z -or'
        subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

        # If the file isn't already gzipped, make it so
        if f.endswith('.tsv'):
            with opener(raw_data_dir + f, 'r') as in_file, gzip.open(raw_data_dir + f + '.gz', 'wb') as zipped:
                zipped.writelines(in_file)

            os.unlink(raw_data_dir + f)


def read_cdr3s_in(path_to_file):
    """
    :param path_to_file:
    :return: counter containing all the CDR3s in the input file and their abundances, and another of nt sequences
    """

    cdr3_dat = coll.Counter()
    nt_dat = coll.Counter()
    line_count = 0
    with opener(path_to_file, 'r') as in_file:

        for line in in_file:
            bits = line.rstrip().split('\t')
            if line_count == 0:
                headers = bits
                params = {}
                for x in range(len(headers)):
                    params[headers[x]] = x

            else:
                v = remove_alleles(bits[params['v_call']])
                j = remove_alleles(bits[params['j_call']])
                cdr3 = bits[params['junction_aa']]
                nt = bits[params['sequence']]
                count = int(bits[params['duplicate_count']])

                # Store all nucleotide seqs + counts, regardless of functionality
                nt_dat[nt] += count

                # Productive rearrangements are stored as v|j|cdr3 strings
                if cdr3:
                    cdr3_dat[v + '|' + j + '|' + cdr3] += count

            line_count += 1

    return cdr3_dat, nt_dat


def remove_alleles(gene_calls):
    """
    :param gene_calls: A V or J gene call, which could include a list of ambiguous gene/allele calls
    :return: Same gene calls, with allele-level information stripped away, sorted alphabetically
    """
    gene_list = list(set([x.split('*')[0] for x in gene_calls.split(',')]))
    gene_list.sort()
    return ','.join(gene_list)


def write_cdr3s_out(dat, outfile_path):
    """
    Output combined CDR3 to a new TSV file
    :param dat: counter containing CDR3s/counts
    :param outfile_path: path to file to write
    :return:
    """

    with opener(outfile_path, 'w') as outfile:

        outfile.write('Blank\tCDR3\tCount\n')

        for cdr3 in dat:
            outfile.write('\t' + cdr3 + '\t' + str(dat[cdr3]) + '\n')

        return


def gini(abundance_list):
    """
    :param abundance_list: list of abundances
    :return: Gini index
    """
    sorted_list = sorted([float(x) for x in abundance_list])
    height, area = 0, 0
    for value in sorted_list:
        height += value
        area += height - value / 2
    gini_area = height * len(abundance_list) / 2
    return (gini_area - area) / gini_area


def jaccard(list1, list2):
    """
    :param list1: list of sequences
    :param list2: another list of sequences
    :return: Jaccard index of list1 and list2
    """
    intersection = list(set(list1) & set(list2))
    union = list(set(list1) | set(list2))
    return len(intersection) / len(union)


def stitch_vj_cdr3(v, j, cdr3, delimiter):
    """
    Turn TCR info into a string (which behaves better in certain data structures than a tuple or list)
    :param v: TRBV(s) used
    :param j: TRBJ(s) used
    :param cdr3: TCR junction (i.e. CDR3 inclusive of C/F)
    :param delimiter: Character used to separate V/J/CDR3
    :return: TCR components stitched together as a string given a particular character
    """
    return v + delimiter + j + delimiter + cdr3


def mwu(list1, list2):
    """
    :param list1: One list of numbers
    :param list2: A second list of numbers
    :return: P value for a Mann-Whitney U test (unpaired, non-parametric)
    """
    return stats.mannwhitneyu(list(list1), list(list2))[1]


def wc(list1, list2):
    """
    :param list1: One list of numbers
    :param list2: A second list of numbers
    :return: P value for a Wilcoxon test (paired, non-parametric)
    """
    return stats.wilcoxon(list(list1), list(list2))[1]


def nest():
    """
    Returns a defaultdict to allow creation of nested dictionaries
    """
    return coll.defaultdict()


def nest_list():
    """
    Returns a defaultdict to allow creation of nested dictionaries
    """
    return coll.defaultdict(list)


def asterisky(pval):
    """
    :param pval: A p-value (float) from a statistical test
    :return: a string (- or *s) based on arbitrary significance thresholds
    """
    if pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    else:
        return '-'


def plot_signif_lines(x1, x2, y, p_val, col):
    """
    Adds significance lines/text to an existing plot (much more manually than plot_3tp_signif)
    :param x1: X position of first group
    :param x2: X position of second group
    :param y: Y axis position
    :param p_val: P value
    :param col: Colour
    :return: Nothing - plot alterations are only output
    """
    if p_val < 0.05:
        ytix = plt.gca().get_ylim()
        lower = ytix[0]
        higher = ytix[-1]
        yrange = higher - lower
        downlength = yrange * 0.005
        offset = 0.1
        text_x = ((x2 - x1) / 2) + x1
        plt.plot([x1, x2], [y, y], color=col, alpha=.7)
        plt.plot([x1, x1], [y, y - downlength], color=col, alpha=.7)
        plt.plot([x2, x2], [y, y - downlength], color=col, alpha=.7)
        plt.annotate(asterisky(p_val), xy=(text_x, y), xytext=(text_x, y),
                     fontsize=15, color=col, va='center', ha='center', alpha=.7)


def plot_3tp_signif(df, population, save_path, y_label, dat_type, sci_y=False):
    """
    Plots a value across the 3 time points (or 3 time point transitions), with colour-coded significance bars
    :param df: Dataframe containing the data to be plotted
    :param population: Population (i.e. column) to be plotted
    :param save_path: Path to final file to save
    :param y_label: Y axis label
    :param dat_type: type of data input, i.e. diversity (for shannon/gini) or jaccard
    :param sci_y: boolean, whether to force scientific notion on Y axis
    :return: Nothing - saved plot is output
    """
    type_x = {'diversity': 'Timepoint', 'jaccard': 'Transition'}
    cols_x = {'diversity': [(1, 2), (2, 3), (1, 3)], 'jaccard': [('1-2', '2-3'), ('2-3', '1-3'), ('1-2', '1-3')]}
    orders = {'diversity': [1, 2, 3], 'jaccard': ['1-2', '2-3', '1-3']}

    p = sns.catplot(x=type_x[dat_type], y=population, hue='Treatment', cut=0, order=orders[dat_type],
                kind='violin', inner='stick', split=True, data=df, hue_order=['Photons', 'Protons'],
                    legend=False, density_norm='count')
    final_ylabel = ' | '.join(population.split('|')[0:2]) + ' ' + y_label
    plt.ylabel(final_ylabel)
    ytix = plt.gca().get_ylim()
    lower = ytix[0]
    higher = ytix[-1]
    yrange = higher - lower
    downlength = yrange * 0.005

    # per timepoint (i.e. Ph1 vs Pr1, Ph2 vs Pr2, Ph3 vs Pr3)
    TP_height = higher - (yrange * .05)
    TP_text_height = higher - (yrange * .04)
    TP_offset = .1
    for i in range(3):

        # Calculate the pval
        tmp = df.loc[df[type_x[dat_type]] == orders[dat_type][i]]
        ph = tmp.loc[tmp.Treatment == 'Photons'][population]
        pr = tmp.loc[tmp.Treatment == 'Protons'][population]

        if sum(ph) != 0 and sum(pr) != 0:
            pval = mwu(ph, pr)

            # Then if that's significant
            if asterisky(pval) != '-':
                # Plot the bars indicating test fields
                plt.plot([i - TP_offset, i + TP_offset], [TP_height, TP_height], color='black', alpha=.7)
                plt.plot([i - TP_offset, i - TP_offset], [TP_height, TP_height - downlength], color='black', alpha=.7)
                plt.plot([i + TP_offset, i + TP_offset], [TP_height, TP_height - downlength], color='black', alpha=.7)
                plt.annotate(asterisky(pval), xy=(i, TP_text_height), xytext=(i, TP_text_height),
                             fontsize=13, color='black', va='center', ha='center', alpha=.7)

    # per treatment (e.g. Ph1 vs Ph2, Ph2 vs Ph3 ...)
    p_offset = {'Photons': -TP_offset, 'Protons': TP_offset}
    p_cols = {'Photons': blue, 'Protons': orange}
    p_heights = {'Photons': higher, 'Protons': higher + (yrange * .05)}

    for tr in ['Photons', 'Protons']:
        for t in cols_x[dat_type]:
            # Calculate the pval
            tmp1 = df.loc[(df[type_x[dat_type]] == t[0]) & (df.Treatment == tr)][population]
            tmp2 = df.loc[(df[type_x[dat_type]] == t[1]) & (df.Treatment == tr)][population]

            if t == (1, 3) or t == ('1-2', '1-3'):
                y_offset = yrange * .025
            else:
                y_offset = 0

            if sum(tmp1) != 0 and sum(tmp2) != 0:
                pval = wc(tmp1, tmp2)
                if asterisky(pval) != '-':
                    # Plot the lines
                    x1 = orders[dat_type].index(t[0])
                    x2 = orders[dat_type].index(t[1])

                    plt.plot([x1 + p_offset[tr], x2 + p_offset[tr]],
                             [p_heights[tr] + y_offset, p_heights[tr] + y_offset],
                             color=p_cols[tr], alpha=.7)
                    plt.plot([x1 + p_offset[tr], x1 + p_offset[tr]],
                             [p_heights[tr] + y_offset, p_heights[tr] + y_offset - downlength],
                             color=p_cols[tr], alpha=.7)
                    plt.plot([x2 + p_offset[tr], x2 + p_offset[tr]],
                             [p_heights[tr] + y_offset, p_heights[tr] + y_offset - downlength],
                             color=p_cols[tr], alpha=.7)
                    # Then add the text
                    index_tuple = [orders[dat_type].index(x) for x in t]
                    plt.annotate(asterisky(pval),
                                 xy=(x_pos_3tp_annot(index_tuple) + p_offset[tr], p_heights[tr] + y_offset),
                                 xytext=(x_pos_3tp_annot(index_tuple) + p_offset[tr], p_heights[tr] + y_offset),
                                 fontsize=13, color=p_cols[tr], va='center', ha='center', alpha=.7)

    plt.xlim(-0.5, 2.5)
    plt.subplots_adjust(bottom=0.13, left=0.20, right=0.98, top=0.9)

    if sci_y:
        for ax in p.axes.flatten():
            ax.ticklabel_format(style='scientific', scilimits=(0, 0), axis='y', useMathText=True)

    plt.savefig(save_path, dpi=300)#, bbox_inches='tight')
    plt.close()


def x_pos_3tp_annot(tp_tuple):
    """
    :param tp_tuple: Tuple of first and second timepoint numbers (i.e. 1 or 2 and 2 or 3 respectively)
    :return: Correct x axis position for the corresponding significance annotation (midway between the x axis points)
    """
    a = tp_tuple[0]
    b = tp_tuple[1]
    return a + ((b - a)/2)


def today():
    """
    :return: Today's day, in ISO format
    """
    return datetime.datetime.today().date().isoformat()


def date(text_date):
    """
    :param text_date: a string containing an ISO-formatted date
    :return: a datetime object (with which additions/subtractions can be easily performed)
    """
    return datetime.datetime.strptime(text_date, date_format).date()


def plot_dir(plot_dir_suffix):
    """
    :return: The path to the plots directory subfolder for results plotted on this day (creating it if needed)
    """
    plot_dir = base_plot_dir + today() + '-' + plot_dir_suffix + '/'
    make_check_dir(plot_dir)
    return plot_dir


def make_check_dir(dir_path):
    """
    :param dir_path: A directory to check whether it exists; if not, make it
    """
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    if dir_path[-1] != '/':
        dir_path += '/'

    return dir_path


def get_metadata():
    """
    :return: Pandas dataframe containing the patient cohort's metadata
    """
    print("Getting metadata...")
    meta = pd.read_csv(raw_data_dir + 'metadata.csv')
    meta = meta.set_index([x for x in meta][0], drop=True)
    return meta


def flatten_list(in_list):
    """
    :param in_list: A list of lists (only one level deep)
    :return: A flat list, having unpacked each of the inner list
    """
    return [thing for inner_list in in_list for thing in inner_list]


def plot_jointplot_stats_positive(df, save_path, x_axis, y_axis, x_lab, y_lab):
    """
    Saves jointplots featuring R2/p value stats, for positive values (as x/y axes cut at zero)
    :param df: Pandas dataframe containing columns to be correlated
    :param save_path: Path to file to save plot to
    :param x_axis: Column to plot on X axis
    :param y_axis: Column to plot on Y axis
    :param x_lab: X axis label
    :param y_lab: Y axis label
    :return: Nothing - saves plot to desired position
    """

    sns.jointplot(data=df, x=x_axis, y=y_axis, kind='reg', stat_func=stats.pearsonr)
    plt.xlim(0, None)
    plt.ylim(0, None)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def get_regression_str(x, y):
    """
    :param x: list of X values
    :param y: list of y values
    :return: a string detailing the R squared and P values of the linear regression of X and Y
    """
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r_2 = r_value * r_value
    return "R2: " + "{0:.2f}".format(r_2) + '\nP: ' + format_p_vals(p_value)


def format_p_vals(number):
    """
    :param number: a float between 0 and 1 (a p value)
    :return: string of formatted p value (to 3 sig figs if > 0.001, otherwise in scientific notation)
    """
    if number > 0.001:
        return "{0:.3f}".format(number)
    else:
        return "{:.3E}".format(number)


def round_up(num):
    """
    :param num: numeric (e.g. float) to be rounded up
    :return: int, input number rounded up to the nearest integer
    """
    return int(num) + (num % 1 > 0)


def save_csv(dataframe, path_to_save, df_name):
    """
    :param dataframe: Pandas dataframe to be saved
    :param path_to_save: path to csv being produced
    :param df_name: identifiable name to save the dataframe under
    :return:
    """
    destination = path_to_save + today() + '-' + df_name + '.csv'
    dataframe.to_csv(destination, index=False, encoding='utf-8')
    return destination


def open_previous_data(directory, data_name, data_type):
    """
    :param directory: directory to look for the file in
    :param data_name: the named portion of a file produced by one of the upstream scripts
    :param data_type: 'csv' or 'pkl'
    :return: the most recent CSV file fitting that name, opened as panas dataframe
    """
    if data_type.lower() not in ['csv', 'pkl']:
        raise IOError('Inappropriate data type not in acceptable list (csv/pkl): \"' + data_type + '\"')

    all_hits = [x for x in os.listdir(directory) if x.endswith('.' + data_type) and data_name.upper() in x.upper()]
    if all_hits:
        all_hits.sort()
        most_recent = all_hits[-1]
        if data_type == 'csv':
            return pd.read_csv(directory + most_recent)
        elif data_type == 'pkl':
            return open_pickle(directory + most_recent)
    else:
        raise IOError('Cannot find a file matching \"' + data_name + '\" in ' + directory)


db_dir = '../CDR3-Databases/'
raw_data_dir = '../Raw-Data/'
conv_data_dir = '../Converted-Data/'
base_plot_dir = '../Plots/'
vdjdb_dir = raw_data_dir + 'VDJdb/'
date_format = '%Y-%m-%d'
stitch_character = '|'
blue = (0.2980392156862745, 0.4470588235294118, 0.6901960784313725)
orange = (0.8666666666666667, 0.5176470588235295, 0.3215686274509804)
alt_cols = sns.xkcd_palette(["greyish", "dusty purple", "amber"])
