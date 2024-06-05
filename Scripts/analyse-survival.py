#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
analyse-survival.py

Performs basic survival analyses on the xRT patient cohort.
Additionally also plots some final graphs that are inconvenient to insert elsewhere.

"""


import numpy as np
import collections as coll
import matplotlib.pyplot as plt
import seaborn as sns
import functions as fxn
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter


__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

# Sort plotting parameters
plt.rcParams.update({'font.size': 20, 'font.sans-serif': 'Arial'})

if __name__ == "__main__":

    # Sort out the necessary directories, checking we're in the right one
    fxn.check_scripts_dir()
    today = fxn.today()
    base_plot_dir = fxn.make_check_dir(fxn.base_plot_dir)
    plot_dir = fxn.plot_dir('analyse-survival')

    # Get previous data
    trans_order = ['1-2', '2-3']  #, '1-3']
    tp_data = fxn.open_previous_data(fxn.raw_data_dir, 'flow-data', 'csv')
    delta_data = fxn.open_previous_data(fxn.conv_data_dir, 'deltas', 'pkl')
    meta_data = fxn.get_metadata()

    delta_data = delta_data.loc[delta_data['Transition'].isin(trans_order)]

    convert = {'OS': 'Survival', 'OS-time-mo': 'Time from diagnosis to last follow up (mo)',
               'OS-time-days': 'Time from diagnosis to last follow up (days)'}

    for f in ['OS', 'OS-time-mo', 'OS-time-days']:
        entries = fxn.flatten_list([[x]*2 for x in meta_data[f]])  # Change 2->3 if including 1-3 transitions
        if f == 'OS':
            entries = [str(x).replace('0', 'Alive').replace('1', 'Deceased') for x in entries]
        delta_data[convert[f]] = entries

    g = sns.lmplot(data=delta_data, y='ΔGini', x='Time from diagnosis to last follow up (days)', hue='Survival',
                   col='Transition', palette=fxn.alt_cols)
    plt.ylim(-0.2, 0.6)
    plt.xlim(0, 5000)

    for ax, trans in zip(g.axes.flat, trans_order):
        subset = delta_data[delta_data['Transition'] == trans]
        live = subset[subset['Survival'] == 'Alive']
        died = subset[subset['Survival'] == 'Deceased']
        txt_y = 1.02
        ax.text(0, txt_y, fxn.get_regression_str(live['ΔGini'], live['Time from diagnosis to last follow up (days)']),
                color=fxn.alt_cols[0], transform=ax.transAxes, fontsize=11, fontweight='bold')
        ax.text(.99, txt_y, fxn.get_regression_str(died['ΔGini'], died['Time from diagnosis to last follow up (days)']),
                color=fxn.alt_cols[1], transform=ax.transAxes, fontsize=11,
                horizontalalignment='right', fontweight='bold')

    plt.savefig(plot_dir + 'OStime-columns.svg', dpi=300, bbox_inches='tight')
    plt.close()

    # Plot the wider time-frame clinical values
    to_scatter = coll.defaultdict(list)
    day_lists = []

    for donor in meta_data.index:

        diagnosis = meta_data.loc[donor]['diagnosis']
        treatment = meta_data.loc[donor]['xRT-type']
        alc_xy_list = []
        cea_xy_list = []
        grade_xy_list = []

        for i in range(1, 9):
            tp = 'TP' + str(i)
            alc_x = meta_data.loc[donor]['rel-' + tp + '-date']
            alc_y = meta_data.loc[donor][tp + '-ALC-K/uL']
            grade_y = meta_data.loc[donor][tp + '-grade']
            if alc_x or alc_y:
                alc_xy_list.append((alc_x, alc_y))
                to_scatter['alc'].append([donor, diagnosis, treatment, alc_x, alc_y])
            if grade_y:
                grade_xy_list.append((alc_x, grade_y))
                to_scatter['grade'].append([donor, diagnosis, treatment, alc_x, grade_y])
            if i > 3:  # Only have CEA data from TP4 onwards
                cea_x = meta_data.loc[donor]['rel-' + tp + '-CEA-date']
                cea_y = meta_data.loc[donor][tp + '-CEA']
                if cea_x or cea_y:
                    cea_xy_list.append((cea_x, cea_y))
                    to_scatter['cea'].append([donor, diagnosis, treatment, cea_x, cea_y])
        day_lists.append([donor, diagnosis, treatment, alc_xy_list, cea_xy_list, grade_xy_list])

    tmp_headers = ['Donor', 'Diagnosis', 'Treatment', 'Day']
    value_key = {'alc': 'ALC (K/uL)', 'cea': 'CEA', 'grade': 'Lymphopenia grade'}
    day_lists = fxn.list_to_df(day_lists, ['Donor', 'Diagnosis', 'Treatment', 'alc', 'cea', 'grade'], False)

    for value in ['alc', 'cea', 'grade']:
        to_scatter[value] = fxn.list_to_df(to_scatter[value], tmp_headers + [value_key[value]], False)

        fig = plt.figure(figsize=(11, 6))
        ax = fig.add_subplot(111)

        sns.scatterplot(x='Day', y=value_key[value], hue='Treatment', style='Diagnosis',
                        data=to_scatter[value], s=100, legend='brief')
        for dl in day_lists.index:
            if day_lists.loc[dl]['Treatment'] == 'Photons':
                col = fxn.blue
            elif day_lists.loc[dl]['Treatment'] == 'Protons':
                col = fxn.orange

            # To ensure proper joining the dots need to sort on the x axis (date), omitting missing entries
            to_plot = day_lists.loc[dl][value]
            to_plot = sorted([x for x in to_plot if not np.isnan(x[0]) and not np.isnan(x[1])])

            for i in range(len(to_plot)):
                if len(to_plot) - i != 1:
                    plt.plot((to_plot[i][0], to_plot[i + 1][0]), (to_plot[i][1], to_plot[i + 1][1]),
                             color=col, linestyle='--', alpha=.8, linewidth=1)

        plt.xlabel('Days (relative to nadir)')
        if value == 'cea':
            plt.yscale('log')

        sns.despine(right=True, top=True)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(plot_dir + value + '-over-time.svg', dpi=300, bbox_inches='tight')
        plt.close()

    # ALC , gini, shannon, CEA, T cell %, T reg %s, bars for treatment etc, last follow up dates ...
    threshold = .1
    threshold_top = 0.1

    test = delta_data.loc[delta_data['Transition'] == '2-3']
    test['OS'] = [''] * len(test)
    test.loc[test.loc[test['Survival'] == 'Alive'].index, 'OS'] = 0
    test.loc[test.loc[test['Survival'] == 'Deceased'].index, 'OS'] = 1

    text_thresh = '{0:.2f}'.format(threshold)
    print('\t', text_thresh, '--------------------------------------------------------------===========********')

    test['Clonal'] = [''] * len(test)
    test.loc[test.loc[test['ΔGini'] > threshold_top].index, 'Clonal'] = 1
    test.loc[test.loc[test['ΔGini'] <= threshold].index, 'Clonal'] = 0

    test = test.loc[test['Clonal'].isin([0, 1])]

    chk_key = {0: 'No', 1: 'Yes'}
    test['Oligoclonality increase?'] = [chk_key[x] for x in test['Clonal']]

    num_clonal = sum(test['Clonal'])
    num_notclonal = len(test) - num_clonal

    T = test["Time from diagnosis to last follow up (mo)"]
    E = test["OS"]

    kmf = KaplanMeierFitter()
    kmf.fit(durations=T, event_observed=E)

    ax = plt.subplot(111)
    c = (test["Clonal"] == 0)
    # print(E, c)
    kmf.fit(durations=T[c], event_observed=E[c], label="ΔGini <= " + text_thresh)
    kmf.plot_survival_function(ax=ax, show_censors=True, ci_show=True, color=fxn.alt_cols[0])
    kmf.fit(T[~c], event_observed=E[~c], label="ΔGini > " + text_thresh)
    kmf.plot_survival_function(ax=ax, show_censors=True, ci_show=True, color=fxn.alt_cols[1])  #, at_risk_counts=True)

    ax.text(0.75, 1.02, 'n = ' + str(num_clonal), color=fxn.alt_cols[0],
            transform=ax.transAxes, fontsize=13, fontweight='bold')
    ax.text(0.9, 1.02, 'n = ' + str(num_notclonal), color=fxn.alt_cols[1],
            transform=ax.transAxes, fontsize=13, fontweight='bold')

    test2 = test[['OS', 'Clonal']]
    test2['Time'] = list(test["Time from diagnosis to last follow up (mo)"])

    cph = CoxPHFitter()
    cph.fit(test2, duration_col="Time", event_col='OS')
    cph.print_summary()
    pval = '{0:.2f}'.format(cph._compute_p_values()[0])
    ax.text(0.3, 1.02, 'p = ' + pval, color='r', transform=ax.transAxes, fontsize=13, fontweight='bold')

    plt.xlabel("Time (mo)")
    plt.ylabel("Survival probability")
    plt.savefig(plot_dir + 'Kaplan-Meier-threshold-' + text_thresh + '.svg', dpi=300, bbox_inches='tight')
    plt.close()

    # Plot just the general distribution of Gini deltas
    fig = plt.figure(figsize=(3, 5))
    ax = sns.violinplot(data=test, y='ΔGini', cut=0, color='lightgray', alpha=.5)
    ax.collections[::2][0].set_alpha(0.5)  #
    sns.stripplot(data=test, hue='Oligoclonality increase?', y='ΔGini', size=10, jitter=.35,
                  palette=fxn.alt_cols, hue_order=['No', 'Yes'])
    plt.ylabel('Post nadir ΔGini')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(plot_dir + 'raw-delta-ginis.svg', dpi=300, bbox_inches='tight')
    plt.close()
