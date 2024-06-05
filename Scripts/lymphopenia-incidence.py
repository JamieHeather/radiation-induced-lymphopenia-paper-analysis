#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
lymphopenia-incidence.py

Look at the incidence of lymphopenia in radiation-naive patients in the larger NGH cohort

"""


import matplotlib.pyplot as plt
import seaborn as sns
import functions as fxn


__version__ = '0.3.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

# Sort plotting parameters
# sns.set(style="darkgrid")
plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial'})

if __name__ == "__main__":

    fxn.check_scripts_dir()
    plot_dir = fxn.plot_dir('lymphopenia-incidence')

    csv_file = fxn.raw_data_dir + 'chemoRT-lymphopenia.csv'
    lp = []

    with open(csv_file, 'r') as in_file:
        line_count = 0
        for line in in_file:
            bits = line.rstrip().split(',')
            if line_count == 0:
                bits[0] = 'ID'
                headers = bits
            else:
                # Convert to the appropriate types
                grade_index = headers.index([x for x in headers if 'grade' in x][0])
                alc_index = headers.index([x for x in headers if 'ALC' in x][0])
                bits[grade_index] = int(bits[grade_index])
                bits[alc_index] = float(bits[alc_index])
                lp.append(bits)

            line_count += 1

    lp = fxn.list_to_df(lp, headers, True)

    num_pr = sum(lp['xRT type'] == 'Protons')
    num_ph = sum(lp['xRT type'] == 'Photons')

    p_alc = fxn.mwu(lp.loc[lp['xRT type'] == 'Protons']['ALC at nadir (K/uL)'],
                    lp.loc[lp['xRT type'] == 'Photons']['ALC at nadir (K/uL)'])

    p_grade = fxn.mwu(lp.loc[lp['xRT type'] == 'Protons']['Lymphopenia grade'],
                      lp.loc[lp['xRT type'] == 'Photons']['Lymphopenia grade'])

    fig = plt.figure(figsize=(2.6, 5))
    ax = sns.boxplot(data=lp, x='xRT type', y='Lymphopenia grade', hue='xRT type', fliersize=0)
    ax = sns.stripplot(data=lp, x='xRT type', y='Lymphopenia grade', hue='xRT type', alpha=.5, jitter=.4,
                       linewidth=1, edgecolor='black')
    fxn.plot_signif_lines(0, 1, 4.15, p_grade, 'black')
    plt.xlabel('')
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'box-grade.svg', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(2.6, 5))
    ax = sns.violinplot(data=lp, x='xRT type', y='ALC at nadir (K/uL)', cut=0, hue='xRT type')
    fxn.plot_signif_lines(0, 1, 2.05, p_alc, 'black')
    plt.xlabel('')
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'violin-ALC.svg', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(2.6, 5))
    ax = sns.violinplot(data=lp, x='xRT type', y='ALC at nadir (K/uL)', cut=0, hue='xRT type')
    ax = sns.stripplot(data=lp, x='xRT type', y='ALC at nadir (K/uL)', hue='xRT type', alpha=.2, jitter=.35,
                       linewidth=1, edgecolor='black')
    fxn.plot_signif_lines(0, 1, 2.05, p_alc, 'black')
    plt.xlabel('')
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'violin-ALC-plus-points.svg', dpi=300, bbox_inches='tight')
    plt.close()

    # And a horizontal version
    fig = plt.figure(figsize=(5, 2.6))
    ax = sns.violinplot(data=lp, y='xRT type', x='ALC at nadir (K/uL)', cut=0, orient='h', hue='xRT type')
    plt.plot([2.02, 2.02], [0, 1], color='black', alpha=.7)
    plt.plot([2, 2.02], [0, 0], color='black', alpha=.7)
    plt.plot([2, 2.02], [1, 1], color='black', alpha=.7)
    plt.annotate(fxn.asterisky(p_alc), xy=(1.9, 0.5), xytext=(1.9, 0.5),
                 fontsize=15, color='black', va='center', ha='center', alpha=.7)
    plt.ylabel('')
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'violin-ALC-horiz.svg', dpi=300, bbox_inches='tight')
    plt.close()

    # Get % of each xRT type with each diagnosis
    diagnoses = list(set(lp['Diagnosis']))
    diagnoses.sort()
    prop_diagnoses = []
    for gp in diagnoses:
        for xrt in ['Photons', 'Protons']:
            # For each treatment type, record the percentage of each radiotherapy group that received it
            prop_diagnoses.append([xrt, gp, (len(lp[(lp['xRT type'] == xrt) & (lp['Diagnosis'] == gp)]) /
                                             len(lp[lp['xRT type'] == xrt])) * 100])

    prop_diagnoses = fxn.list_to_df(prop_diagnoses, ['xRT type', 'Diagnosis', 'Percentage'], False)

    fig = plt.figure(figsize=(6, 5))
    sns.barplot(data=prop_diagnoses, x='Diagnosis', y='Percentage', hue='xRT type')
    plt.xlabel('')
    plt.xticks(rotation=20)
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'bar-diagnosis-percents.svg', dpi=300, bbox_inches='tight')
    plt.close()

    # Plot proportion of donors with each grade value as a barplot
    grades = list(set(lp['Lymphopenia grade']))
    grades.sort()
    prop_grades = []
    for g in grades:
        for xrt in ['Photons', 'Protons']:
            prop_grades.append([xrt, g, (len(lp[(lp['xRT type'] == xrt) & (lp['Lymphopenia grade'] == g)]) /
                                             len(lp[lp['xRT type'] == xrt])) * 100])

    prop_grades = fxn.list_to_df(prop_grades, ['xRT type', 'Diagnosis', 'Percentage'], False)

    fig = plt.figure(figsize=(6, 5))
    sns.barplot(data=prop_grades, x='Diagnosis', y='Percentage', hue='xRT type', legend=False)
    plt.xlabel('')
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'bar-grades-percents.svg', dpi=300, bbox_inches='tight')
    plt.close()
