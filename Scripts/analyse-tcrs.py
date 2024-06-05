#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
analyse-tcrs.py

Perform the bulk of the analysis, incorporating the TCRseq and related metadata analysis

Note that this specifically requires version 0.9.0 of seaborn

"""

import warnings
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import os
import collections as coll
import functions as fxn
warnings.filterwarnings("ignore")


__version__ = '0.8.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

# Sort plotting parameters
# sns.set(style="darkgrid")
plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial'})

if __name__ == "__main__":

    # Sort out the necessary directories, checking we're in the right one
    fxn.check_scripts_dir()
    today = fxn.today()
    base_plot_dir = fxn.make_check_dir(fxn.base_plot_dir)
    plot_dir = fxn.plot_dir('analyse-tcrs')

    # Get metadata
    meta = fxn.get_metadata()

    results_headers = ['Donor', 'Timepoint', 'Diagnosis', 'Treatment',
                       'Unique-TCRs', 'Total-TCRs', 'Gini', 'Shannon',
                       'ALC', 'Day']

    # Find all relevant files and read in to nested dict/pandas df - making use of previous versions if saved
    print("Reading in TCR data ...")
    # Check the converted TSVs are present; if not, generate
    fxn.convert_tcr_data()

    raw_dat_pkls = [x for x in os.listdir(fxn.conv_data_dir) if x.endswith('-raw-tcr-dict.pkl')]
    raw_dat_pkls.sort()
    results_pkls = [x for x in os.listdir(fxn.conv_data_dir) if x.endswith('-tcr-results.pkl')]
    results_pkls.sort()

    if raw_dat_pkls and results_pkls:

        # Takes the most recent entries, requiring pkls for both
        print("\tFound existing TCR raw data/results pickles, reading in:")
        print('\t\t' + fxn.conv_data_dir + raw_dat_pkls[-1])
        print('\t\t' + fxn.conv_data_dir + results_pkls[-1])
        raw_dat = fxn.open_pickle(fxn.conv_data_dir + raw_dat_pkls[-1])
        results = pd.read_pickle(fxn.conv_data_dir + results_pkls[-1])

    else:

        # Otherwise read data in fresh (and make a new pickle)
        print("\tReading in raw TCR data from AIRR Community standard format files")
        all_data_files = [x for x in os.listdir(fxn.conv_data_dir) if x.endswith('.tsv.gz')]
        all_data_files.sort()

        raw_dat = {}
        results = []
        day_vals = pd.DataFrame(index=list(meta.index), columns=['d1', 'g1', 's1', 'd2', 'g2', 's2', 'd3', 'g3', 's3'])

        for fl in all_data_files:
            sample = fl.split('.')[0]
            donor, treatment, tp, day, pad = sample.split('_')
            timepoint = int(tp[-1])
            day = int(day.replace('d', ''))

            # If donor entry isn't present in top level (non-default-)dict, create it
            if donor not in raw_dat:
                raw_dat[donor] = coll.defaultdict(fxn.nest_counter)

            raw_dat[donor][tp], nt_dat = fxn.read_cdr3s_in(fxn.conv_data_dir + fl)

            gini = fxn.gini(list(raw_dat[donor][tp].values()))
            shannon = stats.entropy(list(raw_dat[donor][tp].values()), base=2)

            # Then calculate all the actual metrics!
            day_vals.loc[donor, 'd' + str(timepoint)] = day
            day_vals.loc[donor, 'g' + str(timepoint)] = gini
            day_vals.loc[donor, 's' + str(timepoint)] = shannon

            temp_results = [donor, timepoint,
                            meta.loc[donor, 'diagnosis'], meta.loc[donor, 'xRT-type'],
                            len(list(raw_dat[donor][tp].values())), sum(raw_dat[donor][tp].values()),
                            gini, shannon,
                            float(meta.loc[donor, 'TP' + str(timepoint) + '-ALC-K/uL']),
                            day]

            results.append(temp_results)

        # Get together the all-important plotting dataframe
        results = pd.DataFrame(results)
        results = results.rename(index=str, columns=dict(list(zip(list(range(len(results_headers))), results_headers))))

        results.to_pickle(fxn.conv_data_dir + today + '-tcr-results.pkl')
        fxn.save_pickle(fxn.conv_data_dir + today + '-raw-tcr-dict.pkl', raw_dat)

    print("Plotting TCR data ...")

    for val in results_headers[4:]:

        # Violin plots
        save_name = plot_dir + 'violinplot-' + val + '.svg'
        fxn.plot_3tp_signif(results, val, save_name, '', 'diversity')

    # Plot an overview of diversity scores taking the actual dates of the time points into account
    day_plotting = coll.defaultdict(list)
    deltas = []

    # Collect the necessary data from the results database
    for donor in list(meta.index):
        temp_day = results.loc[results['Donor'] == donor, 'Day']
        temp_gini = results.loc[results['Donor'] == donor, 'Gini']
        temp_shannon = results.loc[results['Donor'] == donor, 'Shannon']
        temp_alc = results.loc[results['Donor'] == donor, 'ALC']
        if results.loc[results['Donor'] == donor, 'Treatment'].iloc[0] == 'Photons':
            typ = 'photon'
        else:
            typ = 'proton'

        day_plotting['Gini-' + typ].append([[temp_day.iloc[0], temp_day.iloc[1]],
                                            [temp_gini.iloc[0], temp_gini.iloc[1]]])
        day_plotting['Gini-' + typ].append([[temp_day.iloc[1], temp_day.iloc[2]],
                                            [temp_gini.iloc[1], temp_gini.iloc[2]]])
        day_plotting['Shannon-' + typ].append([[temp_day.iloc[0], temp_day.iloc[1]],
                                               [temp_shannon.iloc[0], temp_shannon.iloc[1]]])
        day_plotting['Shannon-' + typ].append([[temp_day.iloc[1], temp_day.iloc[2]],
                                               [temp_shannon.iloc[1], temp_shannon.iloc[2]]])
        day_plotting['ALC-' + typ].append([[temp_day.iloc[0], temp_day.iloc[1]],
                                           [temp_alc.iloc[0], temp_alc.iloc[1]]])
        day_plotting['ALC-' + typ].append([[temp_day.iloc[1], temp_day.iloc[2]],
                                           [temp_alc.iloc[1], temp_alc.iloc[2]]])

        # Also collect the *changes* in the various metrics
        deltas.append([donor, '1-2', results.loc[results['Donor'] == donor, 'Treatment'].iloc[0],
                       temp_gini.iloc[1] - temp_gini.iloc[0], temp_shannon.iloc[1] - temp_shannon.iloc[0],
                       temp_alc.iloc[1] - temp_alc.iloc[0]])
        deltas.append([donor, '2-3', results.loc[results['Donor'] == donor, 'Treatment'].iloc[0],
                       temp_gini.iloc[2] - temp_gini.iloc[1], temp_shannon.iloc[2] - temp_shannon.iloc[1],
                       temp_alc.iloc[2] - temp_alc.iloc[1]])

        # Calculate delta of first to third sample
        deltas.append([donor, '1-3', results.loc[results['Donor'] == donor, 'Treatment'].iloc[0],
                       temp_gini.iloc[2] - temp_gini.iloc[0], temp_shannon.iloc[2] - temp_shannon.iloc[0],
                       temp_alc.iloc[2] - temp_alc.iloc[0]])

    # Then plot both the relative date/diversity score figures
    for div in ['Gini', 'Shannon', 'ALC']:
        save_name = plot_dir + 'time-specific-' + div + '.svg'
        fig = plt.figure(figsize=(11, 7))
        ax = fig.add_subplot(111)

        sns.scatterplot(x='Day', y=div, hue='Treatment', style='Diagnosis', data=results, s=100, legend='brief')
        if div == 'Shannon':
            plt.legend(loc='lower right')

        for p in day_plotting[div + '-photon']:
            plt.plot(p[0], p[1], color='blue')

        for p in day_plotting[div + '-proton']:
            plt.plot(p[0], p[1], color='orange')

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles[1:3] + handles[4:], labels=labels[1:3] + labels[4:])

        sns.despine(right=True, top=True)
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        plt.close()

    # Plot the whole repertoire (unsampled) intra-donor overlaps (via Jaccard)
    jaccards = []
    jaccards_headers = ['Donor', 'Treatment', 'Transition', 'Jaccard']
    for donor in list(raw_dat.keys()):
        dates = list(raw_dat[donor].keys())
        dates.sort()
        for transition in [(0, 1), (1, 2), (0, 2)]:
            jaccards.append([donor, meta.loc[donor]['xRT-type'],
                             '-'.join([str(x + 1) for x in transition]),
                             fxn.jaccard(list(raw_dat[donor][dates[transition[0]]].keys()),
                                         list(raw_dat[donor][dates[transition[1]]].keys()))])

    jaccards = fxn.list_to_df(jaccards, jaccards_headers, False)
    fxn.plot_3tp_signif(jaccards, 'Jaccard', plot_dir + 'jaccards.svg', '', 'jaccard')

    # again for general legends
    fig = plt.figure(figsize=(1, 1))
    ax = sns.barplot(x='Transition', y='Jaccard', hue='Treatment', data=jaccards)
    sns.despine(right=True, top=True, bottom=True, left=True)
    plt.ylabel('')
    plt.xlabel('')
    # plt.yticklabel('')
    ax.axis('off')
    plt.xlim(10, 20)
    plt.ylim(10, 20)
    ax.legend(loc='center left', bbox_to_anchor=(1.1, 1.05))
    save_name = plot_dir + 'barplot-LEGEND.svg'
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close()

    # Comparisons of changes in diversity, using the delta data collected above
    deltas_headers = ['Donor', 'Transition', 'Treatment', 'ΔGini', 'ΔShannon', 'ΔALC']
    deltas = pd.DataFrame(deltas)
    deltas = deltas.rename(index=str, columns=dict(list(zip(list(range(len(deltas_headers))), deltas_headers))))
    deltas.to_pickle(fxn.conv_data_dir + today + '-deltas.pkl')

    summary_params = {'ΔGini': 'gini', 'ΔShannon': 'shannon', 'ΔALC': 'alc'}

    summary_ys = {'ΔGini': 0.52, 'ΔShannon': 3.5, 'ΔALC': 2.2}
    summary_ys_off = {'ΔGini': 0.02, 'ΔShannon': .4, 'ΔALC': 0.2}
    x_off = 0.05
    
    transitions = ['1-2', '2-3']
    deltas = deltas.loc[deltas['Transition'].isin(transitions)]

    for div in summary_params:
        save_name = plot_dir + 'change-' + div + '.svg'
        fig = plt.figure(figsize=(11, 7))
        ax = fig.add_subplot(111)
        sns.catplot(x='Transition', y=div, hue='Treatment', cut=0, kind='violin',
                    inner='stick', split=True, data=deltas, legend=False, density_norm='count')
        for trans in transitions:
            vars()['p' + trans.replace('-', '')] = \
                fxn.mwu(deltas.loc[(deltas.Transition == trans) & (deltas.Treatment == 'Photons')][div],
                        deltas.loc[(deltas.Transition == trans) & (deltas.Treatment == 'Protons')][div])
        fxn.plot_signif_lines(0 - x_off, 0 + x_off, summary_ys[div], p12, 'black')
        fxn.plot_signif_lines(1 - x_off, 1 + x_off, summary_ys[div], p23, 'black')
        for treat in ['Photons', 'Protons']:
            vars()['p_' + treat] = \
                fxn.wc(deltas.loc[(deltas.Transition == '1-2') & (deltas.Treatment == treat)][div],
                       deltas.loc[(deltas.Transition == '2-3') & (deltas.Treatment == treat)][div])
        fxn.plot_signif_lines(0 - x_off, 1 - x_off, summary_ys[div] + summary_ys_off[div], p23, fxn.blue)
        fxn.plot_signif_lines(0 + x_off, 1 + x_off, summary_ys[div] + (2 * summary_ys_off[div]), p23, fxn.orange)
        plt.subplots_adjust(bottom=0.13, left=0.20, right=0.98, top=0.9)
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        plt.close()

    # Showing the relationship between the donor-specific changes
    for div in summary_params:

        vars()['d' + summary_params[div] + '_wide'] = deltas.pivot(index='Donor', columns='Transition', values=div)
        vars()['d' + summary_params[div] + '_wide']['Treatment'] = meta['xRT-type']
        vars()['d' + summary_params[div] + '_wide']['Diagnosis'] = meta['diagnosis']

        save_name = plot_dir + 'donor-specific-change-' + summary_params[div] + '.svg'
        fig = plt.figure(figsize=(5, 5))
        # fig = plt.figure()
        ax = fig.add_subplot(111)
        sns.scatterplot(x='1-2', y='2-3', hue='Treatment', style='Diagnosis',
                        data=vars()['d' + summary_params[div] + '_wide'], s=100, legend='brief')
        ax.axvline(0, color='black', alpha=.5, ls='dotted')
        ax.axhline(0, color='black', alpha=.5, ls='dotted')
        plt.xlabel(div + '(TP1-2, baseline→nadir)')
        plt.ylabel(div + '(TP2-3, nadir→recovery)')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.subplots_adjust(bottom=0.13, left=0.20, right=0.98, top=0.9)
        sns.despine(right=True, top=True)
        plt.savefig(save_name, dpi=300)  #, bbox_inches='tight')
        plt.close()

    # and again for separate legends
    fig = plt.figure(figsize=(1, 1))
    ax = sns.scatterplot(x='1-2', y='2-3', hue='Treatment', style='Diagnosis',
                         data=vars()['d' + summary_params[div] + '_wide'], s=100, legend='brief')

    sns.despine(right=True, top=True, bottom=True, left=True)
    plt.ylabel('')
    plt.xlabel('')
    # plt.yticklabel('')
    ax.axis('off')
    plt.xlim(10, 20)
    plt.ylim(10, 20)
    ax.legend(loc='center left', bbox_to_anchor=(1.1, 1.05))
    save_name = plot_dir + 'donor-specific-change-' + summary_params[div] + '-LEGEND.svg'
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close()
