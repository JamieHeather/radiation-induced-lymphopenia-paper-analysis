#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
iterative-subsample-tcrs.py

Perform similar analysis to some of those in analyse-tcrs.py, but on iteratively subsampled datasets

Particularly important for the sharing and diversity experiments
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import os
import sys
import collections as coll
import functions as fxn
import random


__version__ = '0.2.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

# Sort plotting parameters
# sns.set(style="darkgrid")
plt.rcParams.update({'font.size': 20, 'font.sans-serif': 'Arial'})


def sample_repertoire(counter_to_sample, type_of_sampling, sampling_number):
    """
    :param counter_to_sample: raw_dat dict (with date/donor specified)
    :param type_of_sampling: unique or total, i.e. ignoring or factoring in TCR abundance respectively
    :param sampling_number: sample_to
    :return: An appropriately randomly subsampled counter
    """
    if type_of_sampling == 'total':
        return coll.Counter(random.sample(fxn.flatten_list([[x] * counter_to_sample[x] for x in counter_to_sample]),
                                          sampling_number))

    elif type_of_sampling == 'unique':
        return coll.Counter(random.sample(list(counter_to_sample.keys()), sampling_number))

    else:
        print("Error: unknown sampling type - ", sampling_type)


if __name__ == "__main__":

    fxn.check_scripts_dir()
    today = fxn.today()
    plot_dir = fxn.plot_dir('subsampled-tcrs')

    transitions_index = [(0, 1), (1, 2), (0, 2)]
    transitions_tp = [(1, 2), (2, 3), (1, 3)]

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
        print("Failed to find previously sorted raw data files - please run analysis-tcrs.py first")
        sys.exit()

    meta = fxn.get_metadata()

    num_iterations = 100  # TODO change this number manually to alter number of iterations

    sample_to_dict = {  #'unique': [500, 1000, 2000, 5000],  # TODO If desired users can sample unique TCRs too
                      'total': [1000, 10000, 4000]}
    print("Reading in and subsampling data...")
    for sampling_type in sample_to_dict:

        for sample_to in sample_to_dict[sampling_type]:

            print("Sampling to", str(sample_to), 'total TCRs...')

            # Ensure all donors have sufficient TCRs at each timepoint
            donors = list(raw_dat.keys())
            donors.sort()
            to_remove = []
            for d in donors:
                if sampling_type == 'total':
                    if not all(i >= sample_to for i in [sum(raw_dat[d][x].values()) for x in raw_dat[d]]):
                        to_remove.append(d)
                elif sampling_type == 'unique':
                    if not all(i >= sample_to for i in [len(list(raw_dat[d][x].values())) for x in raw_dat[d]]):
                        to_remove.append(d)

            for d in to_remove:
                donors.pop(donors.index(d))

            subsampled = []
            subsampled_headers = ['Donor', 'Treatment', 'Iteration', 'Timepoint', 'Gini', 'Shannon']

            jaccards = []
            jaccards_headers = ['Donor', 'Treatment', 'Iteration', 'Transition', 'Jaccard']

            for iteration in range(num_iterations):
                if (iteration + 1) % 10 == 0:
                    print(iteration + 1)

                # Store each timepoint's subsampled repertoire for Jaccard calculation
                subsampled_donor = {}

                # Then loop through each donor storing each date's diversity metrics...
                for donor in donors:
                    dates = list(raw_dat[donor].keys())
                    dates.sort()
                    for d in range(len(dates)):
                        tp = d + 1
                        subsampled_donor[tp] = sample_repertoire(raw_dat[donor][dates[d]], sampling_type, sample_to)

                        subsampled.append([donor, meta.loc[donor]['xRT-type'], iteration, tp,
                                           fxn.gini(list(subsampled_donor[tp].values())),
                                           stats.entropy(list(subsampled_donor[tp].values()), base=2)])

                    # ... then the Jaccards of each of the transitions
                    for transition in transitions_tp:
                        jaccards.append([donor, meta.loc[donor]['xRT-type'], iteration,
                                         '-'.join([str(x) for x in transition]),
                                         fxn.jaccard(list(subsampled_donor[transition[0]].keys()),
                                                     list(subsampled_donor[transition[1]].keys()))])

            subsampled = fxn.list_to_df(subsampled, subsampled_headers, False)
            jaccards = fxn.list_to_df(jaccards, jaccards_headers, False)

            # Plot subsampled diversities (total and averaged)
            print("Calculating averages...")
            avgs = []
            for d in donors:
                for tp in range(1, 4):
                    subset = subsampled.loc[(subsampled['Donor'] == d) & (subsampled['Timepoint'] == tp)]
                    out_list = [d, subset['Treatment'].iloc[0], tp,
                                np.mean(subset['Gini']), np.std(subset['Gini']),
                                np.mean(subset['Shannon']), np.std(subset['Shannon'])]
                    avgs.append(out_list)

            avgs = fxn.list_to_df(avgs, ['Donor', 'Treatment', 'Timepoint', 'Gini', 'Gini SD',
                                         'Shannon', 'Shannon SD'], False)

            avgs[['ALC', 'Day']] = results[['ALC', 'Day']]
            avgs.index = pd.to_numeric(avgs.index, errors='coerce')
            avgs = avgs.sort_index()

            for val in ['Gini', 'Shannon']:
                save_name = plot_dir + sampling_type + '-' + 'violinplot-' + str(sample_to) + '-all-' + val + '.svg'
                fxn.plot_3tp_signif(subsampled, val, save_name, '', 'diversity')

                save_name = plot_dir + sampling_type + '-' + 'violinplot-' + str(sample_to) + '-avg-' + val + '.svg'
                fxn.plot_3tp_signif(avgs, val, save_name, '', 'diversity')

            # And jaccards (total and averaged)
            save_name = plot_dir + sampling_type + '-' + 'violinplot-' + str(sample_to) + '-all-jaccards.svg'
            fxn.plot_3tp_signif(jaccards, 'Jaccard', save_name, '', 'jaccard')

            avg_jacs = []
            for donor in donors:
                for transition in transitions_tp:
                    trans = '-'.join([str(x) for x in transition])
                    tmp_data = jaccards.loc[(jaccards['Donor'] == donor) & (jaccards['Transition'] == trans)]
                    avg_jacs.append([donor, trans, tmp_data.iloc[0]['Treatment'], np.mean(tmp_data['Jaccard'])])

            avg_jacs = fxn.list_to_df(avg_jacs, ['Donor', 'Transition', 'Treatment', 'Jaccard'], False)
            save_name = plot_dir + sampling_type + '-' + 'violinplot-' + str(sample_to) + '-avg-jaccards.svg'
            fxn.plot_3tp_signif(avg_jacs, 'Jaccard', save_name, '', 'jaccard')

    # Plot an overview of diversity scores taking the actual dates of the time points into account
    day_plotting = coll.defaultdict(list)
    deltas = []

    print("Calculating deltas...")
    # Collect the necessary data from the avgs database
    for donor in list(meta.index):
        temp_day = avgs.loc[avgs['Donor'] == donor, 'Day']
        temp_gini = avgs.loc[avgs['Donor'] == donor, 'Gini']
        temp_shannon = avgs.loc[avgs['Donor'] == donor, 'Shannon']
        temp_alc = avgs.loc[avgs['Donor'] == donor, 'ALC']
        if avgs.loc[avgs['Donor'] == donor, 'Treatment'].iloc[0] == 'Photons':
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
        deltas.append([donor, '1-2', avgs.loc[avgs['Donor'] == donor, 'Treatment'].iloc[0],
                       temp_gini.iloc[1] - temp_gini.iloc[0],
                       temp_shannon.iloc[1] - temp_shannon.iloc[0],
                       temp_alc.iloc[1] - temp_alc.iloc[0]])
        deltas.append([donor, '2-3', avgs.loc[avgs['Donor'] == donor, 'Treatment'].iloc[0],
                       temp_gini.iloc[2] - temp_gini.iloc[1],
                       temp_shannon.iloc[2] - temp_shannon.iloc[1],
                       temp_alc.iloc[2] - temp_alc.iloc[1]])
        deltas.append([donor, '1-3', avgs.loc[avgs['Donor'] == donor, 'Treatment'].iloc[0],
                       temp_gini.iloc[2] - temp_gini.iloc[0],
                       temp_shannon.iloc[2] - temp_shannon.iloc[0],
                       temp_alc.iloc[2] - temp_alc.iloc[0]])

    # Plot the whole repertoire (unsampled) intra-donor overlaps (via Jaccard)
    jaccards = []
    jaccards_headers = ['Donor', 'Treatment', 'Transition', 'Jaccard']
    for donor in list(raw_dat.keys()):
        dates = list(raw_dat[donor].keys())
        dates.sort()
        for transition in transitions_index:
            jaccards.append([donor, meta.loc[donor]['xRT-type'],
                             '-'.join([str(x + 1) for x in transition]),
                             fxn.jaccard(list(raw_dat[donor][dates[transition[0]]].keys()),
                                         list(raw_dat[donor][dates[transition[1]]].keys()))])

    jaccards = fxn.list_to_df(jaccards, jaccards_headers, False)
    fxn.plot_3tp_signif(jaccards, 'Jaccard', plot_dir + 'jaccards.svg', '', 'jaccard')

    # Comparisons of changes in diversity, using the delta data collected above
    deltas_headers = ['Donor', 'Transition', 'Treatment', 'ΔGini', 'ΔShannon', 'ΔALC']
    deltas = pd.DataFrame(deltas)
    deltas = deltas.rename(index=str, columns=dict(list(zip(list(range(len(deltas_headers))), deltas_headers))))
    deltas.to_pickle(fxn.conv_data_dir + today + '-deltas.pkl')

    # Plotting whether or not diversity metrics impact upon whether a patient's tumor recurred
    recurred = {'Y': list(meta.loc[meta['Recurred'] == 1].index),
                'N': list(meta.loc[meta['Recurred'] == 0].index)}

    infections = {'N': list(meta[meta['First-infection'].isnull()].index),
                  'Y': list(meta[meta['First-infection'].notnull()].index)}

    prior_chemo = {'5FU': list(meta.loc[meta['Prior-chemo-group'] == '5FU'].index),
                   'Platin': list(meta.loc[meta['Prior-chemo-group'] == 'Platin'].index)}

    cea_fields = [x for x in meta if 'CEA' in x and 'date' not in x]
    cea_dat = meta[cea_fields]
    high_cea = {'Y': [], 'N': []}
    for donor in cea_dat.index:
        if len([x for x in cea_dat.loc[donor] if x > 50]) > 0:
            high_cea['Y'].append(donor)
        else:
            high_cea['N'].append(donor)

    recurrers = []  # Store which photon patients saw very large expansions post-nadir
    infect = []
    cea_pos = []
    chemos = []
    for row in avgs.index:
        if avgs.loc[row]['Donor'] in recurred['Y']:
            recurrers.append('Y')
        elif avgs.loc[row]['Donor'] in recurred['N']:
            recurrers.append('N')
        else:
            recurrers.append('')

        if avgs.loc[row]['Donor'] in infections['Y']:
            infect.append('Y')
        else:
            infect.append('N')

        if avgs.loc[row]['Donor'] in high_cea['Y']:
            cea_pos.append('Y')
        else:
            cea_pos.append('N')

        if avgs.loc[row]['Donor'] in prior_chemo['5FU']:
            chemos.append('5FU')
        elif avgs.loc[row]['Donor'] in prior_chemo['Platin']:
            chemos.append('Platin')
        else:
            chemos.append('')

    avgs['Recurred'] = recurrers
    avgs['Infections'] = infect
    avgs['CEA high'] = cea_pos
    avgs['Prior Chemo'] = chemos

    recurrers = []  # Store which photon patients saw very large expansions post-nadir
    infect = []
    cea_pos = []
    chemos = []
    for row in deltas.index:

        if deltas.loc[row]['Donor'] in recurred['Y']:
            recurrers.append('Y')
        elif deltas.loc[row]['Donor'] in recurred['N']:
            recurrers.append('N')
        else:
            recurrers.append('')

        if deltas.loc[row]['Donor'] in infections['Y']:
            infect.append('Y')
        else:
            infect.append('N')

        if deltas.loc[row]['Donor'] in high_cea['Y']:
            cea_pos.append('Y')
        else:
            cea_pos.append('N')

        if deltas.loc[row]['Donor'] in prior_chemo['5FU']:
            chemos.append('5FU')
        elif deltas.loc[row]['Donor'] in prior_chemo['Platin']:
            chemos.append('Platin')
        else:
            chemos.append('')

    deltas['Recurred'] = recurrers
    deltas['Infections'] = infect
    deltas['CEA high'] = cea_pos
    deltas['Prior Chemo'] = chemos

    # TODO tidy these up and get them saving somewhere
    feat_orders = {'Recurred': ['Y', 'N'],
                   'Infections': ['Y', 'N'],
                   'CEA high': ['Y', 'N'],
                   'Prior Chemo': ['5FU', 'Platin']}

    clin_plot_dir = fxn.make_check_dir(plot_dir + 'clinical-comparisons_'
                                       + str(sample_to) + '_' + str(num_iterations) + '/')
    colors = ['purple', 'azure']
    sns.set_palette(sns.xkcd_palette(colors))

    for feature in ['Recurred', 'Infections', 'CEA high', 'Prior Chemo']:
        # AVGS
        for metric in ['Gini', 'Shannon']:
            print('\n------------------\n', feature, metric)
            tmp_ord = feat_orders[feature]
            g = sns.catplot(x=feature, y=metric, data=avgs, kind="violin", order=tmp_ord, col='Timepoint', cut=0,
                            inner='stick', hue=feature, density_norm='count')
            for ax, tp in zip(g.axes.flat, [1, 2, 3]):
                tmp_dat1 = avgs[(avgs[feature] == tmp_ord[0]) & (avgs['Timepoint'] == tp)]
                tmp_dat2 = avgs[(avgs[feature] == tmp_ord[1]) & (avgs['Timepoint'] == tp)]
                print(fxn.mwu(tmp_dat1[metric], tmp_dat2[metric]))
            save_name = clin_plot_dir + 'viol-' + feature + '-' + metric + '.svg'
            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            plt.close()
            tmp_dat1 = avgs[(avgs[feature] == tmp_ord[0])]
            tmp_dat2 = avgs[(avgs[feature] == tmp_ord[1])]
            print('\t\t' + str(fxn.mwu(tmp_dat1[metric], tmp_dat2[metric])))
            save_name = clin_plot_dir + 'viol-by-tp' + feature + '-' + metric + '.svg'
        # DELTAS
        for metric in ['\u0394Gini', '\u0394Shannon']:
            print('\n------------------\n', feature, metric)
            tmp_ord = feat_orders[feature]
            g = sns.catplot(x=feature, y=metric, data=deltas[deltas['Transition'] != '1-3'], kind="violin",
                            order=tmp_ord, col='Transition', cut=0, inner='stick', hue=feature, density_norm='count')
            for ax, tp in zip(g.axes.flat, ['1-2', '2-3']):  #, '1-3']):
                tmp_dat1 = deltas[(deltas[feature] == tmp_ord[0]) & (deltas['Transition'] == tp)]
                tmp_dat2 = deltas[(deltas[feature] == tmp_ord[1]) & (deltas['Transition'] == tp)]
                print(fxn.mwu(tmp_dat1[metric], tmp_dat2[metric]))
            save_name = clin_plot_dir + 'delta-viol-' + feature + '-' + metric.replace('\u0394', 'd') + '.svg'
            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            plt.close()
            tmp_dat1 = deltas[(deltas[feature] == tmp_ord[0])]
            tmp_dat2 = deltas[(deltas[feature] == tmp_ord[1])]
            print('\t\t' + str(fxn.mwu(tmp_dat1[metric], tmp_dat2[metric])))

    for val in ['Infections', 'Recurred']:
        save_name = plot_dir + val.lower() + '.svg'
        fig = plt.figure(figsize=(10, 5))
        sns.catplot(x=val, y='Gini',  data=avgs, order=['Y', 'N'], kind="violin", col='Timepoint', cut=0,
                    inner='stick', palette=fxn.alt_cols, hue=val, density_norm='count')
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        plt.close()

    # Save the final averaged 'results' dict out for downstream use
    remaining_cols = [x for x in results if x not in avgs]
    if (list(avgs['Donor']) == list(results['Donor'])) and (list(avgs['Timepoint']) == list(results['Timepoint'])):
        for rcol in remaining_cols:
            avgs[rcol] = list(results[rcol])
        avgs.to_pickle(fxn.conv_data_dir + today + '-' + str(sample_to) + '-sized-avg-tcr-results.pkl')
