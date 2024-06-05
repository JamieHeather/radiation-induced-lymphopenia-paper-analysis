#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
analyse-flow.py

Finish off the other half of the analysis, processing the flow cytometry data.

"""


import functions as fxn
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
__version__ = '0.2.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

# Sort plotting parameters
plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial'})

if __name__ == "__main__":

    # Sort out the necessary directories, checking we're in the right one
    fxn.check_scripts_dir()
    today = fxn.today()
    plot_dir = fxn.plot_dir('flow-analysis')
    data_file = [x for x in os.listdir(fxn.raw_data_dir) if x == 'flow-data.csv'][0]

    # First get the previous data
    print("Reading in previous data ...")
    results = fxn.open_previous_data(fxn.conv_data_dir, 'tcr-results', 'pkl')

    # Get the data
    line_count = 0
    fields_to_add = ['Treatment', 'Gini', 'Shannon', 'ALC', 'Day',
                     'Diagnosis', 'Unique-TCRs', 'Total-TCRs']
    flowdat = []
    with open(fxn.raw_data_dir + data_file, 'r') as in_file:
        for line in in_file:
            bits = line.rstrip().split(',')
            if line_count == 0:
                bits[0] = 'Donor'
                headers = bits
            else:
                # Convert from string to relevant data types
                bits[2] = int(bits[2])
                bits[3:] = [float(x) for x in bits[3:]]
                # Then pull out the relevant fields
                tmp_results = results.loc[(results.Donor == bits[0]) & (results.Timepoint == bits[2])]
                # tmp_vals = [tmp_results.iloc[[0]][x][0] for x in fields_to_add]
                tmp_vals = [tmp_results.iloc[0][x] for x in fields_to_add]
                flowdat.append(bits + tmp_vals)
            line_count += 1

    flowdat = pd.DataFrame(flowdat)
    flowdat = flowdat.rename(index=str, columns=dict(list(zip(list(range(len(headers + fields_to_add))),
                                                              headers + fields_to_add))))

    # Loop through the flow populations and plot
    populations = [x for x in headers if '|' in x or 'cells' in x or 'CD' in x]
    populations.sort()
    percent_plot_dir = fxn.make_check_dir(plot_dir + 'percentages/')
    for population in populations:
        save_name = percent_plot_dir + 'pc-vp-' + population.replace('|', '_').replace(' ', '_') + '.svg'
        fxn.plot_3tp_signif(flowdat, population, save_name, '(%)', 'diversity')

    plt.close('all')

    # Want to convert those into absolute values, going off the ALC
    abs_convert = {
                    'B cells': 'ALC',
                    'CD3- CD19-': 'ALC',
                    'CD3- CD19-|NK cells': 'CD3- CD19-',
                    'T cells': 'ALC',
                    'T cells|NKT cells': 'T cells',
                    'CD4+': 'T cells',
                    'CD4+|Tcm': 'CD4+',
                    'CD4+|Tem': 'CD4+',
                    'CD4+|Temra': 'CD4+',
                    'CD4+|Tn': 'CD4+',
                    'CD8+': 'T cells',
                    'CD8+|Tcm': 'CD8+',
                    'CD8+|Tem': 'CD8+',
                    'CD8+|Temra': 'CD8+',
                    'CD8+|Tn': 'CD8+',
                    'CD4+|Treg': 'CD4+'
                   }

    # Need to ensure we process in the right order (mostly sequential calculations relying on previous values)
    abs_order = ['T cells', 'T cells|NKT cells', 'B cells', 'CD4+', 'CD4+|Tcm', 'CD4+|Tem', 'CD4+|Temra', 'CD4+|Tn',
                 'CD4+|Treg', 'CD8+', 'CD8+|Tcm', 'CD8+|Tem', 'CD8+|Temra', 'CD8+|Tn',
                 'CD3- CD19-', 'CD3- CD19-|NK cells']

    cols = ['Donor', 'Timepoint', 'Diagnosis', 'Treatment', 'ALC']
    absdat = flowdat.loc[:, cols]

    for popn in abs_order:
        print(popn)
        if abs_convert[popn] == 'ALC':
            new_col = flowdat[popn] * (1000 * absdat[abs_convert[popn]])
        else:
            new_col = flowdat[popn] * absdat[abs_convert[popn]]
        # Then stick on to growing dataframe
        new_col = new_col.rename(popn)
        absdat = pd.concat([absdat, new_col], axis=1)
        cols.append(popn)

    # Loop through the flow populations and plot
    abs_plot_dir = fxn.make_check_dir(plot_dir + 'absolute/')
    for population in abs_order:
        print(population)
        save_name = abs_plot_dir + 'abs-vp-' + population.replace('|', '_').replace(' ', '_') + '.svg'
        fxn.plot_3tp_signif(absdat, population, save_name, '(calc\'d cells / μl)', 'diversity', True)

    plt.close('all')

    # Get previous 'delta' data (i.e. changes in diversity/ALC
    deltas_pkls = [x for x in os.listdir(fxn.conv_data_dir) if x.endswith('-deltas.pkl')]
    if len(deltas_pkls) == 0:
        print("Error: cannot find a \'deltas\' pkl!")
        sys.exit()

    deltas_pkls.sort()
    div_deltas = pd.read_pickle(fxn.conv_data_dir + deltas_pkls[-1])

    # Get a sorted list of donors for which we have flow data
    donors = list(set(flowdat['Donor']))
    donors.sort()

    # Then loop through the various population fields and calculate their change for each transition
    deltas = []
    deltas_pc = []
    delta_popn_names = ['Δ' + x for x in abs_order]
    delta_headers = list(div_deltas) + delta_popn_names
    for d in donors:
        print(d)
        tmp_donor_divdat = div_deltas.loc[div_deltas.Donor == d]
        tmp_donor_absdat = absdat.loc[absdat.Donor == d]
        tmp_donor_pc = flowdat.loc[flowdat.Donor == d]

        # Loop through transitions and gather/calculate the appropriate data to act as the base info for new rows
        for t in [(1, 2), (2, 3)]:  #, (1,3)]:
            tmp_divdat = tmp_donor_divdat.loc[tmp_donor_divdat.Transition == str(t[0]) + '-' + str(t[1])]
            out_dat = [tmp_divdat.iloc[[0]][x][0] for x in tmp_divdat]
            out_dat_pc = [tmp_divdat.iloc[[0]][x][0] for x in tmp_divdat]

            # Then calculate the new transitions to add to this
            for p in abs_order:
                v1 = float(tmp_donor_absdat.loc[tmp_donor_absdat.Timepoint == t[0]][p])
                v2 = float(tmp_donor_absdat.loc[tmp_donor_absdat.Timepoint == t[1]][p])
                out_dat.append(v2 - v1)
                v1pc = float(tmp_donor_pc.loc[tmp_donor_pc.Timepoint == t[0]][p])
                v2pc = float(tmp_donor_pc.loc[tmp_donor_pc.Timepoint == t[1]][p])
                out_dat_pc.append(v2pc - v1pc)
            deltas.append(out_dat)
            deltas_pc.append(out_dat_pc)

    deltas = pd.DataFrame(deltas)
    deltas = deltas.rename(index=str, columns=dict(list(zip(list(range(len(delta_headers))), delta_headers))))

    deltas_pc = pd.DataFrame(deltas_pc)
    deltas_pc = deltas_pc.rename(index=str, columns=dict(list(zip(list(range(len(delta_headers))), delta_headers))))

    delta_plot_dir = fxn.make_check_dir(plot_dir + 'deltas/')
    for p in delta_popn_names:
        save_name = delta_plot_dir + 'change-' + p.replace('|', '_').replace(' ', '_') + '.svg'
        fxn.plot_3tp_signif(flowdat, population, save_name, '', 'diversity')
        fig = plt.figure(figsize=(5, 5))
        # ax = fig.add_subplot(111)
        # ax = fig.add_subplot(111)
        p = sns.catplot(x='Transition', y=p, hue='Treatment', cut=0, kind='violin', density_norm='count',
                        inner='stick', split=True, data=deltas, hue_order=['Photons', 'Protons'], legend=False)
        for pax in p.axes.flatten():
            pax.ticklabel_format(style='scientific', scilimits=(0, 0), axis='y', useMathText=True)
        plt.subplots_adjust(bottom=0.13, left=0.20, right=0.98, top=0.9)
        plt.savefig(save_name, dpi=300)  #, bbox_inches='tight')
        plt.close()

    absdat[list(flowdat)[-7:]] = flowdat[list(flowdat)[-7:]]

    # Plot the relation between each population and the various population level metrics
    for d in list(div_deltas)[3:]:

        # Against deltas
        # delta_cols_plot_dir = fxn.make_check_dir(plot_dir + 'comps-delta-cols/')
        delta_cols_pc_plot_dir = fxn.make_check_dir(plot_dir + 'comps-delta-cols-pc/')
        delta_v_alc_dir = fxn.make_check_dir(plot_dir + 'comps-delta-v-alc/')

        for p in delta_popn_names:

            # Deltas split out by transition
            save_name = delta_cols_pc_plot_dir + 'delta-col-pc-' + d + '-' + p.replace('|', '_') + '.svg'
            mins, maxs = fxn.get_limits(deltas_pc)
            mins['\u0394Gini'] = -0.21
            maxs['\u0394Gini'] = 0.6
            trans_order = ['1-2', '2-3']  #, '1-3']

            g = sns.lmplot(data=deltas_pc, x=p, y=d, hue='Treatment', col='Transition', scatter_kws={'clip_on': False},
                           hue_order=['Photons', 'Protons'], col_order=trans_order, truncate=False)
            plt.xlim(mins[p], maxs[p])
            plt.ylim(mins[d], maxs[d])
            for ax, trans in zip(g.axes.flat, trans_order):
                ax.axvline(0, color='black', alpha=.5, ls='dotted')
                ax.axhline(0, color='black', alpha=.5, ls='dotted')
                subset = deltas_pc[deltas_pc['Transition'] == trans]
                ph = subset[subset['Treatment'] == 'Photons']
                pr = subset[subset['Treatment'] == 'Protons']
                ax.text(0.01, 1.02, fxn.get_regression_str(ph[d], ph[p]), color=fxn.blue,
                        transform=ax.transAxes, fontsize=10, fontweight='bold')
                ax.text(0.99, 1.02, fxn.get_regression_str(pr[d], pr[p]), color=fxn.orange,
                        transform=ax.transAxes, fontsize=10, horizontalalignment='right', fontweight='bold')

            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            plt.close()

            # save_name = delta_cols_plot_dir + 'delta-col-pc-' + d + '-' + p.replace('|', '_') + '.svg'
            # mins, maxs = fxn.get_limits(deltas)
            # trans_order = ['1-2', '2-3']  # , '1-3'
            #
            # g = sns.lmplot(data=deltas_pc, x=p, y=d, hue='Treatment', col='Transition',
            #                hue_order=['Photons', 'Protons'], col_order=trans_order)
            # plt.xlim(mins[p], maxs[p])
            # plt.ylim(mins[d], maxs[d])
            # for ax, trans in zip(g.axes.flat, trans_order):
            #     ax.axvline(0, color='r', alpha=.5)
            #     ax.axhline(0, color='r', alpha=.5)
            #     subset = deltas_pc[deltas_pc['Transition'] == trans]
            #     ph = subset[subset['Treatment'] == 'Photons']
            #     pr = subset[subset['Treatment'] == 'Protons']
            #     ax.text(0.01, 1.02, fxn.get_regression_str(ph[d], ph[p]), color=fxn.blue,
            #             transform=ax.transAxes, fontsize=10)
            #     ax.text(0.7, 1.02, fxn.get_regression_str(pr[d], pr[p]), color=fxn.orange,
            #             transform=ax.transAxes, fontsize=10)
            #
            # plt.savefig(save_name, dpi=300, bbox_inches='tight')
            # plt.close()

        # Specifically check the relationship between change in diversity metrics and change in ALC
        if d != '\u0394ALC':
            save_name = delta_v_alc_dir + d + '.svg'
            g = sns.lmplot(data=deltas_pc, x='\u0394ALC', y=d, hue='Treatment', col='Transition',
                           hue_order=['Photons', 'Protons'], col_order=trans_order, truncate=False)
            plt.xlim(mins['\u0394ALC'], maxs['\u0394ALC'])
            plt.ylim(mins[d], maxs[d])
            for ax, trans in zip(g.axes.flat, trans_order):
                ax.axvline(0, color='black', alpha=.5, ls='dotted')
                ax.axhline(0, color='black', alpha=.5, ls='dotted')
                subset = deltas_pc[deltas_pc['Transition'] == trans]
                ph = subset[subset['Treatment'] == 'Photons']
                pr = subset[subset['Treatment'] == 'Protons']
                ax.text(0.01, 1.02, fxn.get_regression_str(ph[d], ph[p]), color=fxn.blue,
                        transform=ax.transAxes, fontsize=10, fontweight='bold')
                ax.text(0.99, 1.02, fxn.get_regression_str(pr[d], pr[p]), color=fxn.orange,
                        transform=ax.transAxes, fontsize=10, horizontalalignment='right', fontweight='bold')

            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            plt.close()

    percent_div_plot_dir = fxn.make_check_dir(plot_dir + 'percent-V-div/')
    percent_baseline_plot_dir = fxn.make_check_dir(plot_dir + 'percent-V-div-baseline-combined/')

    maxes = {'Gini': 1, 'Shannon': fxn.round_up(max(flowdat['Shannon'])), 'ALC':  fxn.round_up(max(flowdat['ALC']))}
    mins = {'Gini': 0, 'Shannon': 5, 'ALC':  fxn.round_up(min(flowdat['ALC']))}

    for tp in range(1, 4):
        tp_temp = flowdat[flowdat['Timepoint'] == tp]
        for flow_pop in [x for x in flowdat if '|' in x or 'cells' in x or x.startswith('CD')]:
            for metric in ['Gini', 'Shannon', 'ALC']:
                save_name = percent_div_plot_dir + str(tp) + '-' + metric + '-' + flow_pop.replace('|', '_') + '.svg'
                g = sns.lmplot(data=tp_temp, x=flow_pop, y=metric, scatter_kws={'clip_on': False},
                               hue='Treatment', hue_order=['Photons', 'Protons'], truncate=False)
                ax = g.axes[0, 0]
                ph = tp_temp[tp_temp['Treatment'] == 'Photons']
                pr = tp_temp[tp_temp['Treatment'] == 'Protons']
                ax.text(0.05, .95, fxn.get_regression_str(ph[metric], ph[flow_pop]), color=fxn.blue,
                        transform=ax.transAxes, fontsize=12, fontweight='bold')
                ax.text(0.99, .95, fxn.get_regression_str(pr[metric], pr[flow_pop]), color=fxn.orange,
                        transform=ax.transAxes, fontsize=12, horizontalalignment='right', fontweight='bold')
                plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                plt.gca().set_xlim(xmin=0)
                plt.ylim(0, maxes[metric])
                plt.savefig(save_name, dpi=300, bbox_inches='tight')
                plt.close()

    for flow_pop in [x for x in flowdat if '|' in x or 'cells' in x or x.startswith('CD')]:
        for metric in ['Gini', 'Shannon', 'ALC']:
            # Plot all separately
            save_name = percent_div_plot_dir + 'all-' + metric + '-' + flow_pop.replace('|', '_') + '.svg'
            sns.set_palette("deep")
            g = sns.lmplot(data=flowdat, x=flow_pop, y=metric, hue='Treatment', col='Timepoint',
                           hue_order=['Photons', 'Protons'], truncate=False)

            for ax, tp in zip(g.axes.flat, list(range(1, 4))):
                subset = flowdat[flowdat['Timepoint'] == tp]
                ph = subset[subset['Treatment'] == 'Photons']
                pr = subset[subset['Treatment'] == 'Protons']
                ax.text(0.01, 1.02, fxn.get_regression_str(ph[flow_pop], ph[metric]), color=fxn.blue,
                        transform=ax.transAxes, fontsize=10, fontweight='bold')
                ax.text(0.99, 1.02, fxn.get_regression_str(pr[flow_pop], pr[metric]), color=fxn.orange,
                        transform=ax.transAxes, fontsize=10, horizontalalignment='right', fontweight='bold')
                # print('\t'.join([flow_pop, metric,
                #                  'photon', str(list(ph[flow_pop])), str(list(ph[metric])),
                #                  'proton', str(list(pr[flow_pop])), str(list(pr[metric]))]))

            plt.gca().set_xlim(xmin=0)
            plt.ylim(0, maxes[metric])
            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            plt.close()

            # And whole populations at baseline combined together
            save_name = percent_baseline_plot_dir + metric + '-' + flow_pop.replace('|', '_') + '.svg'
            tp_temp = flowdat[flowdat['Timepoint'] == 1]
            sns.set_palette(sns.dark_palette("purple"))
            g = sns.lmplot(data=tp_temp, x=flow_pop, y=metric, truncate=False)
            ax = g.axes[0, 0]
            ax.text(0.05, .95, fxn.get_regression_str(tp_temp[metric], tp_temp[flow_pop]),
                    transform=ax.transAxes, fontsize=12)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.gca().set_xlim(xmin=0)
            plt.ylim(mins[metric], maxes[metric])
            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            plt.close()

    plt.close('all')

    # Save the dataframes for later or separate use
    fxn.save_csv(flowdat, fxn.conv_data_dir, 'flow-data')
    fxn.save_csv(deltas, fxn.conv_data_dir, 'flow-deltas')
