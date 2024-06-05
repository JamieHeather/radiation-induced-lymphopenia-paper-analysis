#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
analyse-metadata.py

Runs some basic plotting of the metadata related to the patients
"""


import matplotlib.pyplot as plt
import seaborn as sns
import functions as fxn


__version__ = '0.1.1'
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
    plot_dir = fxn.plot_dir('analyse-metadata')

    # Get metadata
    meta = fxn.get_metadata()

    fig = plt.figure(figsize=(5, 5))
    sns.countplot(data=meta, hue='Sex', x='xRT-type')
    plt.ylabel('Count')
    plt.xlabel('Radiation type')
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'sexes.svg', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(5, 5))
    sns.violinplot(data=meta, x='xRT-type', y='Age-at-nadir', cut=0, inner='stick', hue='xRT-type')
    fxn.plot_signif_lines(0, 1, 88, fxn.mwu(meta[meta['xRT-type'] == 'Protons']['Age-at-nadir'],
                                            meta[meta['xRT-type'] == 'Photons']['Age-at-nadir']), 'black')
    plt.ylabel('Patient age (at nadir)')
    plt.xlabel('Radiation type')
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'ages.svg', dpi=300, bbox_inches='tight')
    plt.close()

    # Plotting which groups of patients had what prior chemotherapy
    fig = plt.figure(figsize=(5, 5))
    sns.countplot(hue='Prior-chemo', data=meta, x='xRT-type')
    plt.ylabel('Count')
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'prior-chemo.svg', dpi=300, bbox_inches='tight')
    plt.close()

    meta['Prior-chemo-regimen'] = meta['Prior-chemo-regimen'].replace('', 'None')

    fig = plt.figure(figsize=(12, 5))
    ax = sns.countplot(x='Prior-chemo-regimen', data=meta, hue='xRT-type',
                       order=['FOLFOX', 'FOLFIRINOX', 'FOLFIRINOX+Losartan', 'Carboplatin+Taxol',
                              'Gemcitabine+Cisplatin', 'None'])
    plt.ylabel('Count')
    plt.xticks(rotation=23)
    ax.set_xticklabels(['FOLFOX', 'FOLFIRINOX', 'FOLFIRINOX\n+Losartan', 'Carboplatin\n+Taxol',
                        'Gemcitabine\n+Cisplatin', 'None'])
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'prior-chemo-regimen.svg', dpi=300, bbox_inches='tight')
    plt.close()

    chemo_gp_order = ['5FU', 'Platin', 'None']
    prop_chemo_gp = []
    for gp in chemo_gp_order:
        for xrt in ['Photons', 'Protons']:
            # For each treatment type, record the percentage of each radiotherapy group that received it
            prop_chemo_gp.append([xrt, gp, (len(meta[(meta['xRT-type'] == xrt) & (meta['Prior-chemo-group'] == gp)]) /
                                            len(meta[meta['xRT-type'] == xrt])) * 100])

    prop_chemo_gp = fxn.list_to_df(prop_chemo_gp, ['xRT type', 'Prior chemo type', 'Percentage'], False)

    fig = plt.figure(figsize=(5, 5))
    sns.countplot(x='Prior-chemo-group', data=meta, hue='xRT-type', order=chemo_gp_order)
    plt.ylabel('Count')
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'prior-chemo-group.svg', dpi=300, bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(5, 5))
    sns.barplot(x='Prior chemo type', y='Percentage', data=prop_chemo_gp, hue='xRT type', order=chemo_gp_order)
    sns.despine(right=True, top=True)
    plt.savefig(plot_dir + 'prior-chemo-percentage.svg', dpi=300, bbox_inches='tight')
    plt.close()

    # Number of patients per treatment group
    save_name = plot_dir + 'treatment-counts.svg'
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    sns.countplot(x='diagnosis', hue='xRT-type', data=meta, legend=False)
    plt.xlabel('')
    plt.ylabel('Number of patients')
    plt.xticks(rotation=20)
    sns.despine(right=True, top=True)
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close()

    # Correlation between treatment length and time between baseline and nadir samples
    save_name = plot_dir + 'treatment-duration-V-nadir-delay.svg'
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111)
    # sns.scatterplot(x='TP1-2-d', y='Rx-d', hue='xRT-type', style='diagnosis', data=meta, s=100, legend='brief')
    g = sns.relplot(x='TP1-2-d', y='Rx-d', style='xRT-type', hue='diagnosis', data=meta, s=70)
    g.set(xlim=(0, None))
    # plt.legend(loc='right')
    plt.xticks(rotation=20)
    sns.despine(right=True, top=True)
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close()

    # Correlation between treatment length and time between nadir and recovery samples
    save_name = plot_dir + 'treatment-duration-V-recovery-delay.svg'
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111)
    sns.relplot(x='TP2-3-d', y='Rx-d', style='xRT-type', hue='diagnosis', data=meta, s=70)
    # plt.legend(loc='right')
    plt.xticks(rotation=20)
    sns.despine(right=True, top=True)
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close()

    save_name = plot_dir + 'treatment-duration-swarm.svg'
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    sns.swarmplot(x='diagnosis', y='Rx-d', hue='xRT-type', data=meta)
    plt.xlabel('')
    plt.ylabel('Treatment duration (days)')
    plt.xticks(rotation=20)
    sns.despine(right=True, top=True)
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.close()
