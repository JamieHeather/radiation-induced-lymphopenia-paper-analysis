#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
cluster-specificities.py

Perform rudimentary antigen prediction by clustering TCRs with published known specificity TCRs from VDJdb

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import Levenshtein as lev
import collections as coll
import functions as fxn
from graph_tool.all import *
import itertools as it
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
import math


__version__ = '0.2.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def rnd_col():
    return [x / 256 for x in list(np.random.choice(range(256), size=3))] + [1]


def split_genes(poss_genes):
    """
    :param poss_genes: a str containing either one or more genes, split by commas
    :return: a list of those genes split by commas (or a list of one if only one)
    """
    sep_genes = poss_genes.split(',')
    if isinstance(sep_genes, str):
        sep_genes = [sep_genes]
    return sep_genes


def tidy_genes(gene_str):
    """
    :param gene_str: name of a TCR gene
    :return: gene name without allele, sorted if ambiguous
    """
    sep_genes = split_genes(gene_str)
    gene_bits = list(set([x.split('*')[0] for x in sep_genes]))
    gene_bits.sort()
    return ','.join(gene_bits)


def collapse_hla(hla_str):
    """
    :param hla_str: HLA field from VDJdb
    :return: two digit resolution gene name(s), sorted if ambiguous
    """
    sep_genes = split_genes(hla_str)
    gene_bits = [x.split(':')[0] for x in sep_genes]
    gene_bits = list(set(gene_bits))
    gene_bits.sort()
    return ','.join(gene_bits)


def tcr_match(v_j_cdr3_1, v_j_cdr3_2, edit_distance):
    """
    :param v_j_cdr3_1: a TCR str in the format "TRBVxx|TRBJxx|C----CDR3----F"
    :param v_j_cdr3_2: another TCR str in the format "TRBVxx|TRBJxx|C----CDR3----F"
    :param edit_distance: int describing allowed CDR3 edit distance
    :return: true/false, for whether or not these two TCRs are a match within the defined criteria
    """
    # Easiest case: if they're the same, they're the same
    if v_j_cdr3_1 == v_j_cdr3_2:
        return True
    # Otherwise we need to break out the parts
    v1, j1, cdr3_1 = v_j_cdr3_1.split('|')
    v2, j2, cdr3_2 = v_j_cdr3_2.split('|')
    # Check for V/J gene matches (allowing one of ambiguous calls to match
    gene_matches = coll.defaultdict()
    gene_matches['v'] = [x for x in v1.split(',') if x in v2.split(',')]
    gene_matches['j'] = [x for x in j1.split(',') if x in j2.split(',')]
    # Found a V and a J match
    if gene_matches['v'] and gene_matches['j']:
        if cdr3_1 == cdr3_2:
            return True
        else:
            dist = lev.distance(cdr3_1, cdr3_2)
            if dist <= edit_distance:
                return True
            else:
                return False
    else:
        return False


def get_rnd_col(col_type):
    """
    :return: A random RGB colour (i.e. a list of 3x decimals in 255-space
    """
    if col_type == 'rgb':
        return [x/256 for x in list(np.random.choice(range(256), size=3))]
    elif col_type == 'bw':
        return [x/256 for x in list(np.random.choice(range(150, 220), size=1))] * 3
    else:
        raise ValueError("Unknown colour type: " + col_type)


def jamkey_3tp(clone_df, col_col, out_path):
    """
    Generate 'Jamkey' style clonal flow plots (basically a Sankey plot made by me, Jamie)
    :param clone_df: A donor-specific dataframe of clones (rows) with 3x float columns (TP1-3), & a categorical column
    :param col_col: Str of the categorical column in the df that dictates which clones should be highlighted (non-NaN)
    :param out_path: Str of the path to save the plot in
    :return: Nothing, just plot

    """
    hw = 0.2            # Half width value of the rectangles (clone sizes) to be plotting
    pad = 0.01          # Padding to insert between rectangles and polygons (clonal transitions)
    fig, ax = plt.subplots(figsize=(3.5, 4))
    start = {1: 0, 2: 0, 3: 0}
    end = {1: 0, 2: 0, 3: 0}
    zero = min(clone_df[['TP1', 'TP2', 'TP3']].min())

    for row in range(len(clone_df)):
        category = clone_df.iloc[row][col_col]
        if isinstance(category, str):
            rnd_col, alpha = get_rnd_col('bw'), 0.4
        else:
            raise IOError('Non string epitope category!')
        for tp in [1, 2, 3]:
            freq = clone_df.iloc[row]['TP' + str(tp)]
            if freq > zero:
                ax.add_patch(Rectangle((tp-hw, start[tp]), hw * 2, freq, facecolor=rnd_col, edgecolor='black', lw=0.5))
                end[tp] += freq

        # And add the connecting polygons
        for change in [(1, 2), (2, 3)]:
            x1, x2 = [change[0] + hw + pad] * 2
            x3, x4 = [change[1] - hw - pad] * 2
            y1 = start[change[0]]
            y2 = end[change[0]]
            y3 = end[change[1]]
            y4 = start[change[1]]
            # Account for clones absent in a particular sample
            freq1 = clone_df.iloc[row]['TP' + str(change[0])]
            freq2 = clone_df.iloc[row]['TP' + str(change[1])]
            if freq1 == zero:
                x1, x2 = [change[0]+0.5] * 2
                # y1, y2 = [y3 - ((y3 - y4) / 2)] * 2
                if y3 >= y1:
                    y1, y2 = [y3 - ((y3 - y1) / 2)] * 2
                else:
                    y1, y2 = [y1 - ((y1 - y3) / 2)] * 2
            if freq2 == zero:
                x3, x4 = [change[1]-0.5] * 2
                if y3 >= y1:
                    y3, y4 = [y3 - ((y3 - y1) / 2)] * 2
                else:
                    y3, y4 = [y1 - ((y1 - y3) / 2)] * 2
            if len(list(set([y1, y2, y3, y4]))) > 1:
                poly = [(x1, y1), (x2, y2), (x3, y3), (x4, y4),]
                ax.add_patch(Polygon(poly, facecolor=rnd_col, alpha=0.8))

        for tp in [1, 2, 3]:
            start[tp] = end[tp]

    ymax = (math.ceil(max(clone_df[['TP1', 'TP2', 'TP3']].sum())*1e6) / 1e6) * 1.02
    plt.ylim(0, ymax)
    plt.xlim(0.7, 3.3)
    plt.xticks([1, 2, 3])
    ax.set_xticklabels(['TP1', 'TP2', 'TP3'])
    # plt.ylabel('Frequency')
    sns.despine(right=True, top=True)
    plt.subplots_adjust(left=0.25, right=0.98, top=0.97, bottom=0.09)
    plt.savefig(out_path, dpi=300)  #, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    # go through each donor, and pick out the top X clones from each time point
    # then calculate the frequency of each of those clones in each time point in that donor

    ########################## CHANGABLE PARAMETERS ##############################################################
    top_x = 100                 # number of top-ranked TCRs to include per time point per donor
    edit_dist = 1               # levenshtein distance within which to cluster CDR3s (with matching V/Js)
    vdjdb_scores = ['2', '3']   # list of VDJdb scores to include (
    epip_threshold = 0.9        # how uniform clusters must be to be considered 'specific'
    ##############################################################################################################

    plt.rcParams.update({'font.size': 20, 'font.sans-serif': 'Arial'})
    plot_dir = fxn.plot_dir('trajectoriesSO-xRT' + str(top_x) + '-VDJdb' + ''.join(vdjdb_scores) +
                            '-edit' + str(edit_dist) + '-Ag' + str(int(epip_threshold * 100)))

    raw_tcrs = fxn.open_previous_data(fxn.conv_data_dir, 'raw-tcr-dict', 'pkl')

    dat = []
    tps = ['TP1', 'TP2', 'TP3']

    for donor in raw_tcrs:
        # Pick the top X TCRbs per time point
        tops = coll.defaultdict(list)
        for tp in tps:
            tops[tp] = raw_tcrs[donor][tp].most_common(top_x)
        all_tcrs = list(set([x[0] for x in tops['TP1']] +
                            [y[0] for y in tops['TP2']] +
                            [z[0] for z in tops['TP3']]))

        # Then for each TCR figure out its frequency in those respective time points
        counts = dict(zip(tps, [sum(raw_tcrs[donor][x].values()) for x in tps]))
        for tcr in all_tcrs:
            freqs = [raw_tcrs[donor][tp][tcr]/counts[tp] for tp in tps]
            normed_freqs = [x/max(freqs) for x in freqs]
            dat.append([donor, tcr] + freqs)

    dat = fxn.list_to_df(dat, ['Donor', 'TCR', 'TP1', 'TP2', 'TP3'], False)
    dat = dat.sample(frac=1)    # Shuffle the order

    donors = list(set(dat['Donor']))
    donors.sort()

    dat['Offset'] = np.nan
    top_freq_dat = ''
    dat_building_dir = fxn.make_check_dir(plot_dir + 'top-clones-' + str(top_x) + '/')

    for d in donors:
        d_dat = dat.loc[dat['Donor'] == d][['TP1', 'TP2', 'TP3']]
        # Find the lowest non-zero value used in this donor
        min_d_val = d_dat.mask(d_dat == 0).min().min()
        # Offset every value for this donor
        dat.loc[dat['Donor'] == d, 'Offset'] = min_d_val / 2
        # Plot the proportion of each time point accounted for in this dataframe
        sum_df = pd.DataFrame(d_dat.sum(axis=0))
        sum_df.reset_index(inplace=True)
        sum_df.columns = ['Timepoint', 'Sum Freq']

        plt.figure(figsize=(3, 3))
        sns.barplot(data=sum_df, x='Timepoint', y='Sum Freq')
        dat_building_dir = fxn.make_check_dir(plot_dir + 'top-clones-' + str(top_x) + '/')
        plt.ylim(0, 1)
        plt.savefig(dat_building_dir + d + '.svg', bbox_inches='tight')
        plt.close()

        # Save these sums to a running tally
        sum_df['Donor'] = d
        if len(top_freq_dat) == 0:
            top_freq_dat = sum_df
        else:
            top_freq_dat = pd.concat([top_freq_dat, sum_df])

    top_freq_dat.reset_index(inplace=True)

    sns.catplot(data=top_freq_dat, x='Timepoint', y='Sum Freq', col='Donor', kind='bar')

    for tp in tps:
        dat[tp] = dat[tp] + dat['Offset']

    dat['TP1-TP2'] = np.log2(dat['TP2'] / dat['TP1'])
    dat['TP2-TP3'] = np.log2(dat['TP3'] / dat['TP2'])

    # Plot total cohort overview
    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=dat, x='TP1-TP2', y='TP2-TP3', hue='Donor', alpha=.2, s=20)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.xlabel('log2(TP2/TP1)')
    plt.ylabel('log2(TP3/TP2)')
    plt.ylim(-10, 13)
    plt.xlim(-10, 13)
    plt.savefig(plot_dir + 'initial-plot.svg', bbox_inches='tight')
    plt.close()

    # read in vdjdb data (placing the 'vdjdb.slim.txt' file at Raw-Data/VDJdb/)
    # Write out a file of all TCRs, plus separate alpha/beta chains
    vdjdb_all = []
    tcr_counter = coll.defaultdict(fxn.nest_list)
    with fxn.opener(fxn.raw_data_dir + 'VDJdb/vdjdb.slim.txt', 'r') as in_file:
        line_count = 0
        for line in in_file:
            bits = line.rstrip().split('\t')
            chain = bits[0]
            if line_count == 0:
                headers = bits

            # Select human TCRbs with V/J/CDR3 information
            elif bits[2] == 'HomoSapiens' and chain == 'TRB' and bits[7] and bits[8] and bits[1]:
                # # Again, get rid of ambiguous gene calls
                # if ',' in bits[7] or ',' in bits[8]:
                #     continue

                # Screen out dubious CDR3s
                if len(bits[1]) < 8 or bits[1][-1] not in ['F', 'W']:
                    continue

                # Then record the tcr information
                v = tidy_genes(bits[7])
                j = tidy_genes(bits[8])
                cdr3 = bits[1]
                tcr = '|'.join([v, j, cdr3])

                # And the epitope information
                pep = bits[3]
                hla = collapse_hla(bits[11])
                epip = '|'.join([hla, pep])
                source = '|'.join(bits[4:6])
                score = bits[15]
                out_row = [tcr, epip, v, j, cdr3, hla, pep, source, score]
                vdjdb_all.append(out_row)
                tcr_counter[tcr][epip].append(out_row)

            line_count += 1

    vdjdb_all = fxn.list_to_df(vdjdb_all,
                               ['TCR', 'Epitope', 'V', 'J', 'CDR3', 'HLA', 'Peptide', 'Source', 'Score'], False)

    # Go through and eliminate cross-reactive (or miss-annotated) TCRs
    vdjdb = []
    for tcr in tcr_counter:
        # Skip if >1 different epitope  # TODO make more clever, e.g. look for closely related epitopes
        if len(tcr_counter[tcr]) > 1:
            continue
        else:
            epip = list(tcr_counter[tcr].keys())[0]
            if len(tcr_counter[tcr][epip]) == 1:
                vdjdb.append(tcr_counter[tcr][epip][0])
            # Otherwise pick the highest scoring entry
            else:
                max_score = 0
                for entry in tcr_counter[tcr][epip]:
                    if int(entry[-1]) > max_score:
                        chosen = entry
                        max_score = int(entry[-1])
                if max_score == 0:
                    chosen = entry
                vdjdb.append(chosen)

    vdjdb = fxn.list_to_df(vdjdb, ['TCR', 'Epitope', 'V', 'J', 'CDR3', 'HLA', 'Peptide', 'Source', 'Score'], False)

    # Go through the CDRs in both the donors and in VDJdb and cluster together
    vdjdb_filt = vdjdb.loc[vdjdb['Score'].isin(vdjdb_scores)]
    print(len(vdjdb_filt))
    # vdjdb_filt = vdjdb.loc[vdjdb['Score'].isin(['2', '3'])]

    # Remodel the patient data to match the VDJdb data, and stick together to allow cross-comparison
    dat_match = dat[['TCR']]
    dat_match['Epitope'] = dat['Donor']
    dat_match[['V', 'J', 'CDR3']] = dat_match['TCR'].str.split('|', expand=True)
    dat_match['HLA'] = '?'
    dat_match['Peptide'] = '?'
    dat_match['Source'] = dat['Donor']
    dat_match['Score'] = '0'

    vdjdb_filt = pd.concat([dat_match, vdjdb_filt], ignore_index=True)
    print(len(vdjdb_filt))

    vdjdb_combs = it.combinations(vdjdb_filt.index, 2)
    edge_freqs = coll.defaultdict(list)
    vertex_freqs = coll.Counter()

    # Then generate a network, in which matching TCRs (same V/J, CDR3 <= X edits away) are joined together by an edge
    g = Graph(directed=False)
    v_dict = coll.defaultdict()
    epip_dict = coll.defaultdict()
    v_tcr = g.new_vertex_property('string')
    v_ag = g.new_vertex_property('string')
    v_s = g.new_vertex_property('string')
    v_shape = g.new_vertex_property('string')
    v_col = g.new_vertex_property('vector<double>')

    for tcrx1, tcrx2 in vdjdb_combs:
        tcr1 = vdjdb_filt.loc[tcrx1]
        tcr2 = vdjdb_filt.loc[tcrx2]
        match = tcr_match(tcr1['TCR'], tcr2['TCR'], edit_dist)
        if match:
            print(tcrx1, tcrx2)

            if tcr1['TCR'] not in v_dict:
                v_dict[tcr1['TCR']] = g.add_vertex()
                v_tcr[v_dict[tcr1['TCR']]] = tcr1['TCR']
                v_s[v_dict[tcr1['TCR']]] = tcr1['Source']
                v_ag[v_dict[tcr1['TCR']]] = tcr1['Epitope']
            if tcr2['TCR'] not in v_dict:
                v_dict[tcr2['TCR']] = g.add_vertex()
                v_tcr[v_dict[tcr2['TCR']]] = tcr2['TCR']
                v_s[v_dict[tcr2['TCR']]] = tcr2['Source']
                v_ag[v_dict[tcr2['TCR']]] = tcr2['Epitope']

            # Vertex source (antigen or donor)
            if tcr1['Epitope'] not in epip_dict:
                epip_dict[tcr1['Epitope']] = rnd_col()
            if tcr2['Epitope'] not in epip_dict:
                epip_dict[tcr2['Epitope']] = rnd_col()
            # Vertex shape
            if tcr1['Epitope'].startswith('HLA'):
                v_shape[v_dict[tcr1['TCR']]] = 'circle'
            else:
                v_shape[v_dict[tcr1['TCR']]] = 'triangle'
            if tcr2['Epitope'].startswith('HLA'):
                v_shape[v_dict[tcr2['TCR']]] = 'circle'
            else:
                v_shape[v_dict[tcr2['TCR']]] = 'triangle'
            #
            v_col[v_dict[tcr1['TCR']]] = epip_dict[tcr1['Epitope']]
            v_col[v_dict[tcr2['TCR']]] = epip_dict[tcr2['Epitope']]
            g.add_edge(v_dict[tcr1['TCR']], v_dict[tcr2['TCR']])

    # Remove unwanted edges
    remove_parallel_edges(g)
    remove_self_loops(g)

    print("\t\t\t\t\t\t\t\t\t\tPlotting...")
    sz_plt = (8000, 8000)
    sz_vertex = 5

    # top_x = 500
    # edit_dist = 1
    # vdjdb_scores = ['0', '1', '2', '3']
    save_prefix = '-'.join(['dat' + str(top_x), 'score' + ''.join(vdjdb_scores), 'edit' + str(edit_dist), ''])

    graph_draw(g, output=plot_dir + save_prefix + 'tcr-labelled-noParallel.pdf',
               vertex_text=v_tcr, vertex_size=sz_vertex, vertex_fill_color=v_col,
               geometry=sz_plt, output_size=sz_plt, vertex_font_family='sans')
    graph_draw(g, output=plot_dir + save_prefix + 'tag-labelled-noParallel.pdf',
               vertex_text=v_ag, vertex_size=sz_vertex,
               vertex_fill_color=v_col, geometry=sz_plt, output_size=sz_plt)
    graph_draw(g, output=plot_dir + save_prefix + 'noParallel.pdf', vertex_size=sz_vertex,
               vertex_fill_color=v_col, geometry=sz_plt, output_size=sz_plt)

    # Try to figure out the clusters that contain donor-derived sequences
    processed_vertices = []          # Nodes (TCRs) that have already been processed
    found_tcrs = coll.defaultdict()  # xRT TCRs with their potential assigned specificities based on clustering

    # Loop through all TCRs (vertices)
    for v in v_dict:
        # Skip it it's already been processed
        if v in processed_vertices:
            continue
        # Otherwise check to see if it's an xRT patient-derived TCRs
        if v_ag[v_dict[v]][0] not in ['T', 'R']:
            continue

        # Get all the TCRs in the cluster, iterating over each vertex until no new vertices are added
        clustered = coll.Counter([v])
        last_size = 1
        finished = False
        # TODO could marginally increase speed here by skipping over nodes that have already had their neighbours added
        while not finished:
            current_nodes = clustered.keys()
            tmp_clustered = coll.Counter()
            for node in current_nodes:
                tmp_clustered += coll.Counter([v_tcr[x] for x in g.iter_out_neighbours(v_dict[node])])
            clustered += tmp_clustered
            if last_size == len(clustered):
                finished = True
            else:
                last_size = len(clustered)

        # Determine whether or not this is an eligible cluster (i.e. any VDJdb-derived TCRs within)
        epip_associations = [v_ag[v_dict[x]] for x in clustered if v_ag[v_dict[x]].startswith('HLA')]
        if not epip_associations:
            continue

        # If so, determine whether there's an obvious consensus specificity
        epip_count = coll.Counter(epip_associations)
        if epip_count.most_common(1)[0][1] / sum(epip_count.values()) < epip_threshold:
            continue

        # TODO could maybe add a threshold here, so that there needs to be X many VDJdb derived associations per cluster
        putative_ag = epip_count.most_common(1)[0][0]
        # If so, apply that specificity to all xRT-derived clustered TCRs (ct) and store the call
        for ct in clustered:
            if v_ag[v_dict[ct]][0] in ['T', 'R']:
                found_tcrs[ct] = putative_ag

        # Finally add clustered to processed_vertices and continue to next cluster
        processed_vertices += clustered.keys()

    # Plot! First add the inferred specificities, antigens, and HLAs to the dataframe
    dat[['Epitope', 'HLA', 'Peptide', 'Source', 'Species', 'Antigen']] = np.nan
    for tcr in found_tcrs:
        for d in donors:
            if len(dat.loc[(dat['TCR'] == tcr) & (dat['Donor'] == d)]) == 1:
                source = coll.Counter(vdjdb_filt.loc[vdjdb_filt['Epitope']
                                                     == found_tcrs[tcr]]['Source']).most_common(1)[0][0]
                antigen, species = source.split('|')
                hla, peptide = found_tcrs[tcr].split('|')
                dat.loc[dat.loc[(dat['TCR'] == tcr) & (dat['Donor'] == d)].index, 'Epitope'] = found_tcrs[tcr]
                dat.loc[dat.loc[(dat['TCR'] == tcr) & (dat['Donor'] == d)].index, 'HLA'] = hla
                dat.loc[dat.loc[(dat['TCR'] == tcr) & (dat['Donor'] == d)].index, 'Peptide'] = peptide
                dat.loc[dat.loc[(dat['TCR'] == tcr) & (dat['Donor'] == d)].index, 'Source'] = source
                dat.loc[dat.loc[(dat['TCR'] == tcr) & (dat['Donor'] == d)].index, 'Species'] = species
                dat.loc[dat.loc[(dat['TCR'] == tcr) & (dat['Donor'] == d)].index, 'Antigen'] = antigen

    fig = plt.figure(figsize=(5, 5))
    sns.scatterplot(data=dat, x='TP1-TP2', y='TP2-TP3', alpha=.1, color='gray')
    sns.scatterplot(data=dat.loc[dat['TCR'].isin(found_tcrs)], x='TP1-TP2', y='TP2-TP3', alpha=.9, hue='Species')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.savefig(plot_dir + 'specificities-over-trajectories.svg', bbox_inches='tight')
    plt.close()

    # Plot donor-specific plots
    for d in donors:
        donor_dat = dat.loc[dat['Donor'] == d]
        for variable in ['Species', 'HLA', 'Epitope']:
            # Plot trajectories
            fig = plt.figure(figsize=(5, 5))
            sns.scatterplot(data=dat, x='TP1-TP2', y='TP2-TP3', alpha=.05, color='lightgray', s=5)
            sns.scatterplot(data=donor_dat, x='TP1-TP2', y='TP2-TP3', alpha=.2, color='black')
            sns.scatterplot(data=donor_dat.loc[dat['TCR'].isin(found_tcrs)], x='TP1-TP2', y='TP2-TP3',
                            alpha=.9, hue=variable)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.savefig(plot_dir + d + '-' + variable + '-specificities-over-trajectories.svg', bbox_inches='tight')
            plt.close()

        # And plot just HLA specific counts
        hla_order = [x[0] for x in coll.Counter(donor_dat.loc[dat['TCR'].isin(found_tcrs)]['HLA']).most_common()]
        pre_df = list(donor_dat.loc[dat['TCR'].isin(found_tcrs)]['HLA'])
        if pre_df:
            hla_df = fxn.list_to_df(pre_df, ['HLA'], False)
            fig = plt.figure(figsize=(5, 5))
            gg = sns.countplot(data=hla_df, x='HLA', order=hla_order)
            gg.set_xticklabels(labels=hla_order, rotation=80)
            plt.xlabel('')
            plt.savefig(plot_dir + 'hla-' + d + '-associated-counts.svg', bbox_inches='tight')
            plt.close()

    # And separate plot for each of the major antigen origin species
    for sp in [x[0] for x in coll.Counter(dat.loc[dat['TCR'].isin(found_tcrs)]['Species']).most_common(6)]:
        sns.scatterplot(data=dat, x='TP1-TP2', y='TP2-TP3', alpha=.1, color='gray')
        fig = plt.figure(figsize=(5, 5))
        sns.scatterplot(data=dat.loc[(dat['TCR'].isin(found_tcrs)) & (dat['Species'] == sp)],
                        x='TP1-TP2', y='TP2-TP3', alpha=.9, hue='Species')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(plot_dir + 'specificities-over-trajectories-' + sp + '.svg', bbox_inches='tight')
        plt.close()

    jamkey_dir = fxn.make_check_dir(plot_dir + 'jamkeys')

    # Fix the inappropriately labelled human antigen entries (where species and antigen are mixed)
    wrong_ag = list(set(dat.loc[dat['Antigen'] == 'HomoSapiens']['Species']))
    for wag in wrong_ag:
        old_nam = 'HomoSapiens|' + wag
        new_nam = wag + '|HomoSapiens'
        dat.loc[dat['Source'] == old_nam, 'Antigen'] = wag
        dat.loc[dat['Source'] == old_nam, 'Species'] = 'HomoSapiens'
        dat.loc[dat['Source'] == old_nam, 'Source'] = new_nam

    # Finally run the Samkey-style flow plotting code
    for d in donors:
        print(d)
        tmp = dat.loc[dat['Donor'] == d]
        tmp = tmp.sort_values(by=['TP1', 'TP2', 'TP3'])  #, ascending=False)
        jamkey_3tp(tmp, 'Species', jamkey_dir + 'clusters-' + d + '.svg')
