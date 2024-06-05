# -*- coding: utf-8 -*-

"""
run-analysis.py
Run each of the relevant scripts in turn, in the right order, to produce the data and plots for analysis, in order.
"""

variables_to_keep = ['script', 'variables_to_keep', 'var_to_check'] + dir()
import os

scripts_order = ['illustrate-nadir.py',
                 'lymphopenia-incidence.py',
                 'analyse-metadata.py',
                 'analyse-tcrs.py',
                 'iterative-subsample-tcrs.py',
                 'analyse-flow.py',
                 'analyse-survival.py',
                 'cluster-specificities.py',
                 'functions.py',
                 'immunoseq2airr.py'
                 ]

# First check all scripts are there
any_missing = False
for script in scripts_order:
    if script not in os.listdir(os.getcwd()):
        print("Cannot find script \'" + script + "\' in this directory.")
        any_missing = True

if any_missing:
    raise IOError("Cannot proceed due to missing required scripts.")
else:
    print("Found all necessary scripts, proceeding to analysis.")
    import functions as fxn

# Ensure all the necessary folders are present
base_plot_dir = fxn.make_check_dir(fxn.base_plot_dir)
raw_data_dir = fxn.make_check_dir(fxn.raw_data_dir)
conv_data_dir = fxn.make_check_dir(fxn.conv_data_dir)

# Check all required input files are present in raw data directory
required_csv = ['chemoRT-lymphopenia.csv', 'flow-data.csv', 'metadata.csv']
any_missing = False
for csv in required_csv:
    if csv not in os.listdir(raw_data_dir):
        print("Cannot find data file \'" + csv + "\' in the Raw-Data directory.")
        any_missing = True

num_tsvs = len([x for x in os.listdir(raw_data_dir) if 'tsv' in x])
if num_tsvs != 60:
    print("Incorrect number of repertoire TSV files in Raw-Data directory. Should be 60, but there\'s " + str(num_tsvs))
    any_missing = True

if any_missing:
    raise IOError("Cannot proceed due to missing required input files.")
else:
    print("Found right number of input files, proceeding to analysis.")

# Then run all those scripts (except for the functions script and the immunoseq2airr script,
# which are both used by other scripts and not directly
visual_spacer = '-' * 50
for script in scripts_order[:-2]:
    print('\n'.join([visual_spacer, script, visual_spacer, '\n']))
    exec(compile(open(script, "rb").read(), script, 'exec'))
    # In between scripts remove any unused variables, to prevent clashes downstream
    for x in dir():
        if x not in variables_to_keep:
            del globals()[x]
    for x in dir():
        if x not in variables_to_keep:
            del locals()[x]

    visual_spacer = '-' * 50
