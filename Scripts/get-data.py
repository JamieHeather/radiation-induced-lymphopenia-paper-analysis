# -*- coding: utf-8 -*-

"""
get-data.py

Takes raw project data downloaded from Adaptive (in v2 format) and converts it into an AIRR-seq compliant format,
  using immunoseq2airr.py (version 0.5.1 or greater).

"""

import gzip
import os
import functions as fxn
import subprocess

__email__ = 'jheather@mgh.harvard.edu'
__version__ = '0.1.1'
__author__ = 'Jamie Heather'


if __name__ == "__main__":
    out_conv_dir = fxn.make_check_dir(fxn.conv_data_dir)
    all_in_files = [x for x in os.listdir(fxn.raw_data_dir) if '.tsv' in x]
    all_in_files.sort()

    fxn.make_check_dir(fxn.conv_data_dir)

    for f in all_in_files:
        print("Processing raw data file " + f + " ...")

        base = f.split('.')[0]
        cmd = "python3 immunoseq2airr.py -i " + fxn.raw_data_dir + f + " -o " + out_conv_dir + base + " -z -or"
        subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)

        # If the file isn't already gzipped, make it so
        if f.endswith('.tsv'):
            with open(fxn.raw_data_dir + f, 'r') as in_file, gzip.open(
                    fxn.raw_data_dir + f + '.gz', 'wt') as zipped:
                zipped.writelines(in_file)

            os.unlink(fxn.raw_data_dir + f)
