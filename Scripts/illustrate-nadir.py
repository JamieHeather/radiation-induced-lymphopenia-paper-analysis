#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
illustrate-nadir.py

Mock up a general schematic of ALC over the course of the three major timepoint samples

"""

import functions as fxn
import matplotlib.pyplot as plt
import seaborn as sns
import random

__version__ = '0.2.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def jitter(numb, val):
    """ Returns a number with random jitter within a given value, for strip chart plotting """
    return numb + random.uniform(-val, val)


# Sort plotting parameters
plt.rcParams.update({'font.size': 20, 'font.sans-serif': 'Arial'})

dat = []
iterations = 250
default_start = 1.5
default_end = 1.35
day_interval = 1
nadir = 0.3

# First have 10 days stable ALC
pre_range = .2
for d in range(-110, -100, day_interval):
    dat.append([d+1, default_start-pre_range])
    dat.append([d+1, default_start+pre_range])
# Then decrease steadily to a nadir
difference = default_start - nadir
daily_drop = difference / 100
current_alc = default_start
drop_range = 0.3
daily_range_gain = (drop_range - pre_range) / 100
current_range = pre_range
for d in range(-100, 0, day_interval):
    current_alc -= daily_drop
    dat.append([d+1, current_alc-current_range])
    dat.append([d+1, current_alc+current_range])
    current_range += daily_range_gain
# Then increase to non-lymphopenic levels
difference = default_end - nadir
daily_gain = difference / 100
gain_range = 0.4
daily_range_gain = (gain_range - current_range) / 100
for d in range(0, 100, day_interval):
    current_alc += daily_gain
    dat.append([d+1,current_alc-current_range])
    dat.append([d+1,current_alc+current_range])
    current_range += daily_range_gain
# Then have another 10 days of stability
post_range = 0.5
for d in range(100, 110, day_interval):
    dat.append([d+1, default_end-current_range])
    dat.append([d+1, default_end+current_range])


dat = fxn.list_to_df(dat, ['Day', 'ALC'], False)
dat['col'] = [' '] * len(dat)

plot_dir = fxn.plot_dir('illustrate-nadir')

fig = plt.figure(figsize=(6, 2))
g = sns.lineplot(x="Day", y="ALC", data=dat, palette=['red'], linewidth=2.5, hue='col', legend=False)
plt.ylim(0, 2)
plt.xlabel('')
plt.setp(g.set_xticks([]))
plt.setp(g.set_xticklabels([]))
sns.despine()
plt.savefig(plot_dir + 'schematic.svg', dpi=300, bbox_inches='tight')
plt.close()
