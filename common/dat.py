# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 11:46:23 2014

@author: alexey

data analysis tools (dat)
"""
import os
import re
import fcntl
import tables
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
from contextlib import contextmanager


# %%
from IPython.display import display, HTML

from pygments import highlight
from pygments.lexers import SqlLexer
from pygments.formatters import HtmlFormatter


def format_column(name):
    """Converts CamelCase and other weird column name formats to pothole case.
    """
    name = name.replace(' ', '_').replace('(', '_').replace(')', '')
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def print_sql(sql_query, style='default'):
    """Print SQL code with syntax highlighting.

    To view all available color styles, run::

        from pygments.styles import get_all_styles
        print(sorted(get_all_styles()))
    """
    formatter = HtmlFormatter(style=style)
    display(HTML('<style type="text/css">{}</style>{}'.format(
        formatter.get_style_defs('.highlight'),
        highlight(sql_query, SqlLexer(), formatter))))


# %%
def print2(a, b, *args, x=60):
    template = '{:%d}{}' % x
    formatted_template = template.format(a, b)
    for arg in args:
        formatted_template += ' ' + str(arg)
    print(formatted_template)


def print_full(x):
    """Print the entire Dataframe / Series
    """
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    print(x)
    pd.reset_option('display.max_columns')
    pd.reset_option('display.max_rows')


def remove_duplicates(seq, keep_null=True):
    if keep_null:
        seen = set()
    else:
        seen = set([np.nan, None, ''])
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


@contextmanager
def open_hdf5_exclusively(filename, mode='a'):
    try:
        fd = os.open(filename, os.O_CREAT | os.O_RDWR)
        fcntl.lockf(fd, fcntl.LOCK_EX)
        f = tables.openFile(filename, mode)
        yield f
    except:
        raise
    finally:
        f.close()


# %% PLOTTING
def make_venn_diagram(data_sets, set_labels, plot_title, output_filename=None):
    """Function for making a Venn diagram of pre-defined shape and font size
    """
    plt.figure(figsize=(12, 12))
    plt.title(plot_title, fontsize=20)

    if len(data_sets) == 2:
        v = venn2(data_sets, set_labels=set_labels)
        c = venn2_circles(data_sets, linestyle='dashed')
    elif len(data_sets) == 3:
        v = venn3(data_sets, set_labels=set_labels)
        c = venn3_circles(data_sets, linestyle='dashed')
    else:
        raise Exception('Wrong number of mutation sets: {}!'.format(len(data_sets)))

    # Resize set counts
    for label_id in ['10', '01', '11', '100', '010', '001', '110', '101', '011', '111']:
        v.get_label_by_id(label_id).set_fontsize(18)

    # Resize set labels
    for label_id in ['A', 'B', 'C']:
        v.get_label_by_id(label_id).set_fontsize(20)

    # Change circle line style and width
    c[0].set_ls('dotted')
    c[1].set_ls('dashed')
    c[2].set_lw(1.0)

    if output_filename is not None:
        # Typical output filename:
        # '/home/kimlab1/strokach/working/elaspic/reports/14-11-xx/venn_ddg_ddgh2o_dtm.png'
        plt.savefig(output_filename, dpi=150, bbox_inches='tight')


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input),
    to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y
