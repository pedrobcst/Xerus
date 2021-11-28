# Copyright (c) 2021 Pedro B. Castro
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
from pathlib import Path
import os, sys
project_path = str(Path(os.path.dirname(os.path.realpath(__file__))).parent) + os.sep # so convoluted..
if project_path not in sys.path:
    sys.path.append(project_path)
from peakutils import baseline
import seaborn as sns
from Xerus.settings.mplibcfg import Mplibcfg
from matplotlib import pyplot as plt
import pandas as pd


def remove_baseline(data: pd.DataFrame, poly_degree: int = 8, plot: bool = False) -> pd.DataFrame:
    """
    Remove baseline (background) from XRD pattern using polynomial interpolation (PeakUtils)

    Parameters
    ----------
    data : pd.DataFrame
        A pandas DataFrame containing the experimental data
    poly_degree : int, default: 8
        Polynomial degree.

    Returns
    -------
    pd.DataFrame
        Returns a dataframe with new intensity where new_intensity = old_intensity - baseline
        Also plot the data, baseline and difference

    """
    ## style configurations
    sns.set()
    sns.set_style("darkgrid", {"patch.edgecolor": "black", "axes.edgecolor": "black"})
    sns.set_context("talk")
    c = Mplibcfg()
    c.largs['fontweight'] = 'normal'
    c.pargs['markersize'] = 6
    c.largs['fontsize'] = 24
    c.lgdnargs['fontsize'] = 16

    ## Get baseline
    _base = baseline(data.int, deg=poly_degree)
    y_diff = data.int - _base

    if plot:
        ## Plot
        fig, ax = plt.subplots(figsize=(12,8))

        # Make plots.
        ax.plot(data.theta, data.int, label='Exp. Data', **c.pargs)
        ax.plot(data.theta, _base, 'red', lw=3, label='Background', **c.pargs)
        ax.plot(data.theta, y_diff, 'black', lw=2, label='Exp data - Background', **c.pargs)

        # Set styles
        ax.tick_params(**c.targs)
        ax.set_ylabel("Intensity (a.u)", **c.largs)
        ax.set_xlabel(r"2$\theta$ (deg.)", **c.largs)

        # Legend
        plt.legend(bbox_to_anchor=(1,1), **c.lgdnargs)

    # Modify data
    data.iloc[:, 1] = y_diff

    return data

def standarize_intensity(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standarize a XRD-pattern intensity by diving by Imax
        New intensity = Intensity / Intensity.max()
    Parameters
    ----------
    df : pd.DataFrame,
        A pandas dataframe containing the experimental data

    Returns
    -------
    pd.DataFrame
        Returns a standarized pandas dataframe.
    """

    df['int'] = df['int'] / df['int'].max()
    return df