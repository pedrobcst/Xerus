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
import os
import sys
from pathlib import Path
from typing import List, Union

import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from matplotlib import pyplot as plt
from plotly.graph_objects import Figure
from Xerus.engine.gsas2utils import get_plot_data
from Xerus.settings.mplibcfg import Mplibcfg
from Xerus.utils.tools import formula_to_latex
from Xerus.utils.tools import plotly_add as to_add

# sys.path.append(GSAS2_BIN)
# import GSASIIscriptable as G2sc  # Ignore Warning


def plot_all_gpx(path: str, outfolder: str = "fig/", fmt: str = ".jpeg", engine: str = "matplotlib") -> None:
    """
    Auxiliary function to plot all gpx files in a folder.
    It will create a folder outside the the gpx folder files
    We are assuming that the folder structure is as follows:
    working_folder/gsas_files/*.gpx
    It will save at
    working_folder/outfolder/figure_name.extensions

    Parameters
    ----------
    path : str
        Path to GPX files to be plotted
    outfolder : str, default: "fig/"
        Path to folder to save the plots
    fmt : str, default: ".jpeg"
        File extension for saving the plots
    engine : str, default: "matplotlib"
        Defines plotting backend ("matplotlib"/"plotly")

    Returns
    -------
    None
    """

    if engine not in ["plotly", "matplotlib"]:
        raise NotImplementedError("Plotting backend not implemented.")

    gpx_files = [os.path.join(path, file) for file in os.listdir(path) if file.endswith(".gpx") and ".bak" not in file]
    # print(gpx_files)
    print("PATH:", path)

    folder = os.path.join(*os.path.dirname(path).split(os.sep)[:-1], outfolder)
    print("FOLDER:", folder)
    if not os.path.isdir(folder):
        os.mkdir(folder)
    for gpx in gpx_files:
        try:
            plot_gpx(gpx, outfolder=folder, fmt=fmt, engine=engine)
            plt.close()
        except KeyError:
            continue


def plot_gpx(path: str, offset: float = 0.2, outfolder: str = "curr", fmt: str = ".jpg", save: bool = True,
             engine: str = "matplotlib", type: str = "bragg", s_info: List[pd.DataFrame] = None, **kwargs) -> Union[Figure, None]:
    """
    Plot a .gpx file using matplotlib or plotly.


    Parameters
    ----------
    path : str
        Path to GSASII .gpx file.
    offset : float, default: 0.2
        Set the offset between the plots (Exp data / Simulated patterns etc)
    outfolder : str, default: "curr"
        Set the folder for exporting the generated plot(s). "curr" defaults to current folder
    fmt : str, default: ".jpg"
        Set the file extension to export the plot (mplib)
    save : bool, default: True
        If this is True, the figure will be saved to `outfolder`
    engine : str, default: "matplotlib"
        Set which backend to use ("matplotlib"/"plotly")
    type : str, default: "bragg"
        Define in which way to show the simulated patterns
            - "bragg": Will show the simulated patterns as bragg positions ticks "|". Supported by both matplotlib and plotly
            - "simul": Plot whole simulated patterns. Supported only for "plotly" backend now.
    s_info : List[pd.DataFrame], default: None
        In case that `type` = `simul`, provide a list of dataframes containing the simulated patterns.
        This is done automatically when using through the XRay API.

    Returns
    -------
    A figure or None.
    """
    if engine == "matplotlib":
        plot_gpx_mplib(path=path,
                       offset=offset,
                       outfolder=outfolder,
                       save=save,
                       fmt=fmt)
        return None
    elif engine == "plotly":
        if type == "bragg":
            fig = plot_gpx_plotly(path=path,
                                  offset=offset,
                                  outfolder=outfolder,
                                  save=save)
        else:
            fig = plot_gpx_plotly_simul(path=path, offset=offset,outfolder=outfolder,save=save,s_info=s_info, **kwargs)
        return fig
    else:
        raise NotImplementedError("Not implemented.")


def plot_gpx_mplib(path: str, offset: float = 0.2, outfolder: str = "curr", fmt: str = ".jpg",
                   save: bool = True) -> None:
    """
    Parse a GPX Project and plot the results using matplotlib.
    Can also export the image into a desired folder.

    Parameters
    ----------
    path : str
        Path to .gpx file to be plotted
    offset : float, default: 0.2
        Offset between experimental patterns and simulations
    outfolder : str, default: "curr"
        Folder to save the figure if `save` is True
    fmt : str, default: ".jpg"
        Image format to be saved.
    save : bool, default: True
        If `save` is True a image will be saved with `fmt` at `outfolder`.

    Returns
    -------
    None
    """
    import gc

    ## style configurations
    sns.set()
    sns.set_style("darkgrid", {"patch.edgecolor": "black", "axes.edgecolor": "black"})
    sns.set_context("talk")
    c = Mplibcfg()
    c.largs['fontweight'] = 'normal'
    c.pargs['markersize'] = 6
    c.largs['fontsize'] = 24
    c.lgdnargs['fontsize'] = 16

    ## Load up GSAS II related information
    bragg, x, yobs, ycalc, ydiff, rwp, wt_map = get_plot_data(path, offset)

    ## prepare figure
    fig, ax = plt.subplots(figsize=(12, 8))
    # plot data, ycalc y diff
    ax.plot(x, yobs, 'o-', color="black", fillstyle="none", label="Y$_{obs}$", markevery=3, **c.pargs)
    ax.plot(x, ycalc, "-", color="red", label="Y$_{calc}$", **c.pargs)
    ax.plot(x, ydiff, '-', color="blue", label="Y$_{obs}$ - Y$_{calc}$")

    # plot bragg positions
    for i, df in enumerate(bragg):
        bragg_positions = df.loc[:, "2-theta"]
        phase_name = df.phase_name.iat[0]  # get phase name
        wtfrac = wt_map[phase_name]  # get wt% for phase name
        formula = phase_name.split("-")[0]  # get the formula
        spg = phase_name.split("-")[1:]  # get the spg.
        print(phase_name)
        if len(spg) > 1:
            spg = "-".join(spg)
        else:
            try:
                spg = spg[0]
            except IndexError:
                spg = "None"  # for old files..
        formula = formula_to_latex(formula)
        phase_name = formula + "-" + spg
        y_dummy = np.zeros(len(bragg_positions)) - (i + 2) * offset
        ax.plot(bragg_positions, y_dummy, "|", ms=24, label=phase_name + " wt:{0:.2f}%".format(wtfrac * 100))

    # axis and labels
    ax.tick_params(**c.targs)
    ax.set_xlabel(r"2$\theta$ (deg.)", **c.largs)
    ax.set_ylabel("Intensity (arb. units)", **c.largs)
    rwp_format = "$Rwp$ = {0:.2f}%".format(rwp)
    l = ax.legend(**c.lgdnargs)
    l.set_title(rwp_format, prop={'size': c.lgdnargs['fontsize']})
    l.get_title().set_position((0, 0))
    l._legend_box.align = "left"
    l.set_bbox_to_anchor((1, 1))
    ax.set_yticks([])
    fig_name = os.path.basename(path).replace(".gpx", fmt)
    dpi = 300  #

    if outfolder == "curr":  # if current folder
        outfolder = os.getcwd()
    if save:
        fig.savefig(os.path.join(outfolder, fig_name), dpi=dpi, bbox_inches="tight")
        plt.close()
        fig.clf()
        del(fig)
        del(ax)
        gc.collect()


def plot_gpx_plotly(path: str, offset: float = 0.3, outfolder: str = "curr", save: bool = True) -> Union[Figure, None]:
    """
    Parse a GPX Project and plot the results using plotly.
    Can also export the image into a desired folder.

    Parameters
    ----------
    path : str
        Path to .gpx file to be plotted
    offset : float, default: 0.2
        Offset between experimental patterns and simulations
    outfolder : str, default: "curr"
        Folder to save the figure if `save` is True.
        "curr" will save to current folder.
    save : bool, default: True
        If `save` is True a image will be saved with HTML at `outfolder`.

    Returns
    -------
    plotly.graph_objects.Figure object
    """
    bragg, x, yobs, ycalc, ydiff, rwp, wt_map = get_plot_data(path, offset)

    # Merge data into dataframes for easy plotting
    df_yobs = pd.DataFrame(data=zip(x, yobs), columns=['theta', 'Y'])
    df_yobs['label'] = 'Yobs'
    df_ycalc = pd.DataFrame(data=zip(x, ycalc), columns=['theta', 'Y'])
    df_ycalc['label'] = 'Ycalc'
    df_ydiff = pd.DataFrame(data=zip(x, ydiff), columns=['theta', 'Y'])
    df_ydiff['label'] = 'Ydiff'

    df = pd.concat([df_yobs, df_ycalc, df_ydiff], axis=0)
    df['rwp'] = rwp

    # Concatanete extra phase vals.
    if len(wt_map) > 1:
        for key, val in wt_map.items():
            df[key] = val

    # Setup legend "title' and hover information
    title = "Rwp: {:.2f}%".format(rwp)
    extra = ['rwp'] + [key for key in wt_map.keys() if len(wt_map) > 1]

    # Create plot
    fig = px.scatter(data_frame=df,
                     x='theta',
                     y='Y',
                     color='label',
                     hover_data=extra,
                     hover_name="label",
                     labels={
                         "Y": "Intensity (a.u)",
                         "theta": r"$\text{2}\theta \text{ (deg.)}$",
                         "label": title
                     },
                     title=os.path.basename(path),
                     width=1000,
                     height=600)

    # Setup hover template for main plot
    t = fig.data[0].hovertemplate
    hover_t = "<b>%{hovertext}</b><br><br>2theta=%{x}<br>Rwp=%{customdata[0]}" + t.split("rwp=%{customdata[0]}")[-1]
    fig.data[0].hovertemplate = hover_t
    fig.data[1].hovertemplate = hover_t
    fig.data[2].hovertemplate = hover_t

    # Set up mode of ycalc and ydiff
    # Line
    fig.data[1].mode = "lines"
    fig.data[2].mode = "lines"

    # Line color
    fig.data[1].marker['color'] = "red"
    fig.data[2].marker['color'] = "blue"

    # Set up open markers for exp data
    fig.data[0].marker['color'] = "black"
    fig.data[0].marker['symbol'] = "circle-open"

    # plot bragg positions
    for i, df in enumerate(bragg):
        hover_base = "theta: %{x}<br>Phase: %{text}<br>"
        colors = sns.color_palette(n_colors=len(bragg))
        bragg_positions = df.loc[:, "2-theta"]
        size = df.loc[:, 'F_calc']
        size = 100 * size/size.max()
        phase_name = df.phase_name.iat[0]  # get phase name
        if len(wt_map) > 1:  # unnecessary if only one phase
            wtfrac = wt_map[phase_name]  # get wt% for phase name
            cdata = [wtfrac] * len(bragg_positions)
            ht = hover_base + "Frac. Pct:%{customdata}<extra></extra>"
        else:
            ht = hover_base + "<extra></extra>"  # used to get rigged of the extra name
            cdata = [[]]
        y_dummy = np.zeros(len(bragg_positions)) - (i + 2) * offset  # set up position of braggs to plot
        marker = {
            "symbol": "line-ns",
            "line_width": 1,
            "size": 12,
            "line": {"color": colors.as_hex()[i]}
        }
        fig.add_scatter(x=bragg_positions,
                        y=y_dummy,
                        mode="markers",
                        text=[phase_name] * len(bragg_positions),
                        name=phase_name, marker=marker,
                        customdata=cdata,
                        hovertemplate=ht)
        fig_name = os.path.basename(path).replace(".gpx", ".html")

        fig.update_layout(modebar_add=to_add)
        if outfolder == "curr":  # if current folder
            outfolder = os.getcwd()
        if save:
            fig.write_html(os.path.join(outfolder, fig_name), include_mathjax="cdn")

    return fig


def plot_gpx_plotly_simul(path: str, s_info: List[pd.DataFrame], offset: float = 0.3, outfolder: str = "curr", save: bool = False, width = 1000, height = 600) -> Union[Figure, None]:
    """
    Parse a GPX Project and plot the results using Plotly.
    This function just plots the full simulation of the patterns instead of using the bragg position ticks "|" as indication
    Can also export the image into a desired folder.

    Parameters
    ----------
    path : str
        Path to .gpx file to be plotted
    s_info: List[pd.DataFrame]
        A list of pandas dataframes containing the simulation information of each pattern
    offset : float, default: 0.2
        Offset between experimental patterns and simulations
    outfolder : str, default: "curr"
        Folder to save the figure if `save` is True
    save : bool, default: True
        If `save` is True a image will be saved with HTML at `outfolder`.

    Returns
    -------
    plotly.graph_objects.Figure object
    """
    bragg, x, yobs, ycalc, ydiff, rwp, wt_map = get_plot_data(path, offset)

    # Merge data into dataframes for easy plotting
    df_yobs = pd.DataFrame(data=zip(x, yobs), columns=['theta', 'Y'])
    df_yobs['label'] = 'Yobs'
    df_ycalc = pd.DataFrame(data=zip(x, ycalc), columns=['theta', 'Y'])
    df_ycalc['label'] = 'Ycalc'
    df_ydiff = pd.DataFrame(data=zip(x, ydiff), columns=['theta', 'Y'])
    df_ydiff['label'] = 'Ydiff'

    df = pd.concat([df_yobs, df_ycalc, df_ydiff], axis=0)
    df['rwp'] = rwp

    # Concatanete extra phase vals.
    if len(wt_map) > 1:
        for key, val in wt_map.items():
            df[key] = val

    # Setup legend "title' and hover information
    title = "Rwp: {:.2f}%".format(rwp)
    extra = ['rwp'] + [key for key in wt_map.keys() if len(wt_map) > 1]

    # Create plot
    fig = px.scatter(data_frame=df,
                     x='theta',
                     y='Y',
                     color='label',
                     hover_data=extra,
                     hover_name="label",
                     labels={
                         "Y": "Intensity (a.u)",
                         "theta": r"$\text{2}\theta \text{ (deg.)}$",
                         "label": title
                     },
                     title=os.path.basename(path),
                     width=width,
                     height=height)

    # Setup hover template for main plot
    t = fig.data[0].hovertemplate
    hover_t = "<b>%{hovertext}</b><br><br>2theta=%{x}<br>Rwp=%{customdata[0]}" + t.split("rwp=%{customdata[0]}")[-1]
    fig.data[0].hovertemplate = hover_t
    fig.data[1].hovertemplate = hover_t
    fig.data[2].hovertemplate = hover_t

    # Set up mode of ycalc and ydiff
    # Line
    fig.data[1].mode = "lines"
    fig.data[2].mode = "lines"

    # Line color
    fig.data[1].marker['color'] = "red"
    fig.data[2].marker['color'] = "blue"

    # Set up open markers for exp data
    fig.data[0].marker['color'] = "black"
    fig.data[0].marker['symbol'] = "circle-open"

    # plot bragg positions
    for i, df in enumerate(s_info):
        hover_base = "theta: %{x}<br>Phase: %{text}<br>"
        colors = sns.color_palette(n_colors=len(bragg))
        theta = df.loc[:, "theta"]
        intensity = df.loc[:, "int"]
        intensity = intensity/intensity.max() # scale down
        intensity = intensity - (i+2) * offset # move down
        phase_name = df.phase_name.iat[0]  # get phase name
        if len(wt_map) > 1:  # unnecessary if only one phase
            wtfrac = wt_map[phase_name]  # get wt% for phase name
            cdata = [wtfrac] * len(theta)
            ht = hover_base + "Frac. Pct:%{customdata}<extra></extra>"
        else:
            ht = hover_base + "<extra></extra>"  # used to get rigged of the extra name
            cdata = [[]]
       # y_dummy = np.zeros(len(theta)) - (i + 2) * offset  # set up position of braggs to plot
        marker = {
            "symbol": "line-ns",
            "line_width": 1,
            "size": 12,
            "line": {"color": colors.as_hex()[i]}
        }
        fig.add_scatter(x=theta,
                        y=intensity,
                        mode="lines",
                        text=[phase_name] * len(theta),
                        name=phase_name,
                        customdata=cdata, marker=marker,
                        hovertemplate=ht)
        fig_name = os.path.basename(path).replace(".gpx", ".html")

        fig.update_layout(modebar_add=to_add)
        if outfolder == "curr":  # if current folder
            outfolder = os.getcwd()
        if save:
            fig.write_html(os.path.join(outfolder, fig_name), include_mathjax="cdn")

    return fig
