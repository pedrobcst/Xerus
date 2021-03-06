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
from typing import Any, List, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from plotly.graph_objects import Figure
from Xerus.utils.tools import plotly_add as to_add


def make_plot_step(runs_info: List[Tuple[pd.DataFrame, list, float, Any, Any]], topn: List[pd.DataFrame],
                   sample_name: str,
                   outfolder: os.PathLike, solver: str = "box") -> List[Figure]:
    """
    Plot each step of the correlation + peak removal in a Plotly interactive Figure.

    Note: Even if we plot the "filter" version, the bragg positions shown will always be related to the unfiltered
    version.

    Parameters
    ----------
    runs_info : A convoluted list of tuples generated by the correlation analysis process.
    topn : The top N patterns to plot
    sample_name : Sample name. Used for filename purposes only.
    outfolder : Folder to save
    solver: Used solver.

    Returns
    -------
    A list of plotly figures for each step.
    """
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    delta = runs_info[0][2]
    raw_data = runs_info[0][0]
    reference_index = set(raw_data.index.tolist())
    rem_indexes = [tuple(reference_index.difference(set(tup[0].index.tolist()))) for tup in
                   runs_info]  # Set up the removed index for plotting purposes.
    if solver == "box":
        final_exp = []
    else:
        # These are the residuals
        final_exp = [d[0] for d in runs_info]
    acc_patterns = []

    # Read up top g simulated patterns
    for n, df in enumerate(topn):
        patterns = [pd.read_csv(path) for path in df.simulated_files]
        filenames = df.filename.tolist()
        crystal_system = df.crystal_system.tolist()
        spacegrouo = df.spacegroup.tolist()
        spacegroup_number = df.spacegroup_number.tolist()
        # Save Cij and Rwp for Rietv. Method
        if solver == "rietveld":
            cij = df.Cij.tolist()
            rwp = df.rwp.tolist()
            gpx = df.gpx.tolist()

        for k, simulation in enumerate(patterns):
            simulation.int = simulation.int / simulation.int.max()  # normalized
            simulation.int = simulation.int - 1 * (k + 1)
            simulation['name'] = filenames[k]
            simulation['spg'] = spacegrouo[k]
            simulation['spgn'] = spacegroup_number[k]
            simulation['csys'] = crystal_system[k]
            # Set Riet method.
            if solver == "rietveld":
                simulation['cij'] = cij[k]
                simulation['rwp'] = rwp[k]
                simulation['gpx'] = gpx[k]
        acc_patterns.append(patterns)

    if solver == "box":
        try:
            for n, (removed, tup) in enumerate(zip(rem_indexes, runs_info)):
                print(len(removed))
                exp2_data = raw_data.copy()
                exp2_data.int = exp2_data.int / exp2_data.int.max()
                exp2_data.loc[removed, 'int'] = None
                final_exp.append(exp2_data)
                list_of_patterns = acc_patterns[n]
                for pattern in list_of_patterns:
                    pattern.loc[removed, 'int'] = None
        except KeyError:
            return "Failed to plot the step data.."

    exp_datan = raw_data.copy()
    exp_datan.int = exp_datan.int / exp_datan.int.max()
    f = []

    for n in range(len(runs_info)):
        if solver == "rietveld":
            gpx = topn[n].gpx.iat[0]
            ycalc = gpx.histograms()[0].getdata("ycalc")
            ycalc = ycalc / ycalc.max()
            theta = final_exp[n].theta.tolist()
            fit_df = pd.DataFrame(data=zip(theta,ycalc), columns=['theta', 'int'])
        dplot = final_exp[n]
        bragg = runs_info[n][1]
        y_bragg = np.zeros(len(bragg)) - 0.2
        error_bragg = np.zeros(len(bragg)) + delta / 2
        list_of_patterns = acc_patterns[n]

        if solver == "box":
            fig = px.line(data_frame=dplot, x='theta', y='int', height=768, width=1280)
            fig.add_scatter(x=exp_datan['theta'],
                            y=exp_datan['int'],
                            mode='lines',
                            line=dict(color="black", width=2, dash="dash"),
                            name="Removed")

            fig.add_scatter(x=bragg,
                            y=y_bragg,
                            error_x=dict(array=error_bragg),
                            mode="markers",
                            name="Bragg")

        elif solver == "rietveld":
            fig = px.line(data_frame=fit_df, x='theta', y='int', height=768, width=1280)
            if n > 0:
                dplot['int'] = dplot['int'] / dplot.int.max() - 0.3
                fig.add_scatter(x=dplot['theta'],
                                y=dplot['int'],
                                mode='lines',
                                line=dict(color="blue"),
                                name="Residual")
            fig.add_scatter(x=exp_datan['theta'],
                            y=exp_datan['int'],
                            mode="markers",
                            marker=dict(color="black", symbol="circle-open"),
                            name="Exp. Data")

        for sim in list_of_patterns:
            fig.add_scatter(x=sim["theta"],
                            y=sim["int"],
                            mode="lines",
                            name=sim["name"].iat[0], customdata=sim.loc[:, ['spg', 'spgn', 'csys']],
                            text=sim['name'],
                            hovertemplate=
                            "<b>%{text}</b><br>" +
                            "Spacegroup: %{customdata[0]}<br>" +
                            "Spacegroup N: %{customdata[1]}<br>" +
                            "Crystal System %{customdata[2]}<br>" +
                            "<extra></extra>"
                            )

        fig.update_traces(connectgaps=False)
        if solver == "box":
            name = "Exp. Data"
        else:
            name = "Fit"
        fig.data[0]['name'] = name
        fig.data[0]['line']['color'] = "red"
        fig.data[0]['showlegend'] = True
        fig.update_layout(
            title="N = " + str(n),
            font=dict(
                family="Arial",
                size=12
            ),
            modebar_add=to_add
        )
        fig.update_layout(legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ))
        final_name = sample_name + "_N=" + str(n) + "_" + "run.html"
        outpath = os.path.join(outfolder, final_name)
        fig.write_html(outpath)
        f.append(fig)


def make_plot_all(runs_info: List[Tuple[pd.DataFrame, list, float, Any, Any]], topn: List[pd.DataFrame],
                  sample_name: str,
                  outfolder: os.PathLike) -> Figure:
    """
    Make ONE plot of all patterns tried in the correlation analysis run.

    Parameters
    ----------
    runs_info : A convoluted list of tuples generated by the correlation analysis process.
    topn : The top N patterns to plot
    sample_name : Sample name. Used for filename purposes only.
    outfolder : Folder to save

    Returns
    -------
    ONE Figure with all simulated (tried) patterns for the run.

    """
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)

    raw_data = runs_info[0][0]
    acc_patterns = []
    k = 1
    for n, df in enumerate(topn):
        patterns = [pd.read_csv(path) for path in df.simulated_files]
        filenames = df.filename.tolist()
        crystal_system = df.crystal_system.tolist()
        spacegrouo = df.spacegroup.tolist()
        spacegroup_number = df.spacegroup_number.tolist()
        for j, simulation in enumerate(patterns):
            simulation.int = simulation.int / simulation.int.max()  # normalized
            simulation.int = simulation.int - 0.3 * k
            simulation['name'] = filenames[j]
            simulation['spg'] = spacegrouo[j]
            simulation['csys'] = crystal_system[j]
            simulation['spgn'] = spacegroup_number[j]
            simulation['n'] = n + 1
            k += 1
        acc_patterns.append(patterns)
    raw_data.int = raw_data.int / raw_data.int.max()
    raw_data['name'] = 'Exp. Data'
    raw_data['spg'] = 'Undef.'
    raw_data['csys'] = 'Undef.'
    raw_data['n'] = 0
    hover_html = "<b>%{text}</b><br>" + \
                 "Run Number:%{customdata[0]}<br>" + \
                 "Spacegroup:%{customdata[1]}<br>" + \
                 "Spacegroup N:%{customdata[2]}<br>" + \
                 "Crystal system:%{customdata[3]}<br><extra></extra>"
    fig = px.scatter(data_frame=raw_data,
                     x="theta",
                     y="int",
                     height=768,
                     width=1280)

    for pattern_list in acc_patterns:
        for sim in pattern_list:
            fig.add_scatter(
                x=sim["theta"],
                y=sim["int"],
                mode="lines",
                name=sim["name"].iat[0],
                text=sim['name'],
                customdata=sim.loc[:, ['n', 'spg', 'spgn', 'csys']],
                hovertemplate=hover_html,
            )
            fig.data[0]['name'] = 'Exp. Data.'
            fig.data[0].marker.color = "black"
            fig.data[0].marker.symbol = 'circle-open'
            fig.data[0]['showlegend'] = True

    n_colors = len(fig.data) - 1
    colors = sns.color_palette(palette="dark", n_colors=n_colors).as_hex()
    for i, entry in enumerate(fig.data[1:]):
        entry['line']['color'] = colors[i]

    fig.update_layout(modebar_add=to_add)

    final_name = sample_name + "All_N_Runs.html"
    outpath = os.path.join(outfolder, final_name)
    fig.write_html(outpath, include_plotlyjs="cdn")
    print(outpath)
    return fig
