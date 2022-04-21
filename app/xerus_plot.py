from copy import deepcopy
from typing import List

import pandas as pd
import plotly.express as px
from matplotlib import pyplot as plt
from peakutils import baseline
from plotly.graph_objects import Figure
from Xerus.readers.datareader import DataReader
from Xerus.settings.mplibcfg import Mplibcfg
from Xerus.utils.tools import plotly_add as to_add

from ui_utils import make_simulation_df


def plot_read_data(data: str, format: str, poly_degree: int = 8, remove_base: bool = True) -> Figure:
    df_new, _, _, _, = DataReader().read_data(data, fmt=format)
    _y = df_new.int.copy()
    background = baseline(_y, deg = poly_degree)
    df_new['background'] = background
    df_new['int_new'] = df_new.int - background
    fig = px.scatter(data_frame = df_new, x='theta', y = 'int', labels={"theta": "Theta", "int": "Intensity"}, template="presentation")
    fig.data[0].mode = "lines"
    fig.data[0].marker['color'] = 'purple'
    fig.data[0].marker['symbol'] = 'circle-open'
    fig.data[0].marker['size'] = 6
    fig.data[0].name = "Exp. Data"
    fig.data[0].showlegend = True

    if remove_base:
        fig.add_scatter(
            x=df_new.theta,
            y=df_new.background,
            name="Background",
            marker = {
                "symbol": "line-ns",
                "line_width": 1,
                "size": 12,
                "line": {"color": "blue"}
            }
        )
        fig.add_scatter(
            x=df_new.theta,
            y=df_new.int_new,
            name="Data - Background",
            marker = {
                "symbol": "line-ns",
                "line_width": 1,
                "size": 12,
                "line": {"color": "orange"}
            }
        )
        fig.data[2].marker['color'] = "orange"
        fig.data[1].marker['color'] = "blue"
    return fig


def plot_highest_correlated(
        data: str,
        format: str,
        cif_info: pd.DataFrame,
        top: int = 30,
        offset: float = 0.3,
        width: int = 1280,
        height: int = 968,
        filter_columns: List[str] = None,
        **kwargs,
    ):
        """
        Plot the top N highest correlated patterns with given experimental data using Plotly.
        Parameters
        ----------
        top : Number of patterns to plot (int), defaults to 30
        offset : Spacing between the plots, defaults to 0.3 (yoffset = y - (ki + 1) * offset) where ki is the kth pattern
        imult : Intensity multiplier. (yintensity = yoffset * imult)
        width : Figure width
        height : Figure height
        kwargs : Keyword arguments to be passed to px.line

        Returns
        -------
        self
        """
        # Set up the "simulated" dataframe in a easy way for plotting
        dummy_data = [
            ("name", "Exp.Data"),
            ("provider", "Xray machine"),
            ("mid", 0),
            ("spacegroup", "?"),
            ("spacegroup_number", "?"),
            ("crystal_system", "?"),
            ("Cij", 1),
        ]
        if filter_columns:
            cif_info = cif_info[cif_info.filename.isin(filter_columns)]
        exp_data, _, _, _, = DataReader().read_data(data, fmt=format)
        simulations = make_simulation_df(cif_info)
        col_order = simulations[0].columns
        for col, val in dummy_data:
            exp_data[col] = val
        exp_data = exp_data[col_order]
        patterns_sorted = sorted(simulations, key=lambda c: c.Cij.iat[0], reverse=True)
        top_patterns = patterns_sorted[:top]
        for i, pattern in enumerate(top_patterns):
            pattern['int'] = (pattern.int) / (pattern.int.max()) - offset*(i+1)
        # Patterns scaling based on exp data
        SCALER_MAX = exp_data.int.max()
        exp_data.int = exp_data.int / SCALER_MAX
        plotly_df = pd.concat([exp_data] + top_patterns)
        hover_data = ["spacegroup", "crystal_system", "spacegroup_number", "Cij"]
        figure = px.line(
            plotly_df,
            x="theta",
            y="int",
            color="filename",
            width=width,
            height=height,
            hover_data=hover_data,
            **kwargs,
        )
        config = dict({"scrollZoom": True})
        ## some default plotly options ##
        hoverlabel = dict(font_size=14, font_family="Arial")
        legend = {"font": dict(size=14)}

        # Edit First and Second Plot
        figure.data[0]["mode"] = "markers"
        figure.data[0].marker.color = "black"
        figure.data[0].marker.symbol = "circle-open"
        figure.data[0].name = "Exp. Data"
        figure.update_layout(hoverlabel=hoverlabel, legend=legend, modebar_add=to_add)
        
        return figure

