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
from pymatgen import Composition
from typing import List, Tuple
import json
import re
import pandas as pd

# Regex for finding numbers
number = r'-?\d+\.?\d*'
plotly_add = [
    "drawline",
    "togglespikelines",
    "eraseshape",
    "hovercompare",
    "hoverclosest",
    "drawopenpath"
]
def blockPrinting(func):
    """
    Decorator for blocking prints of functions
    Credits for Fowler @ stackoverflow
    This is used to supress GSAS-II massive refinement messages as it does not have an option for verbose.

    Parameters
    ----------
    func : A function to suppress the print outputs

    Returns
    -------

    """
    def func_wrapper(*args, **kwargs):
        # keep old stdout
        old_stdout = sys.stdout
        # block all printing to the console
        sys.stdout = open(os.devnull, 'w')
        # call the method in question
        value = func(*args, **kwargs)
        # enable all printing to the console
        sys.stdout = old_stdout
        # pass the return value of the method back
        return value
    return func_wrapper

def formula_to_latex(formula: str) -> str:
    """
    Transforms a formula string such as Er3Ni2, HoB2 etc into a LaTeX format
    eg. Er3Ni2 -> Er$_{3}$Ni$_{2}$ etc.
    Parameters
    ----------
    formula : A chemical formula in correct format (eg. HoB2 not Hob2 or hob2 or hoB2)

    Returns
    -------
    A LaTeX-style string of the chemical formula. Ie: HoB2 -> HoB$_{2}$
    """
    ans = ''
    inputs = re.findall(r'\D+', formula.replace('.', ''))
    numbers = re.findall(number, formula)
    if len(numbers) == 0:
        return formula
    for piece, n in zip(inputs, numbers):
        ans += piece + "$_{" + n + "}$"
    if len(inputs) > len(numbers):
        ans += inputs[-1]
    return ans


def create_folder(path: str) -> None:
    """
    Helper function for creating a folder in case it does not exist.
    Parameters
    ----------
    path : Folder to create

    Returns
    -------
    Nothing
    """
    
    if not os.path.isdir(path):
        try:
            os.mkdir(path)
        except OSError: # catch problem due to paralellization
            pass


def load_json(path: str) -> object:
    """
    Opens a JSON file and returns a dictionary.
    Parameters
    ----------
    path : Path to a JSON file

    Returns
    -------
    The dictionary representation of the file
    """
    with open(path, "r") as fp:
        return json.load(fp)

def save_json(dic: dict, outpath: str) -> None:
    """
    Saves a dictionary into a JSON object
    Parameters
    ----------
    dic : A dictionary to be saved
    outpath : Path to be saved

    Returns
    -------
    Nothing, saves a JSON file into a path.
    """
    with open(outpath, "w") as fp:
        json.dump(dic, fp, indent=4)


def group_data(df: pd.DataFrame, column_name: str, reset: bool=True) -> List[pd.DataFrame]:
    """
    Auxiliary function for grouping by a columns in a pandas dataframe and getting back those groups as a list of dataframes.
    Parameters
    ----------
    df : A dataframe
    column_name : The column(s) to groupby.
    reset : To reset the index or not. Defaults to True.

    Returns
    -------
    A list of dataframes grouped by the desired column.
    """
    group = df.groupby(by=column_name)
    dfs = []
    for key in group.groups:
        data_1 = group.get_group(key)
        if reset:
            data_1 = data_1.reset_index().drop(labels=["index"], axis=1)
        dfs.append(data_1)
    return dfs


def normalize_formula(formula: str) -> str:
    """
    Truncate a formula composition to highest value.
    Ie: HoB2.6 -> HoB3.
    formula : A chemical formula

    Returns
    -------
    A truncated standarized formula that was truncated up.
    """
    import numpy as np
    comp = Composition(formula).fractional_composition.as_dict()
    return ''.join([str(key) + str(np.round(value, 1)) for key, value in comp.items()])


def make_offset(patterns: List[pd.DataFrame], offset: int = 10, imult: float = 1, normalize: bool = True) -> None:
    """
    Helper function to make offset for simulated patterns.
    Parameters
    ----------
    patterns : A list of pd.DataFrames containing the patterns
    offset : Distance between each pattern
    imult : Intensity multiplier
    normalize : To normalize by maximum or not

    Returns
    -------
    Nothing. Modifies the list of dataframes inplace.
    """
    for i, df in enumerate(patterns):
        if normalize:
            df['int'] = df['int'] / df['int'].max()
        df['int'] = df['int'] - (i + 1) * offset
        df['int'] = df['int'] * imult
        df.sort_values(by='theta', inplace=True)


def sort_and_reset(df: pd.DataFrame, by: str ="Cij") -> pd.DataFrame:
    """
    Auxiliary function to sort and reset index of a dataframe.

    Parameters
    ----------
    df : A pandas dataframe
    by : The column to sort by and then reset. Defaults to "Cij"

    Returns
    -------
    A pd.DataFrame
    """
    df = df.sort_values(by=by, ascending=False)
    df = df.reset_index(drop=True)
    return df
