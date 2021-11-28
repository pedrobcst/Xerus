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

# This file defines the data extensions readable by the XRay project
# The main purpose is to be able to create extensive readers to support other file formats
# Keep in mind the that the file extension is limited by the GSAS-II: It should be in the format that is supported by it
# Here, this is just to provide the xray data so it the correlation between the simulated patterns and the experimental
# data can be calculated.
# Add new supported logic for reading files here
#
import pandas as pd
import numpy as np
import os


class DataExtension:
    """
    Defines a DataExtension object

    Parameters
    ----------
    xrd_df : A Pandas DataFrame with two columns: "theta" -> 2Theta angle. "int" -> Intensity
    tmin : Minimum measured angle of experimental pattern
    tmax : Maximum measured angle of experimental pattern
    step : Two-theta measurement step.


    """
    def __init__(self, xrd_df: pd.DataFrame = None, tmin: float = None, tmax: float = None, step: float = None):
        self.xrd_df = xrd_df
        self.tmin = tmin
        self.tmax = tmax
        self.step = step


class RASFile(DataExtension):
    """
    Read a Rigaku .RAS File (Japanase Version)
    Probably works for other .RAS versions files.
    """
    def __init__(self, fmt="ras"):
        super().__init__()
        self.format = fmt

    def get_theta(self, s: str) -> float:
        """
             Helper function to get the intensity out of a string
             Parameters
             ----------
             s : A string where in the format "theta intensity"

             Returns
             -------
             The two-theta angle converted into float.
        """
        return float(s.split(' ')[0])


    def get_int(self, s: str) -> float:
        """
        Helper function to get the intensity out of a string
        Parameters
        ----------
        s : A string where in the format "theta intensity"

        Returns
        -------
        The intensity converted into float.
        """
        return float(s.split(' ')[1])

    def read_data(self, file: os.PathLike) -> None:
        """
        Implements the logic for reading a .RAS Rigaku file.
        WARNING: This is for a SINGLE measurement file ONLY. It wont work if your file have multiple patterns.
        Please seperate different measurements into different files, or modify this function to accomodate such files.
        Parameters
        ----------
        file : Path to a Rigaku .RAS file

        Returns
        -------
        None

        """
        xrd = pd.read_table(file, encoding="shift-jis")
        xrd.columns = ['rd']
        xrd['ok'] = xrd.rd.apply(lambda s: s[0].isdigit())
        xrd = xrd[xrd['ok']]
        xrd = xrd.reset_index(drop=True)
        xrd['theta'] = xrd.rd.apply(self.get_theta)
        xrd['int'] = xrd.rd.apply(self.get_int)
        xrd = xrd.drop(labels=['rd', 'ok'], axis=1)
        xrd['filename'] = file
        self.xrd_df = xrd
        self.tmin = xrd.theta.min()
        self.tmax = xrd.theta.max()
        self.step = np.round(xrd.theta.iat[1] - xrd.theta.iat[0], 4)

class XYFile(DataExtension):
    """
    Reads XYFile
    A XYFile is defined as a file where the two-theta and intensity data are seperated by a space.
    We assume that XYFile has no header.
    If your file has a header, please remove it.
    """
    def __init__(self, fmt="xy"):
        super().__init__()
        self.format = fmt

    def read_data(self, file: os.PathLike) -> None:
        """
        Implements the logic for reading a .xy file
        Parameters
        ----------
        file : The path for a .xy file

        Returns
        -------
        None

        """
        xrd = pd.read_table(file, sep=" ", names=["theta", "int"])
        xrd['filename'] = os.path.basename(file)
        self.xrd_df = xrd
        self.tmin = np.round(xrd.theta.min(), 2)
        self.tmax = np.round(xrd.theta.max(), 2)
        self.step = np.round(xrd.theta.iat[1] - xrd.theta.iat[0], 4)
        print("Sucessfuly read datafile {}".format(file))


class CSVFile(DataExtension):
    """
    Reads a Comma Separeted Value File.
    If your file has a header, please remove it.
    """

    def __init__(self, fmt="csv"):
        super().__init__()
        self.format = fmt

    def read_data(self, file: str) -> None:
        """
        Implements the logic for reading a .txt Rigaku Rint-TTR 3 file
        Parameters
        ----------
        file : The path for a CSV file

        Returns
        -------
        None

        """
        xrd = pd.read_csv(file, names=["theta", "int"])
        xrd['filename'] = os.path.basename(file)
        self.xrd_df = xrd
        self.tmin = np.round(xrd.theta.min(), 2)
        self.tmax = np.round(xrd.theta.max(), 2)
        self.step = np.round(xrd.theta.iat[1] - xrd.theta.iat[0], 4)
        print("Sucessfuly read datafile {}".format(file))

