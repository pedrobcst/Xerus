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

# Implements the data reader
# If your extension is not supported, please define it and the logic for reading at dataextensions.py
# NOTE: If your extension file is not supported by GSASII, please change it into a format that is supported by GSASII.
# In this case, we recommend transforming into a .xy file (twotheta intensity)
from .dataextensions import *
from pandas import DataFrame
from typing import Tuple
import os

readers = {
    "xy": XYFile(),
    "ras": RASFile(),
    "csv": CSVFile()
} # supported formats

class DataReader:
    """
    Implements the data reader to handle differents types of experimental data formats.
    """
    def __init__(self):
        self.readers = readers

    def read_data(self, file: os.PathLike, fmt: str) -> Tuple[DataFrame, float, float, float]:
        """
        Read data for a specific implemented format

        Parameters
        ----------
        file : Path to the datafile.
        fmt : File format

        Returns
        -------
        A tuple containing a (pd.DataFrame: exp data, minimum theta, max theta and two-theta step)
        """
        if fmt not in self.readers.keys():
            print (f"Format {fmt} not implemented")
            raise NotImplementedError

        reader = self.readers[fmt]
        reader.read_data(file)
        return reader.xrd_df, reader.tmin, reader.tmax, reader.step
