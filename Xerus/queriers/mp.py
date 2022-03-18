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
import pymatgen
from pymatgen.io.cif import CifWriter
from Xerus.utils.cifutils import make_combinations
from Xerus.settings.settings import MP_API_KEY
import os
import pandas as pd
from typing import Tuple, List
api_key = MP_API_KEY
dump_folder = "mp_dump/"

def make_cifs(data: pd.DataFrame, symprec: float = 1e-2, folder_path: os.PathLike = dump_folder) -> None:
    """
    Write CIFs from data queried from Materials Project into a folder

    Parameters
    ----------
    data : A pandas dataframe containing the material id, formula and structure information
    symprec : Spglib symprec tolerance to find spacegroups (read pymatgen documentation for more detail)
    folder_path : Folder to save.

    Returns
    -------
    None
    """
    if not os.path.isdir(folder_path):
        os.mkdir(folder_path)
    for mid, name, struc in zip(data['material_id'], data['pretty_formula'], data.structure):
        writer = CifWriter(struc, symprec=symprec)
        filename = folder_path + os.sep + name + "_" + "MP_" + mid + ".cif"

        if name != 'O2':
            writer.write_file(filename)
            print("Writing {}".format(filename))


def querymp(inc_eles: List[str], exc_eles: List = [], max_num_elem:int = 3, api_key: str = api_key,
            combination: bool = True, min_hull: float = 1e-4, write: bool =True,
            folder_path: os.PathLike = dump_folder) -> pd.DataFrame:
    '''

    Parameters
    ----------
    inc_eles : Elements to query (inclusive)
    exc_eles : Elements NOT to query (exclusive)
    max_num_elem : Maximum number of elements
    api_key : API-Key
    combination : boolean. To combine or not the list of elements into a all possible system-types.
    min_hull : Minimum e_above_hull to consider when downloading CIFS
    write : To save the cifs into a folder or not
    folder_path : If write = True, which folder to save.

    Returns
    -------
    Returns a DataFrame with the queried data information with data is available for elements combination.
    '''
    a = pymatgen.MPRester(api_key)
    if not combination:
        qp = {"elements": {"$all": inc_eles, "$nin": exc_eles}, "nelements": {"$lte": max_num_elem}}
        data = a.query(qp, ['pretty_formula', 'structure', 'theoretical', 'material_id', 'e_above_hull'])
        if len(data) > 0:
            data = pd.DataFrame.from_dict(data)
            data = data[~data.theoretical]
            data = data[data.e_above_hull <= min_hull]
            data.reset_index(drop=True, inplace=True)
    else:
        dfs = []
        combations_flat = make_combinations(inc_eles, max_num_elem)
        for comb in combations_flat:
            print('Getting data for the following atoms combinations: {}'.format('-'.join(comb)))
            qp = {"elements": {"$all": comb, "$nin": exc_eles}, "nelements": {"$lte": len(comb)}}
            data = a.query(qp, ['pretty_formula', 'structure', 'theoretical', 'material_id', 'e_above_hull'])
            if len(data) > 0:
                data = pd.DataFrame.from_dict(data)
                data = data[~data.theoretical]
                data.reset_index(drop=True, inplace=True)
                dfs.append(data)
        final = pd.concat(dfs)
        final = final[final.e_above_hull <= min_hull]
        final.reset_index(drop=True, inplace=True)
        return final

    if len(data) == 0:
        return 'No data.'
    if write:
        make_cifs(data, folder_path=folder_path)
    return data
