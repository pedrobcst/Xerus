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
from mp_api.client import MPRester
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


def querymp(inc_eles: List[str], max_num_elem:int = 3, min_hull: float = 1e-4, write: bool =True,
            folder_path: os.PathLike = dump_folder) -> pd.DataFrame:
    '''

    Parameters
    ----------
    inc_eles : Elements to query (inclusive) eg:["Ho", "B"]
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
    api_key: str = MP_API_KEY
    properties = ['formula_pretty',  'material_id',  'structure',  'energy_above_hull',  'theoretical']
    cifs = ['pretty_formula', 'material_id', 'structure', 'e_above_hull', 'theoretical']

    with MPRester(api_key) as mpr:    
        datas = mpr.summary.search( chemsys=inc_eles,
                                    fields = properties, 
                                    theoretical = False, 
                                    energy_above_hull = ( 0, min_hull))
        datadf = [data.dict() for data in datas]
        datadf = pd.DataFrame(datadf)
        datadf = datadf.drop(columns = ['fields_not_requested'])
        datadf = datadf.set_axis(cifs, axis='columns')
        
    if len(datadf) == 0:
        return 'No data.'
    if write:
        make_cifs(datadf, folder_path=folder_path)
    return datadf
