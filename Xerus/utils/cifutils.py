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
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from tqdm import tqdm
import glob
from glob import glob
from itertools import combinations
import shutil
from typing import List, Tuple, Dict
from pymatgen.core import Composition
from zipfile import ZipFile
from typing import List
import re
import json
dump_folders = ["mp_dump/", "cod_dump/", "aflow_dump/"]
test_folder = "queried_cifs/"
r = r'[a-zA-Z]+'


def movecifs(dump_folders: List[os.PathLike] = dump_folders, test_folder: str = test_folder) -> None:
    """
    Auxiliary function to move queried cifs from downloaded folder to the folder for testing
    Parameters
    ----------
    dump_folders : List of download folders
    test_folder : Folder for testing

    Returns
    -------
    Nothing
    """
    if not os.path.isdir(test_folder):
        os.mkdir(test_folder)

    for folder in dump_folders:
        if os.path.isdir(folder):
            files = os.listdir(folder)
            for file in files:
                shutil.move(os.path.join(folder,file), test_folder)
            shutil.rmtree(folder)

def write_cif(filename: str, CIF: str, outfolder: os.PathLike) -> None:
    """
    Auxilia function for writing a CIF string into a folder
    Parameters
    ----------
    filename : Filename to save
    CIF : A CIF string or byte-like CIF string
    outfolder : Folder to save

    Returns
    -------
    Nothing. Writes a CIF to a folder.
    """
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    with open(os.path.join(outfolder,filename), "w") as f:
        f.write(CIF)


def make_combinations(element_list: List[str], max_num_elem: int) -> List[Tuple[str]]:
    """
    Auxiliary function to make a flat list of elements combinations..
    Parameters
    ----------
    element_list : List of elements to be combined
    max_num_elem : Maximum combination lenght.

    Returns
    -------
    A flat list of tuples of each combination. ie:
    ['Ho', 'B'] -> [('Ho'), ('B'), ('Ho', 'B')] (I think)

    """
    comb = [list(combinations(element_list, x)) for x in range(1, max_num_elem + 1)]
    combinations_flat = [item for sublist in comb for item in sublist]
    return combinations_flat

def get_combinations_oxide(iterable, length):
    return [iterable[i:i+length] for i in range(len(iterable) - length + 1)]

def make_combinations_oxide(element_list: List[str], max_num_elem: int) -> List[Tuple[str]]:
    """
    Function to make a list of elements combinations for oxides
    Parameters
    ----------
    element_list : List of elements to be combined
    max_num_elem : Maximum combination lenght.

    Returns
    -------
    A list of each combination. ie:
    ['Fe','La','Sr','Ba','O'] ->[['O', 'Fe'], ['O', 'La'], ['O', 'Sr'], 
    ['O', 'Ba'], ['O', 'Fe', 'La'], ['O', 'La', 'Sr'], ['O', 'Sr', 'Ba'], 
    ['O', 'Fe', 'La', 'Sr'], ['O', 'La', 'Sr', 'Ba'], ['O', 'Fe', 'La', 'Sr', 'Ba']]

    """
    specific_element = 'O'
    # Create a sublist without the specific element
    sublist = [x for x in element_list if x != specific_element]
    combinations_oxide = []
    # Generate combinations of different lengths
    for length in range(1, max_num_elem):
        current_combinations = get_combinations_oxide(sublist, length)
    for combo in current_combinations:
        combinations_oxide.append([specific_element,] + list(combo))
    return combinations_oxide


def rename_multicif(folder_path: str, provider: str = 'COD') -> None:
    """
    Rename CIF in a folder for the conventional way:
    MaterialName_Provider_ProviderID.CIF
    Parameters
    ----------
    folder_path : Folder with cifs to rename
    provider :  Provider

    Returns
    -------
    None
    """
    os_sep = os.sep
    for cif in tqdm(glob(folder_path + "*.cif")):
        try:
            t = CifParser(cif)
            struc = t.get_structures(primitive=False)[0]
            comp = struc.composition.reduced_formula
            splited = cif.split(os.sep)
            path = os_sep.join(splited[:-1]) + os_sep
            mid = splited[-1].replace(".cif", "")
           # print(cif, path, comp, provider, mid)
            name = path + comp + "_" + provider + "_" + mid + ".cif"
            os.rename(cif, name)
        except:
            print("Pymatgen failed to parse the following CIF: {}".format(cif))
            print("The metadata for such CIF will be unobtainable. This CIF will be removed.")
            os.remove(cif)
            continue




def get_provider(filename: str) -> Tuple[str, str, str, str]:
    """
    Helper function to parse provider of CIF.
    Please note that the CIF filename format should be in the following way:
    MaterialName_Provider_ProviderID.cif
    Parameters
    ----------
    filename : CIF Filename

    Returns
    -------
     A tuple containing (provider, mid, formula_name, raw_name)

    """
    raw_name = os.path.basename(filename)
    splits = raw_name.split("_")
    formula_name = splits[0]
    provider = splits[1]
    mid = splits[2].split(".cif")[0]  # hopefully it works!!
    return provider, mid, formula_name, raw_name


def get_elements(name: str) -> List[str]:
    """
    Parse a valid composition string and return a list of containing elements.
    Please note that a composition has to be in a valid pymatgen format
    HoB2 -> correct
    hob2 -> wrong
    Hob2 -> wrong
    etc
    Parameters
    ----------
    name : A composition-like string, eg: HoB2

    Returns
    -------
    A list containing the elements
    HoB2 -> ['Ho', 'B']
    """
    return Composition(name).chemical_system.split("-")


def make_system_types(elements: List[str], size: int = None, ceramic_mode: bool = False) -> List[str]:
    """
    Make "system-types" from a list of elements

    Parameters
    ----------
    elements : A list of elements
    size : Length of combinations to make

    Returns
    -------
    A list of combinations in the format:
    "Ho-B", "B-O" etc.

    """
    if size is None:
        size = len(elements)
    if ceramic_mode:
        return [make_system(standarize(comb)) for comb in make_combinations_oxide(elements, size)]
    else:
        return [make_system(standarize(comb)) for comb in make_combinations(elements, size)]


def make_system(comp: str) -> str:
    """
    Make a system-like string (eg. "Ho-B") from a standarized composition.
    A standarized composition is a chemical formula generated by Pymatgen Composition, ie: Composition(formula).formula
    Parameters
    ----------
    comp : A pymatgen standarized composition (Composition(comp).formula)

    Returns
    -------
    A system type, eg "Ho-B" from Ho1 B1.

    """
    return '-'.join(re.findall(r, comp))


def standarize(comb: List[str]) -> str:
    """
    Standarize a composition from a list of elements
    Parameters
    ----------
    comb : A list of elements eg ['Ho', 'B']

    Returns
    -------
    A pymatgen standarized composition, ie: ['Ho', 'B'] -> Ho1 B1.
    """
    formula = ' '.join(comb)
    return Composition(formula).formula


def get_ciflist(folder: str, cifmt: str = "r", save_json: bool = False) -> List[Dict]:
    """
    Parse more information from cif such as lattice and spacegroup.
    :param folder: path
    :return: dataframe with cif information
    """
    accumulated = []
    from pathlib import Path
    folder_path = Path(folder)
    #used to be for cif in glob(folder + "*.cif")
    for cif in folder_path.glob("*.cif"):
        try:
            cif_pg = CifParser(cif).get_structures(primitive=False)[0]
            try:
                analyzer = SpacegroupAnalyzer(cif_pg)
                space_group = analyzer.get_space_group_symbol()
                lattice_type = analyzer.get_crystal_system()
                space_group_number = analyzer.get_space_group_number()
                pymatgen_status = {"status": 1, "tested": True} #setup flags for testing
                simul_status = {"status": 0, "tested": False}
                gsas_status = {"status": 0, "tested": False}
                ran = False #setup flags for ran
            except TypeError:
                space_group = "Error"
                lattice_type = "Error"
                space_group_number = "Error"
                pymatgen_status = {"status": -1, "tested": True} #setup flags for testing
                simul_status = {"status": 0, "tested": False}
                gsas_status = {"status": 0, "tested": False}
                ran = False
        except ValueError:
            space_group = "Error"
            lattice_type = "Error"
            space_group_number = "Error"
            pymatgen_status = {"status": -1, "tested": True}  # setup flags for testing
            simul_status = {"status": 0, "tested": False}
            gsas_status = {"status": 0, "tested": False}
            ran = False
        partial = get_provider(cif)
        with open(cif, cifmt) as f:
            content = f.read()
        entry = {
            "provider": partial[0],
            "id": partial[1],
            "name": partial[2],
            "filename": partial[3],
            "spacegroup": space_group,
            "spacegroup_number": space_group_number,
            "crystal_system": lattice_type,
            "system_type": make_system(cif_pg.composition.formula),
            "elements": get_elements(partial[2]),
            "num_elem": len(get_elements(partial[2])),
            "pymatgen_status": pymatgen_status,
            "simul_status": simul_status,
            "gsas_status": gsas_status,
            "ran": ran,
            "CIF": content
        }
        accumulated.append(entry)
    if save_json:
        with open(os.path.join(folder, "cif.json"), "w") as fp:
            json.dump(accumulated, fp)
    return accumulated


def get_spgname(cifpath: str) -> str:
    """
     Auxiliary tool for helping in making custom names for GSAS II Wrapper..
    Parameters
    ----------
    cifpath : CIF path to be parsed

    Returns
    -------
    A string with name: Phase-Spacegroup
    """
    try:
            cif_pg = CifParser(cifpath).get_structures(primitive=False)[0]
            try:
                analyzer = SpacegroupAnalyzer(cif_pg)
                space_group = analyzer.get_space_group_symbol()
            except TypeError:
                space_group = "Error"
    except ValueError:
            space_group = "Error"
    name = get_provider(cifpath)[2]
    return name + "-" + space_group.replace("/", "_")

