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
from typing import List, Tuple, Union
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
project_path = str(Path(os.path.dirname(os.path.realpath(__file__))).parent) + os.sep # so convoluted..
if project_path not in sys.path:
    sys.path.append(project_path)
from Xerus.settings.settings import GSAS2_BIN
import os

sys.path.append(GSAS2_BIN)
import GSASIIscriptable as G2sc  # Ignore Warning

def get_cell_density(lstfile: str, n_phases: int) -> List[dict]:
    """
    Parse .lst file to get new cell values with errors and density from a refinement.

    Parameters
    ----------
    lstfile : str
        Path to GSASII .lst file.
    n_phases : int
        Number of phases in the .lst file

    Returns
    -------
    List[dict]
        Returns a List of dictionaries for each phase containing the crystal structure information (a,b,c,V) and associ-
        ated errors.
    """
    results = []
    with open(lstfile, "r") as file:
        data = file.read()
        k2= data.split("names :           a           b           c       alpha        beta       gamma      Volume")
        for i in range(1, n_phases + 1):
            try:
                phasename = k2[i-1].split("Result for phase:")[1].split("\n")[0].replace(" ", "")
                abc_angle_vol = k2[i].split("\n")[1].split("values:")[1]
                err = k2[i].split("\n")[2].split("esds  : ")[1]
                density = k2[i][k2[i].find("Density:"):k2[i].find("g/cm**3")].split("Density:")[-1]
                new_cell = [val for val in abc_angle_vol.split(" ") if len(val) > 0] # strip white space
                new_cell = new_cell[:3] + [new_cell[-1]] # remove angles
                err = [val for val in err.split(" ") if len(val) > 0] # strip white space

                final_results  = dict(name=phasename,
                      a = eval(new_cell[0]),
                      err_a = eval(err[0]),
                      b = eval(new_cell[1]),
                      err_b = eval(err[1]),
                      c = eval(new_cell[2]),
                      err_c = eval(err[2]),
                      V = eval(new_cell[3]),
                      err_V = eval(err[3]),
                      density = eval(density)
                                    )
            except:
                final_results = {"name": "error"}
            results.append(final_results)
    return results


def get_wtfract(lstfile: str, n_phases: int) -> List[float]:
    """
    Parse .lst file to get Weighted fraction percentange of each phase.

    Parameters
    ----------
    lstfile : str
        Path to GSASII .lst file.
    n_phases : int
        Number of phases in the .lst file

    Returns
    -------
    List[int]
        Returns a List of integers for containing the weight fraction (%) for each phase.
    """
    wtfracs = []
    with open(lstfile, "r") as file:
        data = file.read()
        for i in range(1, n_phases + 1):
            try:
                wt = eval(data.split('Weight fraction')[i].split(",")[0].split(":")[-1])
            except:
                wt = -1
            wtfracs.append(wt)
    return wtfracs


def get_plotting_data(proj: G2sc.G2Project, normalize: bool = True, offset: float = 0.3) -> Tuple[
    np.array, np.array, np.array, np.array, float]:
    """
    Helper function for extracting Yobs, Ycal and Residual from a G2Project
    Parameters
    ----------
    proj : G2Project object
        A G2Project object containing the refinement
    normalize : bool, default: True
        To normalize the data by exp_data.max() [between 0,1.0]
    offset : float, default: 0.3,
        Offset between yobs and ydiff. [ydiff = (yobs-ycalc) - offset]

    Returns
    -------
    Tuple(2theta[np.array], yobs[np.array], ycalc[np.array], ydiff[np.array])
        A tuple containing np.arrays of 2theta, yobs, ycalc, ydiff
    """
    pwdr = proj.histograms()[0]  # get pwdr object

    ## Get the data relating to exp data & fitting
    yobs = pwdr.getdata("Yobs")
    ycalc = pwdr.getdata("Ycalc")
    ydiff = pwdr.getdata("residual")
    x = pwdr.getdata("X")
    rwp = pwdr.get_wR()
    if normalize:
        scaler = yobs.max()
        yobs = yobs / scaler
        ycalc = ycalc / scaler  # We normalzie the calc based on obs.
        ydiff = (yobs - ycalc) - offset

    return x, yobs, ycalc, ydiff, rwp


def get_reflections(proj: G2sc.G2Project) -> List[pd.DataFrame]:
    """
    Extract reflections information from a G2Project object

    Parameters
    ----------
    proj : G2Project


    Methods
    -------
    Exports the information into a dummy CSV file and reads from there.

    Returns
    -------
    List[pd.DataFrame]
        A list of dataframes containing the reflection data for each phase in the provided G2Project object.
    """
    f_path = proj.data['Controls']['data']['LastSavedAs']
    save_name = os.path.basename(f_path).replace(".gpx", "")
    save_name = save_name.replace(".", "-")
    pwdr = proj.histograms()[0]  # Get the powder info
    pwdr.Export(save_name, ".csv", "refl")  # export ito to csv
    header_lines = 2  # 0,1,2
    phase_header = 3
    phases_info = len(proj.phases())
    refl_info = open(save_name + ".csv").read()  # read it
    os.remove(save_name + ".csv")  # Remove the dummy
    data = refl_info.split("\n")

    phase_names = data[phase_header + 1: phases_info + phase_header + 1]
    phase_names = [
        entry.split(",")[0].replace('"', "") for entry in phase_names
    ]  # get a phase name
    phase_projections = data[phases_info + phase_header + 2:]
    reflections_header = phase_projections[0]  # save the header
    cname = reflections_header.replace('"', "").split(",")  # transaform into oclumns
    phase_map = {
        float(i): phase_names[i] for i in range(len(phase_names))
    }  # make a map of thep hases

    reflections = parse_projections(phase_projections, cname, phase_map)
    return reflections


def parse_projections(projection_list: list, column_names: list, phase_map: dict) -> List[pd.DataFrame]:
    """
    Helper function for parsing a list of projections from multiphase project of GSAS II

    Parameters
    ----------
    projection_list : list
        A list of data from a .csv file generated by GSASII reflection
    column_names : list
        Column names of csv file
    phase_map : dict
        A map of phase number -> phase name

    Returns
    -------
    List[pd.DataFrame]
        Returns a list of pandas dataframe with the reflections information for each phase.
    """
    phases_projection = []
    phase_data = []
    for line in projection_list[1:]:  # Skip reader
        try:
            if line[0].isdigit():
                phase_data.append(line)
            else:
                if len(phase_data) >= 0:
                    pdata_parsed = [list(map(float, l.split(","))) for l in phase_data]
                    p_dataframe = pd.DataFrame(data=pdata_parsed, columns=column_names)
                    p_dataframe["phase_name"] = p_dataframe["phase #"].map(phase_map)
                    phases_projection.append(p_dataframe)  # We finished parsed the data
                    phase_data = []
        except IndexError:
            pdata_parsed = [list(map(float, l.split(","))) for l in phase_data]
            p_dataframe = pd.DataFrame(data=pdata_parsed, columns=column_names)
            p_dataframe["phase_name"] = p_dataframe["phase #"].map(phase_map)
            phases_projection.append(p_dataframe)  # We finished parsed the data
    return phases_projection


def get_plot_data(path: str, offset: float = 0.3, get_bragg: bool = True) -> Tuple[List[pd.DataFrame], np.array,
                                                                                   np.array, np.array,
                                                                                   np.array, float, dict]:
    """
    Opens a .gpx file and parses all the information needed for plotting:
        - bragg positions
        - 2theta
        - yobs
        - ycalc
        - ydiff
        - rwp
        - wt% of each phase

    Parameters
    ----------
    path : str
        Path to a .gpx file to be parsed
    offset : float, default: 0.3
        Offset of ydiff
    get_bragg : bool, default: True
        To parse the bragg positions of each phase or not.

    Returns
    -------
    Tuple
        Returns a tuple containing bragg positions, 2 theta, yobs, ycalc, ydiff, rwp and wt of each phase.
    """
    proj = G2sc.G2Project(gpxfile=path)
    if get_bragg:
        bragg = get_reflections(proj)  # get the bragg reflections
    else:
        bragg = None
    x, yobs, ycalc, ydiff, rwp = get_plotting_data(proj, normalize=True, offset=offset)
    if len(proj.phases()) == 1:
        wts = [1.0]
    else:
        wts = get_wtfract(path.replace(".gpx", ".lst"), len(proj.phases()))
    pnames = [p.name for p in proj.phases()]
    wt_map = {name: wt for name, wt in zip(pnames, wts)}

    return bragg, x, yobs, ycalc, ydiff, rwp, wt_map


def make_gpx(name: Union[List[str], str], space: Union[List[str], str], path: Union[List[str], str]) -> str:
    """
    Auxiliary function to reconstruct final gpx path out of filename spacegroup and cif path.

    Parameters
    ----------
    name : str
        Sample name
    space : List[str] or str
        Spacegroups of the phases
    path : List[str] or str
        Path to CIF of phases.

    Returns
    -------
    Final .gpx name used for refinement during combination analysis.
    """
    if type(name) == list:
        gpx_name = "_".join([n + '-' + s.replace("/", "-") for (n, s) in zip(name,space)]) + ".gpx"
        _path = os.path.dirname(path[0]).replace("cifs", "gsas2_files")
        gpx_path = os.path.join(_path,gpx_name)
        return gpx_path
    else:
        gpx_name = name + "-" + space.replace("/", "-") + ".gpx"
        _path = os.path.dirname(path).replace("cifs", "gsas2_files")
        gpx_path = os.path.join(_path, gpx_name)
    return gpx_path

def parse_reflections_array(array: np.array) -> Tuple:
    """
    Helper method for parsing reflections without exporting to ".CSV"
    Parse an reflection array of GSASII defined as:

    index	explanation
    0,1,2	h,k,l (float)
    3	(int) multiplicity
    4	(float) d-space, Å
    5	(float) pos, two-theta
    6	(float) sig, Gaussian width
    7	(float) gam, Lorenzian width
    8	(float) F2obs
    9	(float) F2calc
    10	(float) reflection phase, in degrees
    11	(float) intensity correction for reflection, this times F2obs or F2calc gives Iobs or Icalc
    12	(float) Preferred orientation correction
    13	(float) Transmission (absorption correction)
    14	(float) Extinction correction


    Parameters
    ----------
    array : np.array
        A np.array from a reflection list of A G2Pwdr Object ('ReflList')

    Returns
    ----------
    Tuple
        A tuple containing (h,k,l dspace, 20 and Icalc)
    """
    h = array[0]
    k = array[1]
    l = array[2]
    dspace = array[4]
    theta = array[5]
    icalc = array[9]*array[11]
    return (h, k,l, dspace, theta, icalc)


def parse_reflections_gpx(path: str) -> List[pd.DataFrame]:
    """
    Parse the reflections of a GPX File using the array method (without exporting .CSV)

    Parameters
    ----------
    path : str
        Path to a .GPX file.

    Returns
    -------
    List[pd.DataFrame]
        A list of pandas dataframes containing the reflections data for each phase present in the .gpx file
    """
    project = G2sc.G2Project(gpxfile=path)
    pwdr = project.histograms()[0]
    phase_names = list(pwdr.reflections().keys())
    parsed = []
    columns = [
        "h",
        "k",
        "l",
        "d-space",
        "2theta",
        "Icalc",
        "Phase"
    ]
    for phase in phase_names:
        reflection_dict = pwdr.reflections()[phase]['RefList']
        data = pd.DataFrame(data=[list(parse_reflections_array(array)) + [phase] for array in reflection_dict], columns=columns)
        data['Icalc'] = data['Icalc']/data['Icalc'].max() * 100 # normalzie 0 to 100
        parsed.append(data)
    return parsed


