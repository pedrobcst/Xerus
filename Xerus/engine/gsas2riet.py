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
from typing import List, Tuple, Union

import pandas as pd
from Xerus.engine.gsas2utils import get_wtfract, parse_reflections_array
from Xerus.settings.settings import GSAS2_BIN, INSTR_PARAMS
from Xerus.utils.tools import blockPrinting

sys.path.insert(0, os.path.expanduser(GSAS2_BIN))
import GSASIIscriptable as G2sc


def HistStats(gpx: G2sc.G2Project) -> float:
    """
    Taken from GSAS II Tutorial webpage.
    Returns Rwp (%) from G2Project

    Parameters
    ----------
    gpx : a G2Project

    Returns
    -------
    float
        Returns Rwp(%) of a refinement step.
    """
    print(u"*** profile Rwp, " + os.path.split(gpx.filename)[1])
    for hist in gpx.histograms():
        try:
            print("\t{:20s}: {:.2f}".format(hist.name, hist.get_wR()))
        except TypeError:
            return 999999  # Huge rwp
    print("")
    gpx.save()
    return hist.get_wR()


@blockPrinting
def simulate_pattern(filename: str, instr_params:str = INSTR_PARAMS, tmin: float = 10.0, tmax: float = 70.0,
                     step: float = 0.02, scale: float = None, threshold: float = 0.3) -> int:
    """
    Helper function for testing simulation of obtained CIFS from providers.
    Simulates a pattern and check if its able to export details.

    Parameters
    ----------
    filename : str
        Path to a Crystallography Information File (CIF).
    instr_params : str
        Instrument parameters for X-Ray machine.
    tmin : float, default: 10.0
        Minimum 2theta.
    tmax : float, default: 70.0
        Maximum 2theta
    step : float, default: 0.02
        Two-theta step.
    scale : float, default: None
        Scale
    threshold : float, default: 0.03
        Threshold for testing if there is any peaks in the simulation.
        Rarely, some CIFs provides simulations with no peaks, that is, ysim = 1 for every theta.
        In this case, we try to "quickly" check by seeing if there is any flutation (std).
        In case there is none, this will be considered "bad" CIF.
        This is just a placeholder trial, and should be changed in the future)

    Returns
    -------
    int
        Returns simulation status: 1 if success, -1: if failed.
    """
    if GSAS2_BIN not in sys.path:
        sys.path.insert(0, os.path.expanduser(GSAS2_BIN))
    status = 1
    gpx = G2sc.G2Project(filename="dummy.gpx")
    try:
        phase0 = gpx.add_phase(filename, phasename="dummy", fmthint="CIF")
    except:
        status = -1
    hist1 = gpx.add_simulated_powder_histogram("dummy",
                                               iparams=instr_params,
                                               Tmin=tmin,
                                               Tmax=tmax,
                                               Tstep=step,
                                               phases=gpx.phases(),
                                               scale=scale)
    # gpx.histograms()[0].data['Instrument Parameters'][0]['Zero'] = 0
    try:
        gpx.do_refinements([{}])
    except:
        status = -1

    try:
        y = gpx.histogram(0).getdata("ycalc")
        gpx.histogram(0).Export("dummy", ".csv", "refl")  # export reflections to dummy
        if y.std() < threshold:
            status = -1
    except (KeyError, IndexError) as e:
        return -1  
    return 1


def run_gsas_mp2(powder_data: str, cifs: Tuple[str], phasenames: Tuple[str], outfolder: str,
                 instr_params: str = INSTR_PARAMS, return_gpx: bool = False, ori: bool = False) -> Tuple[float, str, List[int]]:
    """
    Run GSAS Refinement based on tutorial for multi-phase.


    Parameters
    ----------
    powder_data : str
        Path for a experimental powder XRD pattern, in a format supported by GSASII
    cifs : array-like
        A list or tuple of path to CIF files.
    phasenames : array-like
        A list or tuple of name for each phase in the CIF files
    outfolder : str
        Name of the folder to save GSASII-related files.
    instr_params : str, default: INSTR_PARAMS (config.conf)
        Instrument paramaeters for X-Ray machine
    return_gpx: bool, default: False
        If True will return the full GPX project instead of filename

    Returns
    -------
    Tuple

    """
    old_stdout = sys.stdout
    # block all printing to the console
    sys.stdout = open(os.devnull, 'w')
    # call the method in question
    # enable all printing to the console
    outfolder_new = os.path.join(outfolder, 'gsas2_files')
    if not os.path.isdir(outfolder_new):
        os.mkdir(outfolder_new)

    gpx_t = '_'.join(phasenames) + ".gpx"
    gpx_name = os.path.join(outfolder_new, gpx_t)
    gpx = G2sc.G2Project(newgpx=gpx_name)
    gpx.data['Controls']['data']['max cyc'] = 8

    hist1 = gpx.add_powder_histogram(datafile=powder_data,
                                     iparams=instr_params)
    phases = []
    hist1.data['Instrument Parameters'][0]['I(L2)/I(L1)'] = [0.5, 0.5, 0]
    for i, cif in enumerate(cifs):
        phases.append(gpx.add_phase(cif, phasename=phasenames[i], histograms=[hist1]))


    refdict0 = {"set": {"Background": {"no. coeffs": 6, "refine": True},
                        "Cell": True,
                        "Scale": True,
                        'Sample Parameters': ['DisplaceX']}}
    # refdict4a = {"set": {'Sample Parameters': ['DisplaceX', 'Shift']}}
    refdict_ori = {"set": {"Pref.Ori.": True}, "phases": phases[0]}
    refdict5c = {"set": {'Instrument Parameters': ['U', 'V', 'W']}}
    dictList = [refdict0,  refdict5c]
    if ori:
        dictList = [refdict0,  refdict_ori, refdict5c]
    gpx.do_refinements(dictList)
    wt = get_wtfract(gpx_name.replace(".gpx", ".lst"), len(cifs))
    rwp = HistStats(gpx)
    if return_gpx:
        sys.stdout = old_stdout
        return rwp, gpx, wt

    sys.stdout = old_stdout
    return rwp, gpx_name, wt

def simulate_spectra(cifpath: str, tmin: float, tmax: float, step: float, outfolder: str,
                         instr_params,
                         simul_suffix: str = '_Simul',
                         refl_suffix: str = '_Refl',
                         sim_folder: str = 'Simul', refl_folder: str = 'Refl') -> Tuple[str, str, bool]:
    """
    Simulate a pattern using GSAS II. This needs to be define inside the scope so it can be paralleized.
    Parameters
    ----------
    cifpath : Path to CIF File
    tmin : min theta
    tmax : max theta
    step : step
    outfolder : folder to save
    instr_params : instrument paremeters to use
    scale : histogram scale
    simul_suffix : simulation files suffix (csv)
    refl_suffix : reflections files suffix (csv)
    sim_folder : simulation folder
    refl_folder : reflections folder

    Returns
    -------
    A Tuple containing (path to simulations, path to reflections, simulation successful)
    """
    old_stdout = sys.stdout
    # block all printing to the console
    sys.stdout = open(os.devnull, 'w')
    # call the method in question
    # enable all printing to the console
    cifname = os.path.basename(cifpath).split(".cif")[0]

    wanted = ['h', 'k', 'l', 'd-space', '2-theta', 'F_calc']
    gpxname = cifname.replace(".", "_") + ".gpx"
    gpx = G2sc.G2Project(filename=gpxname)
    phase0 = gpx.add_phase(cifpath, phasename="dummy", fmthint="CIF")
    hist1 = gpx.add_simulated_powder_histogram("dummy",
                                               iparams=instr_params,
                                               Tmin=tmin,
                                               Tmax=tmax,
                                               Tstep=step,
                                               phases=gpx.phases(),
                                               scale=100.0)
    gpx.set_Controls('cycles', 1)
    # In case there is still any issue with the CIF / Phase that was not caught
    try:
        gpx.do_refinements([{}])
        x = gpx.histogram(0).getdata("x")  # Grab theta
        y = gpx.histogram(0).getdata("ycalc")  # Grab ntesntieis
        simulated_data = pd.DataFrame(data=zip(x, y), columns=["theta", "int"])
        # gpx.histogram(0).Export(cifname, ".csv", "refl")  # export reflections to dummy
    except:
        # self._notrun.append(cifpath)
        sys.stdout = old_stdout
        return False, False, False

    pwdr = gpx.histograms()[0]
    phase = list(pwdr.reflections().keys())[0]
    reflection_dict = pwdr.reflections()[phase]['RefList']
    reflections = pd.DataFrame(data=[parse_reflections_array(array) for array in reflection_dict], columns=wanted)

    # Messy below. Avoid get caught by an OSError due to race conditions of parallel processing.
    if not os.path.isdir(outfolder):
        try:
            os.mkdir(outfolder)
        except OSError:
            pass
    if not os.path.isdir(os.path.join(outfolder, sim_folder)):
        try:
            os.mkdir(os.path.join(outfolder, sim_folder))
        except OSError:
            pass

    if not os.path.isdir(os.path.join(outfolder, refl_folder)):
        try:
            os.mkdir(os.path.join(outfolder, refl_folder))
        except OSError:
            pass
    filename = os.path.basename(cifpath)
    simul_filename = filename.split(".cif")[0] + simul_suffix + ".csv"
    refl_filename = filename.split(".cif")[0] + refl_suffix + ".csv"

    simul_out = os.path.join(outfolder, sim_folder, simul_filename)
    refl_out = os.path.join(outfolder, refl_folder, refl_filename)

    simulated_data.to_csv(simul_out, index=False)
    reflections.to_csv(refl_out, index=False)

    sys.stdout = old_stdout
    return simul_out, refl_out, True


<<<<<<< HEAD
<<<<<<< HEAD
@blockPrinting
=======

>>>>>>> Changed testing to cif reading only.
=======
@blockPrinting
>>>>>>> readd blockprinting decorator
def run_gsas(powder_data: str, cif_name: str, phasename: str,
             limits: Tuple[float, float] = (10., 70.) ,max_cyc: int = 8,
             instr_params: str = INSTR_PARAMS) -> Tuple[float, str]:
    """
    Run GSAS Refinement based on tutorial for ONE phase.
    Used for "testing" CIFs caming from database for refinement.

    Parameters
    ----------
    powder_data : str
        Path for a experimental powder XRD pattern, in a format supported by GSASII.
    cif_name : str
        Path to a CIF file.
    phasename : str
        Name for the phase in the CIF file.
    limits : array-like, default: (10.0, 70.0)
        A tuple containing (min theta, max theta) for refinement.
    max_cyc : int
        Max refinement cycles
    instr_params : str, default: INSTR_PARAMS (config.conf)
        Path to instrument parameters.

    Returns
    -------
    array-like
        Returns a tuple containing (Rwp, gpx_name)
    """
    print("Loading powder data: {} \n".format(powder_data))
    print("Using CIF: {} \n".format(cif_name))
    print("Phasename: {}".format(phasename))
    gpx_name = cif_name.replace(".cif", ".gpx")
    gpx = G2sc.G2Project(newgpx=gpx_name)
    gpx.data['Controls']['data']['max cyc'] = max_cyc
    hist1 = gpx.add_powder_histogram(datafile=powder_data,
                                     iparams=instr_params)
    try:                                
        phase0 = gpx.add_phase(cif_name, phasename=phasename, histograms=[hist1])
        hist1.data['Instrument Parameters'][0]['I(L2)/I(L1)'] = [0.5, 0.5, 0]
    except:
        return -1
    return 1

    # # Set to use iso
    # #for val in phase0.data['Atoms']:
    # #    val[9] = 'I'
    # refdict0 = {"set": {"Background": {"no. coeffs": 6, "refine": True}, "Scale": True}}
    # refdict1 = {"set": {"Cell": True, 'Sample Parameters': ['DisplaceX']}}
    # # refdict4a = {"set": {'Sample Parameters': ['DisplaceX']}}
    # refdict_ori = {"set": {"Pref.Ori.": True}}
    # refdict4b = {"set": {"Atoms": {"all": "XU"}}}
    # refdict5a = {"set": {'Limits': limits}}
    # refdict5c = {"set": {'Instrument Parameters': ['U', 'V', 'W']}}

    # dictList = [refdict0, refdict1, refdict_ori, refdict5a, refdict4b, refdict5c]
    # gpx.do_refinements(dictList)
    # #os.remove(gpx_name)
    # return HistStats(gpx), gpx_name


@blockPrinting
def quick_gsas(powder_data: str, cif_name: str, phasename: str, outfolder: str, max_cyc: int = 8,
               instr_params: str = INSTR_PARAMS, name: bool = False, ori: bool = False) -> Union[Tuple[float, str], str]:
    """
    Run a quick GSAS Refinement based on tutorial for ONE phase.
    To avoid any issue while refining the candidates structures according to the given pattern,
    only background, sample parameters and instrument parameters are refined.

    Parameters
    ----------
    powder_data : str
        Path for a experimental powder XRD pattern, in a format supported by GSASII.
    cif_name : str
        Path to a CIF file.
    phasename : str
        Name for the phase in the CIF file.
    outfolder: str
        Folder to save GSASII-related files
    max_cyc : int
        Max refinement cycles
    instr_params : str, default: INSTR_PARAMS (config.conf)
        Path to instrument parameters.
    name: bool, default: False
        If `name` is True, it this function will only return the gpx_name
    ori: bool, default: False
        If `ori` is set to True, prefered orientation will be also refined.

    Returns
    -------
    array-like
        Returns a tuple containing (Rwp, gpx_name)
    """
    outfolder_new = os.path.join(outfolder, 'gsas2_files')
    if not os.path.isdir(outfolder_new):
        os.mkdir(outfolder_new)
    print("Loading powder data: {} \n".format(powder_data))
    print("Using CIF: {} \n".format(cif_name))
    print("Phasename: {}".format(phasename))
    gpx_t = phasename + ".gpx"
    gpx_name = os.path.join(outfolder_new, gpx_t)
    gpx = G2sc.G2Project(newgpx=gpx_name)
    gpx.data['Controls']['data']['max cyc'] = max_cyc
    hist1 = gpx.add_powder_histogram(datafile=powder_data,
                                     iparams=instr_params)
    phase0 = gpx.add_phase(cif_name, phasename=phasename, histograms=[hist1])
    hist1.data['Instrument Parameters'][0]['I(L2)/I(L1)'] = [0.5, 0.5, 0]

    # Set to use iso
    # for val in phase0.data['Atoms']:
    #     val[9] = 'I'
    # refdict0 = {"set": {"Background": {"no. coeffs": 8, "refine": True},
    #                     "Cell": True,
    #                     "Instrument Parameters": ["Zero"],
    #                     "Scale": True}}
    refdict0 = {"set": {"Background": {"no. coeffs": 6, "refine": True}}}
    refdict1 = {"set": {"Cell": True, 'Sample Parameters': ['DisplaceX']}}

    # refdict4a = {"set": {'Sample Parameters': ['Shift', 'DisplaceX']}}
    # refdict4b = {"set": {"Atoms": {"all": "XU"}}}
    refdict_ori = {"set": {"Pref.Ori.": True}}
    refdict5c = {"set": {'Instrument Parameters': ['U', 'V', 'W']}}
    if ori:
        dictList = [refdict0, refdict1, refdict_ori, refdict5c]
    else:
        dictList = [refdict0, refdict1, refdict5c]
    gpx.do_refinements(dictList)
    # os.remove(gpx_name)
    if name:
        return gpx_name
    else:
        return HistStats(gpx), gpx



def refine_comb(comb: dict, data: str, wf: str) -> dict:
    """
    A auxiliary function defined in scope so it can be run in parallel by Ray
    Run refinement of obtained combinations of cif files.
    Parameters
    ----------
    comb : A combination dictionary
    data : Exp. datafile
    wf : Working folder

    Returns
    -------
    A dictionary containing the refinement results
    """

    cif_paths = [c['full_path'] for c in comb['comb']]  # Grab the cif files in a list
    names = [c['name'] + '-' + c['spacegroup'].replace("/", "-") for c in comb['comb']]
    out = "_".join(names)
    rwp, gpx, wt = run_gsas_mp2(powder_data=data,
                                cifs=cif_paths,
                                phasenames=names,
                                outfolder=wf)
    print('Phase: {0}, Rwp: {1}'.format(out, rwp))
    rdata = dict(rwp=rwp, wt=wt)
    return rdata
