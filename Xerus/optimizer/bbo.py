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
from __future__ import annotations
import os
import sys
from pathlib import Path
import optuna
import shutil
import numpy as np
import pandas as pd
from itertools import combinations
from multiprocessing import Process, Queue
from Xerus.settings.settings import GSAS2_BIN
from typing import List, Union, Tuple
from Xerus.engine.gsas2viz import plot_gpx
from Xerus.engine.gsas2utils import get_wtfract, get_cell_density  # For some test
from Xerus.utils.tools import blockPrinting
from Xerus.settings.mplibcfg import Mplibcfg
from matplotlib import pyplot as plt
from plotly.graph_objects import Figure
import plotly.express as px
sys.path.append(GSAS2_BIN)
import GSASIIscriptable as G2sc  # Ignore Warning


def combine(list_of_phases: List[str]) -> List[str]:
    """
    Auxiliary tool for making a combination of phases "numbers"
    Used for Preferred Orientation and Atomic refinement of multi phase samples.

    Parameters
    ----------
    list_of_phases : A List of phases indexes (automatic)

    Returns
    -------
    A list of string representations for each phase, ie "0", "0,1" etc.
    """
    add = []
    for i in range(2, len(list_of_phases)):
        combined = combinations(list_of_phases, i)
        parsed = [','.join(c) for c in combined]
        for element in parsed:
            add.append(element)
    return list_of_phases + add + [','.join(list_of_phases)]


class OptunaTrial:
    """
    A wrapper for GSAS II Project for using in optuna blackbox optimization
    Creates new GSAS II projects based on inputs of optuna and perform refinement.

    Adapted from
    Ozaki, Y., Suzuki, Y., Hawai, T., Saito, K., Onishi, M., and Ono, K.
    Automated crystal structure analysis based on blackbox optimisation.
    npj Computational Materials 6, 75 (2020).
    https://doi.org/10.1038/s41524-020-0330-9


    Parameters
    ----------
    trial_number : Trial number
    gpx_path : Path to gpx for trial run
    working_folder : working folder
    random_state : optuna random state.
    """

    def __init__(self, trial_number, gpx_path: str, working_folder: str, random_state: int):
        """
        Initialize trial
        """
        basefolder = "optuna/"  # create a optuna folder to keep the files
        copy_path = os.path.join(working_folder, basefolder)  # setup the apth for copying

        # Check if folder exist. If doesnt create.
        if not os.path.isdir(copy_path):
            try:
                os.mkdir(copy_path)
            except OSError: # to deal with parallel issues.
                pass

        # Grab the gpx filename and setup the trial files name
        self.gpx_filename = os.path.basename(gpx_path)
        gpx_trialname = self.gpx_filename.split(".gpx")[0]
        self.trial_name = gpx_trialname + "_seed{0}_trial_{1}.gpx".format(random_state, trial_number)
        self.trial_path = os.path.join(copy_path, self.trial_name)

        # Create a project with a default project name
        shutil.copyfile(gpx_path, self.trial_path)
        self.gpx = G2sc.G2Project(gpxfile=self.trial_path)

        # Add two histograms to the project
        self.hist1 = self.gpx.histograms()[0]
        self.phase0 = self.gpx.phases()[0]
        self.hist1.data['Instrument Parameters'][0]['I(L2)/I(L1)'] = [0.5, 0.5, 0]

        # Set to use iso
       # for val in self.phase0.data['Atoms']:
        #    val[9] = 'I'

    @blockPrinting
    def refine_and_calc_Rwp(self, param_dict: dict, m: str = "rwp") -> float:
        """
        Receives a param dict, refine and return rwp

        Parameters
        ----------
        param_dict : dict
            A dictionary containing the refinement params
        m : str
            "Method". Defines which object should be returned. Defaults to "rwp". For GoF use "gof"
        """
        self.gpx.do_refinements([param_dict])
        if m == "rwp":
            for hist in self.gpx.histograms():
                _, Rwp = hist.name, hist.get_wR()
        else:
            Rwp = self.gpx.data['Covariance']['data']['Rvals']['GOF']
        return Rwp


class BlackboxOptimizer:
    """
    Blackbox optimization for Rietveld Refinement.
    Adapted from
    Ozaki, Y., Suzuki, Y., Hawai, T., Saito, K., Onishi, M., and Ono, K.
    Automated crystal structure analysis based on blackbox optimisation.
    npj Computational Materials 6, 75 (2020).
    https://doi.org/10.1038/s41524-020-0330-9

    to allow for multi phase + extra Rietveld parameters (strain, etc)

    Parameters
    ----------
    gpx_path : Base GPX Project containing experimental information and loaded CIFS
    working_folder : Folder to save the results
    random_state : Random state for optuna
    n_trials : Number of trials to perform blackbox optimization
    n_jobs : Number of jobs
    allow_pref_orient : Allow Pref. Ori. to be considered as a refinement parameter. Defaults to False.
    allow_atomic_params : Allow Atomic positions ('U', 'XU') to be considered as refinement paremeter. Defaults to False
    allow_angle : Allow acute angle refinement. Defaults to False
    allow_strain : Allow strain-related parameters to be considered for refinement. Defaults to False
    allow_broad :  Allow broadening-related parameters to be considered for refinement. Defaults to False
    force_ori : Force to use pref. orient. always. Defaults to False
    m : Parameter to optmize. Defaults to "rwp" for Rwp. For GoF use "gof"
    constraints : Any kind of GSAS-II acceptable refinement constraint dictionary. See GSASII docs for more info.
    """

    def __init__(self,
                 gpx_path: str,
                 working_folder: str,
                 random_state: int,
                 n_trials: int = 100,
                 n_jobs: int = -1,
                 allow_pref_orient: bool = True,
                 allow_atomic_params: bool = False,
                 allow_angle: bool = False,
                 allow_strain: bool = False,
                 allow_broad: bool = False,
                 force_ori: bool = True,
                 m: str = "rwp",
                 constraints = None):


        self.gpx_path = gpx_path
        self.wf = working_folder
        self.random_state = random_state
        self.project = G2sc.G2Project(gpxfile=self.gpx_path)
        self.gpx_name = os.path.basename(gpx_path).replace(".gpx", "")
        self.nphases = len(self.project.phases())
        self.phaselist = combine(
            [str(i) for i in range(self.nphases)])  # used for setting up which phases to refine prefered orientation
        self.ntrials = n_trials
        self.njobs = n_jobs
        self.study = optuna.study.Study
        self.gpx_best = None
        self.best_trial = None
        # get info about exp. pattern
        theta = self.project.histograms()[0].getdata("X")
        self.tmin = np.round(theta.min(), 2)
        self.tmax = np.round(theta.max(), 2)
        self.prefori = allow_pref_orient
        self.atompos = allow_atomic_params
        self.reflim = allow_angle
        self.strain = allow_strain
        self.broad = allow_broad
        self.opt_param = m
        self.force_ori = force_ori
        self.rwp_best = 999999
        self.best_trial_number = None
        self.trials_df = None
        self.lattices_df = []
        self.best_cell_dfs = []
        self.constraints = constraints


    def objective_sp(self, trial: optuna.trial) -> float:
        """
        Objective Function for Optuna
        for one-phase refinement case.

        Parameters
        ----------
        trial : optuna.trial object

        Returns
        -------
        Rwp : float

        """

        ### define search space ###

        #use iso
        use_iso = trial.suggest_categorical('Use Iso', [True, False])
        # max cycle
        max_cyc = trial.suggest_int('Max cycles', 1, 20 + 1)

        # Limits (acute angle)
        if self.reflim:
            limits_lb = trial.suggest_uniform('Limits lower bound', self.tmin, self.tmax - 20)
            limits_ub = trial.suggest_uniform('Limits upper bound', limits_lb + 20, self.tmax)
            limits_refine = trial.suggest_categorical('limits refine', [True, False])
            refdict0 = {'set': {'Limits': [limits_lb, limits_ub]}, 'refine': limits_refine}
        # Background
        else:
            refdict0 = None
        background_type = trial.suggest_categorical(
            'Background type', ['chebyschev',
                                'cosine',
                                'Q^2 power series',
                                'Q^-2 power series',
                                'lin interpolate',
                                'inv interpolate',
                                'log interpolate'])
        no_coeffs = trial.suggest_int('Number of coefficietns', 1, 15 + 1)  # [1, 16)
        background_refine = trial.suggest_categorical('Background refine', [True, False])
        refdict0bg_h = {
            'set': {
                'Background': {
                    'type': background_type,
                    'no. coeffs': no_coeffs,
                    'refine': background_refine
                }
            }
        }

        # Instrument and Sample Paramaters
        instrument_parameters1_refine = []
        for p in ['Zero']:
            if trial.suggest_categorical('Instrument_parameters refine %s' % (p), [True, False]):
                instrument_parameters1_refine.append(p)
        refdict1_h = {'set': {'Cell': True, 'Instrument Parameters': instrument_parameters1_refine}}

        # Strain and Broadening
        # Strain
        if self.strain:
            strain_params = []
            strain_ref = True
            bab = False
            for p in ['Babinet', 'Extinction', 'HStrain']:
                if trial.suggest_categorical('Strain paramaters refine %s' % (p), [True, False]):
                    strain_params.append(p)
                if 'Babinet' in strain_params:
                    bab_choice = trial.suggest_categorical('Babinet parameters', [['BabA'], ['BabU'], ['BabA', 'BabU']])
                    bab_dic = {"Babinet": bab_choice}
                    bab = True
            strain_dic = {p: True for p in strain_params if p != 'Babinet'}
            if bab:
                strain_dic.update(bab_dic)
            if len(strain_dic) == 0:
                strain_ref = False
                refdict_strain = None
            else:
                refdict_strain = {"set": strain_dic}
              #  print("Strain", refdict_strain, strain_ref)
        else:
            refdict_strain = None
            strain_ref = False

        # Broadening
        if self.broad:
            broad_ref = True
            broad_params = []
            for p in ['Size', 'Mustrain']:
                if trial.suggest_categorical('Broadening paramaters refine %s' % (p), [True, False]):
                    broad_params.append(p)
            broad_dic = {p: {"type": "isotropic", "refine": True} for p in broad_params}
            if len(broad_dic) == 0:
                broad_ref = False
                refdict_broad = None
            else:
                refdict_broad = {"set": broad_dic}
            # print("Broadening", refdict_broad, broad_ref)
        else:
            refdict_broad = None
            broad_ref = False

        sample_parameters1_refine = []
        for p in ['DisplaceX', 'DisplaceY', 'Scale']:
            if trial.suggest_categorical('Sample_parameters refine %s' % (p), [True, False]):
                sample_parameters1_refine.append(p)
        refdict1_h2 = {"set": {'Sample Parameters': sample_parameters1_refine}}

        # Instrument Parameters
        instrument_parameters2_refine = []
        for p in ['U', 'V', 'W', 'X', 'Y', 'SH/L']:
            if trial.suggest_categorical('Peakshape_parameters refine %s' % (p), [True, False]):
                instrument_parameters2_refine.append(p)
        refdict2_h = {'set': {'Instrument Parameters': instrument_parameters2_refine}}



        # Pref. Orient.
        if self.prefori:
            pref_ori = trial.suggest_categorical('Ref. Pref. Ori', [True, False])
            if pref_ori:
                refdict_ori = {"set": {"Pref.Ori.": pref_ori}}
            else:
                refdict_ori = None
        else:
            pref_ori = False
            refdict_ori = None

        # Atoms
        if self.atompos:
            refatoms = trial.suggest_categorical('Refine Atomic. Pos', [True, False])
            if refatoms:
                atom_params = trial.suggest_categorical('Atomic. Params', ['U', 'XU'])
                refdict3_h = {'set': {'Atoms': {'all': atom_params}}, "refine": refatoms}
            else:
                refdict3_h = None
        else:
            refatoms = False
            refdict3_h = None
            

        # Limits (wide angle)
        if self.reflim:
            refdict_fin_h = {'set': {'Limits': [self.tmin, self.tmax]}, 'refine': True}
        else:
            refdict_fin_h = None
        # Evaluate

        _params = {
            "max_cyc": [max_cyc, True],
            "use_iso": [use_iso, True],
            "method": [self.opt_param, True],
            "refdict0": [refdict0, self.reflim],
            "refdict0bg_h": [refdict0bg_h, True],
            "refdict1_h": [refdict1_h, True],
            "refdict_strain": [refdict_strain, strain_ref],
            "refdict_broad": [refdict_broad, broad_ref],
            "refdict1_h2": [refdict1_h2, True],
            "refdict_ori": [refdict_ori, pref_ori],
            "refdict2_h": [refdict2_h, True],
            "refidict3_h": [refdict3_h, refatoms],
            "refdict_fin_h": [refdict_fin_h, self.reflim]
        }
        refine_params_list = [v[0] for v in _params.values() if v[1]]


        def evaluate(trial_number, refine_params_list, q):
            """
            Evaluate refinement for given paramaters set
            Parameters
            ----------
            trial_number : int
                Optuna trial number

            refine_params_list : List[dict]
                A list of dictionaries with refinement directives

            q : multiprocessing queue object
                Queue object.

            """
            ERROR_PENALTY = 1e9
            #_Rwp = ERROR_PENALTY
            try:
                # set up a project trial
                project = OptunaTrial(trial_number=trial_number, gpx_path=self.gpx_path,
                                      working_folder=self.wf,
                                      random_state=self.random_state)
                if self.constraints is not None:
                    for cons in self.constraints:
                        project.gpx.add_constraint_raw(cons[0], cons[1])
                max_cyc = refine_params_list[0]  #
                iso = refine_params_list[1]
                m = refine_params_list[2]
                project.gpx.data['Controls']['data']['max cyc'] = max_cyc  # max cycles
                p = project.gpx.phases()[0]
                if iso:
                    for val in p.data['Atoms']:
                        val[9] = 'I'
                # Begin refinement
                for params in refine_params_list[3:]:
                    _Rwp = project.refine_and_calc_Rwp(params,m=m)
                # validate Uiso >= 0 for all phases
                u_iso_list = [atom.uiso for atom in p.atoms()]
                if min(u_iso_list) < 0:
                    # Uiso < 0
                    #   print("UISO ERROR")
                    _Rwp = ERROR_PENALTY
                # Check for large rwp. Meaning that refinement failed so give a penalty
                # try:
                #     if Rwp >= 100:
                #         print("Rwp greater than 100, setting up as {}".format(ERROR_PENALTY))
                #         Rwp = ERROR_PENALTY
                # except:
                #         Rwp = ERROR_PENALTY
                if _Rwp is None:
                    # print(refine_params_list)
                    _Rwp = ERROR_PENALTY
                q.put(_Rwp)
                # print("Trial number {}, Rwp: {}".format(trial_number, _Rwp))

            except Exception as e:
                # Refinement failed
                print(e, file=sys.stderr)
                q.put(ERROR_PENALTY)

        q = Queue()
        p = Process(target=evaluate, args=(trial.number, refine_params_list, q))
        p.start()
        # print(q)
        Rwp = q.get()
        # print(Rwp)
        if self.opt_param == "rwp":
            n = "Rwp"
        else:
            n = "GOF"
        if Rwp < self.rwp_best:
            self.rwp_best = Rwp
            self.best_trial_number = trial.number
            self.gpx_best = self.gpx_name + "_seed{0}_trial_{1}.gpx".format(self.random_state, self.best_trial_number)
        print("Trial: {} - {}: {:.3f}\n Current best {}: {:.3f} in trial number: {}".format(trial.number, n, Rwp, n, self.rwp_best, self.best_trial_number))

        p.terminate()
        p.join()
        return Rwp


    def objective_mp(self, trial: optuna.trial) -> float:
        """
        Objective Function for Optuna
        for multiple-phase refinement case.

        Parameters
        ----------
        trial : optuna.trial object

        Returns
        -------
        Rwp : float
        """

        # use iso
        use_iso = trial.suggest_categorical("Use Iso", [True, False])
        # Max cycles and pref. orientation & atom phases
        max_cyc = trial.suggest_int("Max cycles", 1, 20 + 1)

        # Limits (acute angle)
        if self.reflim:
            limits_lb = trial.suggest_uniform('Limits lower bound', self.tmin, self.tmax - 20)
            limits_ub = trial.suggest_uniform('Limits upper bound', limits_lb + 20, self.tmax)
            limits_refine = trial.suggest_categorical('limits refine', [True, False])
            refdict0 = {'set': {'Limits': [limits_lb, limits_ub]}, 'refine': limits_refine}
        else:
            refdict0 = {}

        
        # Background
        background_type = trial.suggest_categorical(
            'Background type', ['chebyschev',
                                'cosine',
                                'Q^2 power series',
                                'Q^-2 power series',
                                'lin interpolate',
                                'inv interpolate',
                                'log interpolate'])
        no_coeffs = trial.suggest_int('Number of coefficietns', 1, 15 + 1)  # [1, 16)
        background_refine = trial.suggest_categorical('Background refine', [True, False])

        # This is the base dictionary that is needed to do multiphase refinement in GSAS II (according to documentation/tutorials)
        refdict0bg_h = {
            'set': {
                'Background': {
                    'type': background_type,
                    'no. coeffs': no_coeffs,
                    'refine': background_refine
                },
                "Instrument Parameters": ["Zero"],
                "Scale": True,
                "Cell": True

            }
        }

        # Strain and Broadening
        # Strain
        if self.strain:
            strain_params = []
            strain_ref = True
            bab = False
            for p in ['Babinet', 'Extinction', 'HStrain']:
                if trial.suggest_categorical('Strain paramaters refine %s' % (p), [True, False]):
                    strain_params.append(p)
                if 'Babinet' in strain_params:
                    bab_choice = trial.suggest_categorical('Babinet parameters', [['BabA'], ['BabU'], ['BabA', 'BabU']])
                    bab_dic = {"Babinet": bab_choice}
                    bab = True
            strain_dic = {p: True for p in strain_params if p != 'Babinet'}
            if bab:
                strain_dic.update(bab_dic)
            if len(strain_dic) == 0:
                strain_ref = False
                refdict_strain = None
            else:
                phases_strain= trial.suggest_categorical('Strain.Phases', self.phaselist)
                pstr = [int(i) for i in phases_strain.split(",")]
                refdict_strain= {"set": strain_dic, "phases": pstr}
        #    print("Strain", refdict_strain, strain_ref)
        else:
            refdict_strain = None
            strain_ref = False

        # Broadening
        if self.broad:
            broad_ref = True
            broad_params = []
            for p in ['Size', 'Mustrain']:
                if trial.suggest_categorical('Broadening paramaters refine %s' % (p), [True, False]):
                    broad_params.append(p)
            broad_dic = {p: {"type": "isotropic", "refine": True} for p in broad_params}
            if len(broad_dic) == 0:
                broad_ref = False
                refdict_broad = None
            else:
                phases_broad = trial.suggest_categorical('Pref.Broad.Phases', self.phaselist)
                pbrd = [int(i) for i in phases_broad.split(",")]
                refdict_broad = {"set": broad_dic, "phases": pbrd}
        #    print("Broadening", refdict_broad, broad_ref)
        else:
            refdict_broad = None
            broad_ref = False
        # Sample and Instrument parameters

        sample_parameters1_refine = []
        for p in ['DisplaceX', 'DisplaceY', 'Scale']:
            if trial.suggest_categorical('Sample_parameters refine %s' % (p), [True, False]):
                sample_parameters1_refine.append(p)
        refdict1_h2 = {"set": {'Sample Parameters': sample_parameters1_refine}}

        # Pref. Orientation
        if self.prefori:
            if not self.force_ori:
                ori_refine = trial.suggest_categorical('Pref Ori Refine', [True, False])
                # ori_refine = True
            else:
                ori_refine = True
            if ori_refine:
                phases_ori = trial.suggest_categorical('Pref.Ori.Phases', self.phaselist)
                pori = [int(i) for i in phases_ori.split(",")]
                #     print(pori)
            else:
                pori = [0]
            refdict_ori = {"set": {"Pref.Ori.": ori_refine}, "phases": pori}
        else:
            ori_refine = False
            refdict_ori = {}

        # Instr. params
        instrument_parameters2_refine = []
        for p in ['U', 'V', 'W', 'X', 'Y', 'SH/L']:
            if trial.suggest_categorical('Peakshape_parameters refine %s' % (p), [True, False]):
                instrument_parameters2_refine.append(p)
        refdict2_h = {'set': {'Instrument Parameters': instrument_parameters2_refine}}

        # Atoms Positions and Thermal Paramaters
        if self.atompos:
            xu_refine = trial.suggest_categorical("XU_Refine", [True, False])
            if xu_refine:
                 xu_params = trial.suggest_categorical("XU_Opt", ["U", "XU"])
                 phases_xu = trial.suggest_categorical('XU. Phases', self.phaselist)
                 pxu = [int(i) for i in phases_xu.split(",")]
            else:
                 xu_params = "U"
                 pxu = [0]
            refdict3_h = {"set": {"Atoms": {"all": xu_params}}, "refine": xu_refine, "phases": pxu}
        else:
            xu_refine = False
            refdict3_h = {}
        # Wide angle
        # Limits (wide angle)
        refdict_fin_h = {'set': {'Limits': [self.tmin, self.tmax]}, 'refine': True}

        # Evaluate
        _params = {
            "max_cyc": [max_cyc, True],
            "use_iso": [use_iso, True],
            "method": [self.opt_param, True],
            "refdict0": [refdict0, True],
            "refdict0bg_h": [refdict0bg_h, True],
            "refdict_strain": [refdict_strain, strain_ref],
            "refdict_broad": [refdict_broad, broad_ref],
            "refdict1_h2": [refdict1_h2, True],
            "refdict_ori": [refdict_ori, ori_refine],
            "refdict2_h": [refdict2_h, True],
            "refidict3_h": [refdict3_h, xu_refine],
            "refdict_fin_h": [refdict_fin_h, True]
        }
        refine_params_list = [v[0] for v in _params.values() if v[1]]
        # refine_params_list = [max_cyc,
        #                       refdict0bg_h,
        #                       refdict1_h2,
        #                       refdict_ori,
        #                       refdict2_h,
        #                       refdict3_h,
        #                       # refdict_fin_h
        #                       ]

        # print(refine_params_list)
        def evaluate(trial_number, refine_params_list, q):
            """
            Evaluate refinement for given paramaters set
            Parameters
            ----------
            trial_number : int
                Optuna trial number

            refine_params_list : List[dict]
                A list of dictionaries with refinement directives

            q : multiprocessing queue object
                Queue object.
            """

            ERROR_PENALTY = 1e9
            try:
                project = OptunaTrial(trial_number=trial_number, gpx_path=self.gpx_path,
                                      working_folder=self.wf,
                                      random_state=self.random_state)

                if self.constraints is not None:
                    for cons in self.constraints:
                        project.gpx.add_constraint_raw(cons[0], cons[1])

                max_cyc = refine_params_list[0]
                iso = refine_params_list[1]
                m = refine_params_list[2]
                project.gpx.data['Controls']['data']['max cyc'] = max_cyc  # max cycles
                phases = project.gpx.phases()
                # Begin refinement
                # setup isotropic paramaeter. Maybe we have to also allow it to chose which p
                if iso:
                    for phase in phases:
                        for val in phase.data['Atoms']:
                            val[9] = 'I'
                for params in refine_params_list[3:]:
                    Rwp = project.refine_and_calc_Rwp(params, m=m)
                # validate Uiso >= 0 for all phases
                for p in phases:
                    u_iso_list = [atom.uiso for atom in p.atoms()]
                    if min(u_iso_list) < 0:
                        # Uiso < 0
                        Rwp = ERROR_PENALTY

                # Check if refinement fineshed.....
                wtfrac = get_wtfract(project.trial_path.replace(".gpx", ".lst"), self.nphases)
                if min(wtfrac) < 0:  # Some phase has negative %
                    Rwp = ERROR_PENALTY
                q.put(Rwp)

            # Failed refinement
            except Exception as e:
                # Refinement failed
                print(e, file=sys.stderr)
                q.put(ERROR_PENALTY)

        q = Queue()
        p = Process(target=evaluate, args=(trial.number, refine_params_list, q))
        p.start()
        # print(q)
        Rwp = q.get()
        if self.opt_param == "rwp":
            n = "Rwp"
        else:
            n = "GOF"
        if Rwp < self.rwp_best:
            self.rwp_best = Rwp
            self.best_trial_number = trial.number
            self.gpx_best = self.gpx_name + "_seed{0}_trial_{1}.gpx".format(self.random_state, self.best_trial_number)
        print("Trial: {} - {}: {:.3f}\n Current best {}: {:.3f} in trial number: {}".format(trial.number, n, Rwp, n,
                                                                                            self.rwp_best,
                                                                                            self.best_trial_number))
        # print(Rwp)
        p.terminate()
        p.join()
        return Rwp

    def create_study(self, n_startup: int = 20) -> BlackboxOptimizer:
        """
        Create an Optuna study

        Parameters
        ----------
        n_startup : int
            Number of starting trials for Optuna run

        Returns
        -------
        self
            Returns self.
        """
        # Create Optuna study
        self.study = optuna.create_study(study_name=self.gpx_name + '_seed%s' % (self.random_state),
                                         sampler=optuna.samplers.TPESampler(n_startup_trials=n_startup,
                                                                            seed=self.random_state))
        return self

    def optimize(self) -> None:
        """
        Begin Optimization

        """
        if self.nphases == 1:
            print("One phase only found. Running optimizer for a single phase")
            self.study.optimize(self.objective_sp, n_trials=self.ntrials, n_jobs=self.njobs)
        else:
            print("More than one phase found. Running optimizer for more than one phase")
            self.study.optimize(self.objective_mp, n_trials=self.ntrials, n_jobs=self.njobs, show_progress_bar = True)
        self.best_trial = self.study.best_trial
        self.gpx_best = self.gpx_name + "_seed{0}_trial_{1}.gpx".format(self.random_state, self.best_trial.number)
        self.trials_df = self.study.trials_dataframe()
        seed_name = "seed" + str(self.random_state)
        gpxs_folder = os.path.join(self.wf, "optuna")
        to_read = [d for d in os.listdir(gpxs_folder) if (d.endswith('.lst') and seed_name in d)]
        to_read = sorted(to_read, key=lambda e: eval(e.split(".lst")[0].split("_")[-1]))

        for k in range(self.nphases):
            data_ = []
            for i, p in enumerate(to_read):
                lst_path = os.path.join(gpxs_folder, p)
                result = get_cell_density(lst_path, self.nphases)[k]
                result['rwp'] = self.trials_df['value'].iat[i]
                result['trial'] = i
                data_.append(result)
            cells = pd.DataFrame(data_)
            self.lattices_df.append(cells)
        self.best_cell_dfs = get_cell_density(os.path.join(self.wf, "optuna", self.gpx_best.replace(".gpx", ".lst")
                                                           ), self.nphases)

    def plot_best(self, save: bool = False, path: str = None, engine = "matplotlib") -> Union[Figure, None] :
        """
        Plot best refinement from the optimization process.

        Parameters
        ----------
        save : bool, default: True
            If this is True, the figure will be saved to `path`

        path : str, default: None
            In case that save is True, save the figure to the given path.

        engine : str, default = "matplotlib"
            Which plotting engine to use.
            For static plots, use "matplotlib" (default)
            For interactive HTML plots using plotly, use "plotly"

        Returns
        -------
        A Figure or None
            Returns a plotly figure object if plotly is chosen. Otherwise returns None.
        """
        trial_path = os.path.join(self.wf, "optuna", self.gpx_best)
        if not save:
            fig = plot_gpx(trial_path, save=save, engine=engine)
        else:
            fig = plot_gpx(trial_path, save=save, outfolder=path, engine=engine)
        return fig

    def set_optuna_verbosity(self, verbosity: int = 0) -> BlackboxOptimizer:
        """
        Set Optuna verbosity.

        Parameters
        ----------
        verbosity : int, default: 0
            Sets Optuna verbosity. Defaults to 0 (logging.CRITICAL)
            For more information about the options please check Optuna documentation.

        Returns
        -------
        self
        """
        optuna.logging.set_verbosity(verbosity)
        return self

    def rwp_plot(self, limits: Tuple[int, int] = (5, 50)) -> None:
        """
        Plot the evolution of min Rwp in function of trials
        Parameters
        ----------
        limits : Tuple[int, int]
            A tuple with (min, max) values to plot (y-scale) of Rwp.

        Returns
        -------
        None

        """
        _min, _max = limits
        df = self.study.trials_dataframe()
        # df = df[df['value'] <= _max]
        df.reset_index(drop=True, inplace=True)
        firstmin = df['value'].iat[0]
        minvals = [firstmin]
        currmin = firstmin
        for rwp in df['value'][1:]:
            if rwp < currmin:
                currmin = rwp
            minvals.append(currmin)
        fig, ax = plt.subplots(figsize=(12, 8))
        c = Mplibcfg()
        plt.plot(minvals, **c.pargs)
        plt.xlabel("Trials", **c.largs)
        plt.ylabel("Rwp (%)", **c.largs)
        plt.tick_params(**c.targs)
        plt.ylim([_min, _max])
        plt.show()

    def plot_lattice_trials(self, lattice: str = "a", engine: str = "matplotlib",
                            rwp_cutoff: float = 25, phase_number: int = 0) -> None:
        """
        Plot all the lattice parameters obtained by the optimization process.

        Parameters
        ----------
        lattice : str, default: "a"
            Defines which lattice constants to plot. Options: "a", "b" or "c"

        engine : str, default: "matplotlib"
            Set the plotting engine. "Matplotlib" for matplotlib, "plotly" for plotly

        rwp_cutoff : float, default = 25
            The min rwp cutoff to show. For example with rwp = 25, only crystal structures with Rwp(%) < 25 will be plotted.

        phase_number : int, default = 0
            In case of more than one phase refinement, defines which phase to plot. Defaults to 0.

        Returns
        ----------
        None
        """

        from Xerus.settings.mplibcfg import Mplibcfg
        import seaborn as sns

        cells = self.lattices_df[phase_number]
        cells = cells[cells.rwp <= rwp_cutoff]
        cells.reset_index(drop=True, inplace=True)
        best = cells[cells.trial == self.best_trial_number]
        x = "trial"
        y = lattice
        y_error = "err_" + lattice

        if engine == "plotly":
            fig = px.scatter(data_frame=cells,
                             x=x,
                             y=y,
                             error_y=y_error,
                             color="rwp")
            # Add a point for the best trial
            fig.add_scatter(x=best["trial"],
                            y=best[y],
                            error_y=dict(array=np.array(best[y_error])),
                            name="Best trial",
                            customdata=[best['rwp']],
                            hovertemplate="Best trial<br>Rwp: %{customdata[0]}<br>Latt cons:%{y}<br><extra></extra>"
                            )
            fig.update_layout(legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
            ))

            fig.show()
        if engine == "matplotlib":
            ## style configurations
            sns.set()
            sns.set_style("darkgrid", {"patch.edgecolor": "black", "axes.edgecolor": "black"})
            sns.set_context("talk")
            c = Mplibcfg()
            c.largs['fontweight'] = 'normal'
            c.pargs['markersize'] = 12
            c.largs['fontsize'] = 24
            c.lgdnargs['fontsize'] = 20
            fig, ax = plt.subplots(figsize=(6, 6))
            plt.errorbar(cells[x], cells[y], yerr=cells[y_error], fmt="o", capsize=5, **c.pargs, label="_nolegend_")
            plt.errorbar(best[x], best[y], yerr=best[y_error], fmt='o', capsize=5, **c.pargs, label="Best trial")
            plt.legend(**c.lgdnargs)
            plt.tick_params(**c.targs)
            plt.xlabel("$trial$ (number)", **c.largs)
            plt.ylabel(y + " ($\AA$)", **c.largs)
