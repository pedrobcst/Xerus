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
from typing import List, Union, Tuple
import warnings
import os
import sys
from pathlib import Path
import shutil
project_path = str(Path(os.path.dirname(os.path.realpath(__file__))).parent) + os.sep  # so convoluted..
if project_path not in sys.path:
    sys.path.append(project_path)
import pandas as pd
from Xerus.settings.settings import GSAS2_BIN
from Xerus.settings.settings import INSTR_PARAMS
from Xerus.utils.tools import create_folder
from Xerus.optimizer.bbo import BlackboxOptimizer
from Xerus.utils.cifutils import get_ciflist, get_spgname
from .gsas2riet import run_gsas_mp2, quick_gsas
from .gsas2viz import plot_gpx
from .gsas2utils import get_cell_density, get_plot_data, parse_reflections_gpx

sys.path.append(GSAS2_BIN)
warnings.filterwarnings('ignore')
import GSASIIscriptable as G2sc  # Ignore Warning


class GSASRefiner:
    """
    A wrapper for ease of use BlackBox Optimization method.
    Automatically handles creation of gpx files from cifs and experimental data.


    Parameters
    ----------
    expdata : str
        Path to experimental data file. This file has to be in the any of GSASII supported formats.
    cifs : List[str] or str
        A list of path to CIF Files OR a path to a folder containing the CIF files.
        In case of one phase, provide a one-element list, ie:['path/to/MyCIF.cif']
    working_folder : str
        Path to folder to save the analysis and results
    gpxname : str, default: "auto"
        Base-gpx project name. Defaults to "auto": The gpx project filename will just be the expdata with ".gpx" instead.
    """

    def __init__(self, expdata: str, cifs: Union[List[str], str], working_folder: str, gpxname: str = "auto"):
        create_folder(working_folder)
        self.expdata = expdata
        self.wf = working_folder
        self.name = os.path.basename(self.expdata).split(".")[0]

        # Parse cif folder
        if type(cifs) == str:
            ciflist = get_ciflist(cifs)
            self.phasenames = [c['name'] + '-' + c['spacegroup'].replace("/", "_") for c in ciflist]
            self.ciflist = [os.path.join(cifs, c['filename']) for c in ciflist]
        else:
            self.ciflist = cifs
            self.phasenames = [get_spgname(cif) for cif in
                               self.ciflist]

        # In case of no gpx name is given
        if gpxname == "auto":
            self.gpxname = os.path.basename(self.expdata).split(".")[0] + ".gpx"
        else:
            self.gpxname = gpxname

        # Save gpx path
        self.gpxpath = os.path.join(self.wf, self.gpxname)
        # initialize GSAS project
        self.gpx = G2sc.G2Project(newgpx=self.gpxpath)
        ## add patterns & phases
        hist1 = self.gpx.add_powder_histogram(datafile=self.expdata,
                                              iparams=INSTR_PARAMS)
        phases = []
        hist1.data['Instrument Parameters'][0]['I(L2)/I(L1)'] = [0.5, 0.5, 0]
        # add cifs & uiso
        for i, cif in enumerate(self.ciflist):
            phases.append(self.gpx.add_phase(cif, phasename=self.phasenames[i], histograms=[hist1]))
            # # Set to use iso
            # for val in phases[i].data['Atoms']:
            #     val[9] = 'I'

        self._quickgpx = None
        self._quickrwp = None
        self._quickwt = None
        self._qglat = None

        @property
        def quick_gpx(self):
            if self._quickgpx is not None:
                return self._quickgpx
            else:
                return "No quick analysis has been ran. Please run quick_refinement option."

        @quick_gpx.setter
        def quick_gpx(self, value):
            self._quickgpx = value

        @property
        def quick_rwp(self):
            if self._quickrwp is not None:
                return self._quickrwp
            else:
                return "No quick analysis has been ran. Please run quick_refinement option."

        @quick_rwp.setter
        def quick_rwp(self, value):
            self._quickrwp = value

        @property
        def quick_wt(self):
            if self._quickwt is not None:
                return self._quickwt
            else:
                if len(self.ciflist) > 1:
                    return "No quick analysis has been ran. Please run quick_refinement option."
                else:
                    return 100

        @quick_wt.setter
        def quick_wt(self, value):
            self._quickwt = value

    def fit_optimize(self, random_state: int = 42, n_jobs: int = -1, n_trials: int = 200,
                     n_startup: int = 20,
                     allow_pref_orient: bool = True,
                     allow_atomic_params: bool = False,
                     allow_broad: bool = False,
                     allow_strain: bool = False,
                     allow_angle: bool = False,
                     force_ori: bool = False,
                     verbose: str = None,
                     param: str = "rwp",
                     constraints = None,
                     ) -> GSASRefiner:
        """
        Run a blackbox optimization in the current loaded data and crystal structures.
        Please a take a while to note all the possible parameters for optimization.
        Set the options according to your data quality etc.

        Parameters
        ----------
        random_state : int, default: 42
            Sets the Optuna randomstate seed.
        n_jobs : int, default: -1
            Set optuna number of jobs. Note: this might be Deprecated in newest versions of Optuna
        n_trials : int, default: 200
            Set the number of trials to be done
        n_startup : int, default: 20
            Set the number of startup trials of the study.
        allow_pref_orient : bool, default: True
            Allow for pref. orientation. to be included in the refinement parameter search
        allow_atomic_params : bool, default: False
            Allow for atomic parameters (U, XU) to be considered as a refinement parameter
        allow_broad : bool, default: False
            Allow broadening terms to be considered for refinement
        allow_strain : bool, default: False
            Allow for strain terms to be considered for refinement
        allow_angle : bool, default: False
            Allow for acute angle refinement or not
        force_ori : bool, default: False
            If allow_pref is True, this will force the optimizer to always consider pref. orientation (in test)
        verbose : str, default: None
            Sets Optuna verbosity:
            Use "silent" for verbosity 0, otherwise it will be set to 1.
        param : str, default: "rwp"
            Objective goal. Defaults to "rwp". Options:
            "rwp": Rwp
            "gof": Goodness of Fit (GoF)
        constraints : dict-like, default: None
            A GSASII-like constraint dictionary. For details on how to setup any kind of constraints in the refinement
            (ie, not allowing negative occupancies), see GSAS II documentantion.

        Returns
        -------
        GSASRefiner (self)
        """

        self.optim = BlackboxOptimizer(
            gpx_path=self.gpxpath,
            working_folder=self.wf,
            random_state=random_state,
            n_trials=n_trials,
            n_jobs=n_jobs,
            allow_pref_orient=allow_pref_orient,
            allow_atomic_params=allow_atomic_params,
            allow_angle=allow_angle,
            allow_strain=allow_strain,
            allow_broad=allow_broad,
            m=param,
            force_ori=force_ori,
            constraints=constraints)
        if verbose == "silent":  # Make optuna less verbose
            self.optim.set_optuna_verbosity(0)
        else:
            self.optim.set_optuna_verbosity(1)

        self.optim.create_study(n_startup=n_startup).optimize()
        self.best = self.optim.gpx_best
        self.lattice_best = self.optim.best_cell_dfs
        return self

    def plot_result(self, save: bool = False, path: str = None, engine="matplotlib") -> None:
        """
        Plot the result with lowest obtained Rwp.

        Parameters
        ----------
        save : bool, default: False
            True to save to a folder defined by `path`
            False otherwise
        path : str, default: None
            In case that save is True, defines the folder for saving the plot
        engine : str, default: "matplotlib"
            Defines the backend for plotting
            "matplotlib": Uses matplotlib backend
            "plotly": Uses plotly backend

        Returns
        -------
        None
        """
        # Check for engine
        if engine not in ["matplotlib", "plotly"]:
            raise NotImplementedError("Plotting backend not supported.")

        if self.optim:
            fig = self.optim.plot_best(save=save, path=path, engine=engine)
            return fig

    def plot_rwp(self, limits: Tuple[int, int] = (5, 50)) -> None:
        """
        Plot the evolution of Rwp as a function of trials.

        Parameters
        ----------
        limits : tuple (int, int)
            Specificy the rwp min and max limits to show in the evolution plot

        Returns
        -------
        None
        """

        if self.optim:
            self.optim.rwp_plot(limits=limits)

    def quick_refinement(self) -> GSASRefiner:
        """
        Performs a quick refinement.
        This is used to check if the loaded CIFS and experimental data are working correctly and were not chosen wrongly.

        Returns
        -------
        self
        """

        n_phases = len(self.ciflist)
        # multiple phase
        if n_phases > 0:
            rwp, gpx, wt = run_gsas_mp2(powder_data=self.expdata,
                                        cifs=self.ciflist,
                                        phasenames=self.phasenames,
                                        outfolder=self.wf)


            self.quickgpx = gpx
            self.quickrwp = rwp
            self.quickwt = wt
        plot_gpx(self.quickgpx, save=False)

        return self

    def plot_lattice_trials(self, lattice: str = "a", engine: str = "matplotlib",
                            rwp_cutoff: float = 25, phase_number: int = 0) -> None:
        """
        Plot the obtained crystal structures (lattice constants) obtained from all the runs.

        Parameters
        ----------
        lattice : str, default: "a"
            Define which lattice constant to plot. ("a", "b" or "c")
        engine : str, default: "matplotlib"
            Define which plotting backend to use.
            "matplotlib": Uses matplotlib backend
            "plotly": Uses plotly backend
        rwp_cutoff : float, default: 25.0
            Defines the minimum "rwp" to be considered to plot. If all crystal structures are desired, increase this
            number to a large number
        phase_number : int, default: 0
            Defines the "phase" number to be plotted, in case of multiple phase refinement.

        Returns
        -------
        None
        """
        # Check for correct usage of lattice parameter
        if lattice not in ["a", "b", "c"]:
            raise ValueError("Invalid lattice constant parameter.")

        # Check for implemented plotting backends
        if engine not in ["matplotlib", "plotly"]:
            raise NotImplementedError("Invalid plotting backend.")

        self.optim.plot_lattice_trials(lattice=lattice,
                                       engine=engine,
                                       rwp_cutoff=rwp_cutoff,
                                       phase_number=phase_number)

    def export_results(self, outfolder="export") -> None:
        """
        Exported the best (minimum Rwp/GoF) results:
            - Crystal structure information (a,b,c,V,density, phase fractions for each phase, rwp plus erros)
                -> Filename: `self.name` + "_structure".csv
            - Yobs, Ycalc, Ydiff (Rietveld refinement details)
                -> Filename: `self.name` + ".csv"
            - GSASII .lst file
            - Bragg reflections of each phase in reflections folder.
        Parameters
        ----------
        outfolder : str, default: "export"
            Defines the folder where the data will be exported to.

        Returns
        -------
        None

        """
        path = os.path.join(self.wf, outfolder)
        create_folder(path)
        lst_file = os.path.join(self.wf, 'optuna', self.best.replace(".gpx", ".lst"))
        _, theta, yobs, ycalc, ydiff, rwp, wt = get_plot_data(os.path.join(self.wf,'optuna',self.best))
        # make lattice into df
        crystal_df = pd.DataFrame(data=self.lattice_best)
        # export df
        crystal_df['rwp'] = rwp
        crystal_df['wt'] = crystal_df['name'].map(wt)
        # plot df
        plot_data = pd.DataFrame(data=zip(theta, yobs, ycalc, ydiff), columns=['theta', 'yobs', 'ycalc', 'ydiff'])
        # reflections of each phase
        reflections = parse_reflections_gpx(os.path.join(self.wf,'optuna',self.best))

        # file names
        main_name = f"{self.name}.csv"
        crystal_name = f"{self.name}_structure.csv"
        # save to folder
        crystal_df.to_csv(os.path.join(path, crystal_name), index=False)
        plot_data.to_csv(os.path.join(path, main_name), index=False)

        # copy lst file
        lst_name = os.path.basename(lst_file) # Get basename from path
        shutil.copy(lst_file, os.path.join(path, lst_name))
        # make a reflection folder
        refl_path = os.path.join(path, "reflections")
        create_folder(refl_path)
        # save each
        for refl in reflections:
            name = f"{refl['Phase'].iat[0]}_reflections.csv"
            refl.to_csv(os.path.join(refl_path, name), index = False)
        return self

