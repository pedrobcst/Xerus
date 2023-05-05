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

__version__ = '1.1b'

import os
from multiprocessing import Pool
import concurrent.futures
from typing import Any, List, Union

import numpy as np
import pandas as pd
import plotly.express as px
from pymatgen.core import Composition

from Xerus.db.localdb import LocalDB
from Xerus.engine.gsas2riet import refine_comb, simulate_spectra
from Xerus.engine.gsas2utils import make_gpx
from Xerus.engine.gsas2viz import plot_all_gpx, plot_gpx
from Xerus.engine.gsas2wrap import GSASRefiner
from Xerus.readers.datareader import DataReader
from Xerus.settings.settings import GSAS2_BIN, INSTR_PARAMS
from Xerus.similarity.pattern_removal import (combine_pos,
                                              run_correlation_analysis,
                                              run_correlation_analysis_riet)
from Xerus.similarity.visualization import make_plot_all, make_plot_step
from Xerus.utils.cifutils import make_system
from Xerus.utils.preprocessing import remove_baseline, standarize_intensity
from Xerus.utils.tools import (create_folder, group_data, load_json,
                               make_offset, normalize_formula)
from Xerus.utils.tools import plotly_add as to_add
from Xerus.utils.tools import save_json


class XRay:
    """
    Main handler for detecting phases of experimental XRD patterns

    Parameters
    ----------
    name : Sample / Project Name (str)
    working_folder : Folder to save analysis results (str)
    elements : A list of possible elements in the pattern to be analyze. Ie: For HoB2 elements = ['Ho', 'B'].
    exp_data_file : Path to experimental data to be analyzed
    data_fmt : Experimental data format. Note that it has to be supported by the DataExtensions/DataReader
    maxsys : Limit search to a maximum system size. Ie, if maxsys=3 and len(elements) = 4, it searchs only up to ternary systems
    max_oxy : Similar to above but for oxides restriction. For example, if max_oxy = 2, and element list is ['Ho', 'B', 'O'],
    the combinations to be search for oxygen containing systems will be only up to Ho-O and B-O, Ho-B-O will not be searched.
    remove_background: bool, default: True
        Argument to remove background (baseline) of XRD pattern. Defautls to True.
        If set to True, it will remove using the baseline estimation algorithm provided by Peakutils (polynomial)
    poly_degree: int, default: 8
        Polynomial degree to pass to peakutils.baseline for background estimation.
        Only relevant IF remove_background is True
    standarize_int: bool, default: True
        To standarize the intensity of the XRD-pattern by diving by Imax or not.
    use_preprocessed: bool, default: True
        If set to True, for GSASII refinements it will use the intensity-standarized data.
        (Note: It will only use intensity standarized data for refinement, NOT background removed data.)

    """

    def __init__(
        self,
        name: str,
        working_folder: str,
        elements: List[str],
        exp_data_file: str,
        data_fmt: str = "auto",
        maxsys: Union[Any, float] = None,
        max_oxy: int = 2,
        remove_background: bool = False,
        poly_degree: int = 8,
        standarize_int: bool = True,
        use_preprocessed: bool = True,
    ):

        if data_fmt == "auto":
            data_fmt = exp_data_file.split(".")[-1].lower()
            print(f"No datafmt passed. Assuming its {data_fmt}")
        self.exp_data, self.tmin, self.tmax, self.step = DataReader().read_data(
            exp_data_file, data_fmt
        )  # Read up the data
        self.standarize = standarize_int
        self.datafmt = data_fmt
        self.rb = remove_background
        self.pd = poly_degree
        self.exp_data_file = exp_data_file
        self.name = name
        self.working_folder = working_folder  # Set up working folder
        create_folder(working_folder)  # Create folder if doest not exist.
        self.elements = elements  # Set up elements
        self.instr_params = INSTR_PARAMS  # Instr params from your xrday machines
        self.filename = os.path.basename(exp_data_file)
        self.cif_all = None
        self.cif_info = None
        self.cif_notsim = None
        self.cif_notran = None
        self.simulated_gsas2 = []
        self.reflections_gsas2 = []
        self._notrun = []
        self.cif_info = None
        self.chosen_ids = []
        self.modified_df = None
        self.max_oxy = max_oxy
        self._gpx_best = None
        self._cif_best = None
        self._rwp_best = None
        self._results = None
        self.optimizer = None
        self.optgpx = None
        self.optlat = None
        self.optrwp = None
        self.use_preprocessed = use_preprocessed

        # Check up Preprocessing
        if standarize_int:
            print("Standarizing intensity to [0,1]..")
            self.exp_data = standarize_intensity(self.exp_data)
        if use_preprocessed:
            preprocess_name = f"{self.filename.split('.')[0]}_preprocessed.csv"
            self._preprocess_path = os.path.join(self.working_folder, preprocess_name)
            export_data = self.exp_data.copy()
            export_data.drop(["filename"], axis=1, inplace=True)
            export_data.to_csv(self._preprocess_path, index=False, header=False)
            self._old_data = self.exp_data_file
            self._old_fmt = self.datafmt
            print(f"Exported new datafile at {self._preprocess_path}")

        if remove_background:
            print(f"Removing background using polynomial degree: {poly_degree}")
            self.exp_data = remove_baseline(
                self.exp_data, poly_degree=poly_degree, plot=True
            )

        @property
        def rwp(self):
            if self._rwp_best is not None:
                return self._rwp_best
            else:
                return "Please run analyze() function first"

        @rwp.setter
        def rwp(self, new_value):
            self._rwp_best = new_value

        @property
        def best_gpx(self):
            if self._gpx_best is not None:
                return self._gpx_best
            else:
                return "Please run analyze() function first"

        @best_gpx.setter
        def best_gpx(self, new_value):
            self._gpx_best = new_value

        @property
        def cifs_best(self):
            if self._cif_best is not None:
                return self._cif_best
            else:
                return "Please run analyze() function first"

        @cifs_best.setter
        def cifs_best(self, new_value):
            self._cif_best = new_value

        @property
        def results(self):
            if self._results is not None:
                return self._results
            else:
                return "Please run analyze() first"

        @results.setter
        def results(self, new_value):
            self._results = new_value

        if maxsys is None and elements is not None:
            self.maxsys = len(elements)
        else:
            self.maxsys = maxsys

        self.corrinfo = None  ## :)

    def get_cifs(
        self,
        ignore_provider: List[str] = None,
        ignore_comb: List[str] = None,
        ignore_ids: List[str] = None,
    ) -> XRay:
        """
        Get cifs from MongoDB and write to working folder.
        If a certain combination of elements CIF is not present, it will automatically download
        Parameters
        ----------

        ignore_provider : Ignores a certain list of providers. In case of one, needs a list with one element, eg ["AFLOW"]
        ignore_comb : List of possible combinations to be ignored.
        ignore_ids: List of possible unique IDs to be ignored.

        Returns
        -------
        self
        """
        cif_meta, cif_notran, cif_notsim = LocalDB().get_cifs_and_write(
            element_list=self.elements,
            outfolder=self.working_folder,
            maxn=self.maxsys,
            max_oxy=self.max_oxy,
            name = self.name
        )
        self.cif_info = cif_meta
        if ignore_provider is not None:
            self.cif_info = self.cif_info[~self.cif_info.provider.isin(ignore_provider)]
            self.cif_info.reset_index(drop=True, inplace=True)
        if ignore_comb is not None:
            correct = []
            for comb in ignore_comb:
                _correct = "".join(comb.split("-"))
                correct.append(make_system(Composition(_correct).formula))
            self.cif_info = self.cif_info[~self.cif_info.system_type.isin(correct)]
            self.cif_info.reset_index(drop=True, inplace=True)
        if ignore_ids is not None:
            self.cif_info = self.cif_info[~self.cif_info.id.isin(ignore_ids)]
            self.cif_info.reset_index(drop=True, inplace=True)

        self.cif_notran = cif_notran
        self.cif_notsim = cif_notsim
        folder = self.working_folder
        name_out_cif = os.path.join(
            folder, os.path.normpath(folder).replace(os.sep, "_") + "_all_cifs.csv"
        )
        name_out_notran = os.path.join(
            folder, os.path.normpath(folder).replace(os.sep, "_") + "_not_used_riet.csv"
        )
        name_out_notsim = os.path.join(
            folder, os.path.normpath(folder).replace(os.sep, "_") + "_not_used_sim.csv"
        )

        self.cif_info.to_csv(name_out_cif, index=False)
        self.cif_notran.to_csv(name_out_notran, index=False)
        self.cif_notsim.to_csv(name_out_notsim, index=False)
        return self

    def simulate_all(self, n_jobs: int = -1, timeout_value: int = 10):
        """
        Parallel simulation of XRD patterns using Modin (Ray backend)
        Parameters
        ----------
        n_jobs: int, default: -1
            How many workers to use when starting multiprocessing Pool. If none is given, it will use all
        Returns
        -------
        self
        """
        if n_jobs > os.cpu_count():
            n_jobs = os.cpu_count()
        elif n_jobs == -1:
            n_jobs = None

        tmin = self.tmin
        tmax = self.tmax
        step = self.step
        ciflist = self.cif_info.copy()
        self.cif_all = self.cif_info.copy()  # keep a copy of all cifs
        working_folder = self.working_folder
        instr_params = self.instr_params
        print("Simulating {} patterns".format(len(ciflist)))
        args = [
            [f] + [tmin, tmax, step, working_folder, instr_params]
            for f in ciflist.full_path
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
            future_to_args = {executor.submit(simulate_spectra, *arg): arg for arg in args}
            results=[]
            for future in concurrent.futures.as_completed(future_to_args):
                args = future_to_args[future]
                try:
                    result = future.result(timeout=timeout_value)
                    results.append(result)
                except concurrent.futures.TimeoutError:
                    print(f'Task for file {args[0]} timed out after {timeout_value} seconds')
                    results.append(None)
                except Exception as e:
                    print(f'Task for file {args[0]} raised an exception: {e}')
                    results.append(None)
        
        paths = [None] * len(ciflist)
        for i, result in enumerate(results):
            if result is not None:
                paths[i] = result 
        print("Done. Cleaning up GSASII files.")
        # Clean up
        files = [
            file
            for file in os.listdir(os.getcwd())
            if file.endswith(".gpx") or file.endswith(".lst")
        ]
        for file in files:
            # print(f"Cleaning up {file}")
            os.remove(file)
        ciflist["simulated_files"] = [r[0] if r is not None else None for r in paths]
        ciflist["simulated_reflects"] = [r[1] if r is not None else None for r in paths]
        ciflist["sm_ran"] = [r[2] if r is not None else False for r in paths]
        ciflist = ciflist[ciflist["sm_ran"]]
        ciflist.drop(["sm_ran"], axis=1, inplace=True)
        ciflist.reset_index(drop=True, inplace=True)
        ciflist = pd.DataFrame(
            ciflist.to_records(index=False)
        )  # convert back to a pandas df
        self.cif_info = ciflist.copy()
        return self

    def _read_simulations(self) -> XRay:
        """
        Helper function for reading a simulations based on simulated (ran) information
        Returns
        -------
        self
        """
        reflections = []
        simulations = []
        for sim, ref, spg, spgnum, cs, name, provider, mid, fname, cij in zip(
            self.cif_info["simulated_files"],
            self.cif_info["simulated_reflects"],
            self.cif_info["spacegroup"],
            self.cif_info["spacegroup_number"],
            self.cif_info["crystal_system"],
            self.cif_info["name"],
            self.cif_info["provider"],
            self.cif_info["id"],
            self.cif_info["filename"],
            self.cif_info["Cij"],
        ):
            if spg != "None":
                sim_data = pd.read_csv(sim)
                ref_data = pd.read_csv(ref)
                sim_data["name"] = name
                sim_data["provider"] = provider
                sim_data["mid"] = mid
                sim_data["filename"] = fname
                sim_data["spacegroup"] = spg
                sim_data["spacegroup_number"] = spgnum
                sim_data["crystal_system"] = cs
                sim_data["Cij"] = cij
                reflections.append(ref_data)
                simulations.append((sim_data, fname))
        self.simulated_gsas2 = simulations
        self.reflections_gsas2 = reflections
        return self

    def select_cifs(
        self,
        cif_info: Union[pd.DataFrame, str] = "auto",
        save: bool = True,
        normalize: bool = True,
        by_systemtype: bool = True,
    ) -> pd.DataFrame:
        """
        Filter CIFs by correlation plus either stoichiometry + spacegroup or by system type + spacegroup.
        This method is done so, we avoid using alot of patterns of the same structure when searching for phases.
        Try to returns highest correlated pattern with experimental data.


        Parameters
        ----------
        cif_info : pd.DataFrame
            A Pandas dataframe containing query data from database
        save : bool, default: True
            Override cif info. Not used.
        normalize : bool, default: True
            Attempts to "normalize" compositions for groping (stoich+spacegroup) method.
            This case is to try to group off-stoichiometric compositions with same spacegroup that should be the same.
            Not so effective.
        by_systemtype : bool, default: True
            Instead of using stoichiometry, uses "system_type", ie: HoB2, HoB4 etc => Ho-B. Then


        Notes
        ------
        This method is an ATTEMPT at trying to reduce the number of identical structures that comes from different sources.
        Altough it can filter quite well, it fails for 'non-stoichiometric' cases
        A method for grouping identical crystal structures is still required.

        Returns
        -------
        pd.DataFrame
            Returns the filtered dataframe.

        """
        # This needs to be rewritten someday.
        if type(cif_info) == str:
            if cif_info == "auto":
                cif_info = self.cif_info
                if self.cif_info is None:
                    raise TypeError("Report file or no patterns have been simulated.")

        if not by_systemtype:
            if not normalize:
                dfs = group_data(cif_info, column_name=["name", "spacegroup_number"])
            else:

                cif_info["normalized_formula"] = cif_info.name.apply(normalize_formula)
                dfs = group_data(
                    cif_info, column_name=["normalized_formula", "spacegroup_number"]
                )
        else:
            dfs = group_data(cif_info, column_name=["system_type", "spacegroup_number"])
        chid = []
        for df in dfs:
            df.sort_values(by="Cij", ascending=False, inplace=True)
            # Get the id corresponding to the highest correlated.
            # Here we are hoping that those patterns are all the same (fingers-crossed)
            material_id = df.id.iat[0]
            chid.append(material_id)

        df_ = cif_info[cif_info.id.isin(chid)].copy()  # create new df
        if save:
            self.modified_df = df_.copy()
            self.chosen_ids = chid
            return df_
        else:
            return df_

    def calculate_correlations(
        self, select_cifs: bool = True, by_sys: bool = True
    ) -> str:
        """
        Calculate correlations between experimental data and obtained crystal structures
        Parameters
        ----------
        select_cifs : If set True, it will select one crystal structure out of a group of similar (same spacegroup/lattice)
         by taking the one with highest correlation with experiemntal data. This is to avoid a run only having the
         "correct" pattern from many different sources.

        Returns
        -------
        A string with the filename with the correlations with the simulated patterns with the experimental pattern

        """

        exp_data = self.exp_data.copy()
        df_data = [exp_data.loc[:, "int"]]
        filenames = [os.path.basename(self.exp_data.filename.iat[0])]
        self.cif_info.reset_index(drop=True, inplace=True)
        for path in self.cif_info.loc[:, "simulated_files"]:
            df_data.append(pd.read_csv(path).loc[:, "int"])
            filenames.append(os.path.basename(path))
        intensities = pd.concat(df_data, axis=1)  # concat intensities
        intensities.dropna(inplace=True)  # drop anyna just incase
        intensities.columns = filenames  # set the columns for filenames
        Cij = intensities.corr()
        Cij_ = Cij.iloc[1:, 0]
        self.cif_info["Cij"] = np.array(Cij_)

        # Try to filter out "same" CIFs (otherwise we just get a bunch of high Cij from the same phase)
        if select_cifs:
            self.cif_info = self.select_cifs(
                save=False, normalize=True, by_systemtype=by_sys
            )
        folder = self.working_folder
        name_out = os.path.join(
            folder,
            os.path.normpath(folder).replace(os.sep, "_") + "_Correlations_run_1.csv",
        )
        # Export already sorted
        self.cif_info.sort_values(by="Cij", ascending=False).to_csv(name_out)
        print(
            f"""Highest Correlated pattern is {self.cif_info.sort_values(by='Cij', ascending=False).name.iat[0]}, with Cij: {self.cif_info.sort_values(by='Cij', ascending=False).Cij.iat[0]}"""
        )
        self.corrinfo = name_out
        return name_out

    def analyze(
        self,
        n_runs: Any[int, str] = "auto",
        grabtop: int = 3,
        delta: float = 1.3,
        combine_filter: bool = False,
        select_cifs: bool = True,
        plot_all: bool = False,
        ignore_provider: List[str] = ("AFLOW",),
        ignore_comb: List[str] = None,
        is_ceramic: bool = False,
        ignore_ids: List[str] = None,
        solver: str = "box",
        group_method: str = "system_type",
        auto_threshold: int = 10,
        r_ori: bool = False,
        n_jobs: int = -1,
    ) -> pd.DataFrame:
        """
        Search for possible phases in a experimental pattern given the possible elements
        Parameters
        ----------
        n_runs : Number of phases to search, including main phase
        grabtop : At each run, compare how many patterns?
        delta : Width of peak removal (bragg)
        combine_filter : When doing combinations, allows the code to search a grabtop+1 position with a pattern has
        already shown up in a previous run?
        select_cifs : Defaults to True. See calculate_correlation documentantion
        plot_all : To export all refinement plots (defaults to False)
        ignore_provider : A list of providers to ignore. Default: ["AFLOW"], due to the large amount of theoretical crystal structures. Can be manually turned on.
        ignore_comb : A list of combinations to ignore. Eg: ["B-O"], would ignore searching for B-O oxide.
        is_ceramic : make auto-filling of ignore_comb if it is ceramic material
        ignore_ids: A list of possible unique IDs to ignore. Defaults to None.
        solver: Decide which solver to use. Defaults to box method "box". For residual method use "rietveld"
        group_method: Decides how to try to group similiar crystal structures. Defaults to "system_type".
        For stoichmetry based grouping use "stoich".
        auto_threshold: Threshold for when to stop box method when n_runs is set to auto. Defaults to 10 (percent)
        r_ori: Allows the `rietveld` search method to also try to account for texture (pref. orientation). Defaults to False
        n_jobs: How many proccesses to use when simulating / refining combinations (box)
        Returns
        -------
        A pandas dataFrame with the search results.
        """

        # Setup to max CPU count.
        if n_jobs == -1:
            n_jobs = os.cpu_count()

        if solver not in ["box", "rietveld"]:
            raise ValueError
        if group_method not in ["system_type", "stoich"]:
            raise ValueError
        elif group_method == "system_type":
            systype = True
        else:
            systype = False
        if self.use_preprocessed:
            self.exp_data_file = self._preprocess_path
            self.datafmt = "csv"
            print(
                f"Using preprocessed data {self._preprocess_path}. New datafmt is: {self.datafmt}"
            )

        if n_runs == "auto":
            auto = True
            n_runs = 100
        else:
            auto = False
            
        if is_ceramic:
            to_ignore = make_system_types([ele for ele in elements if ele != "O"])
            if ignore_comb:
                ignore_comb.extend(to_ignore)
            else:
                ignore_comb = to_ignore

        # Get the cifs, simulate the patterns, run correlation (first phase)
        self.get_cifs(
            ignore_provider=ignore_provider,
            ignore_comb=ignore_comb,
            ignore_ids=ignore_ids,
        ).simulate_all(n_jobs=n_jobs).calculate_correlations(
            select_cifs=select_cifs, by_sys=systype
        )

        # Plot highest correlated first.
        self.plot_highest_correlated()
        # Remove the patterns and get the information of each run etc:

        if solver == "rietveld":
            if n_runs == "auto":
                raise (
                    NotImplementedError,
                    "n_runs = `auto` not available for residual method.",
                )
            runs, topn, topnfilter = run_correlation_analysis_riet(
                experimental_data=self.exp_data_file,
                patterns_data=self.cif_info,
                number_of_runs=n_runs,
                datafmt=self.datafmt,
                grabtop=grabtop,
                working_folder=self.working_folder,
                allow_ori=r_ori,
            )

        else:
            runs, topn, topnfilter = run_correlation_analysis(
                experimental_data=self.exp_data_file,
                patterns_data=self.cif_info,
                delta=delta,
                number_of_runs=n_runs,
                auto=auto,
                datafmt=self.datafmt,
                grabtop=grabtop,
                working_folder=self.working_folder,
                remove_background=self.rb,
                poly_degree=self.pd,
                auto_threshold=auto_threshold,
            )
            if auto:
                n_runs = len(runs)

        topnfilter = [pd.DataFrame(data) for data in topnfilter]
        if n_runs == 1:
            df_final = topn[0]
            df_final["gpx_path"] = df_final.apply(
                lambda row: make_gpx(row["name"], row.spacegroup, row.full_path), axis=1
            )
            df_final.sort_values(by="rwp", inplace=True)
            df_final.reset_index(drop=True, inplace=True)
            df_final.to_csv(
                os.path.join(self.working_folder, "Refinement Results.csv"), index=False
            )
            df_final.drop(["bragg"], axis=1, inplace=True)
            print("Analysis complete")
            best_rwp = df_final["rwp"].iat[0]
            best_gpx = df_final["gpx_path"].iat[0]
            best_cifs = df_final["filename"].iat[0]
            print(f"Best result {best_cifs}, Rwp: {best_rwp:.2f}%")
            self.rwp = best_rwp
            self.gpx_best = best_gpx
            self.cifs_best = best_cifs
            self.results = df_final
            self.plot_result(0, engine="plotly", mode="simul", save=True)
            print(
                f"Saved final plot result as {os.path.basename(best_gpx).replace('gpx', 'html')}"
            )
        else:

            make_plot_all(runs, topn, self.name, self.working_folder)
            make_plot_step(runs, topn, self.name, self.working_folder, solver=solver)
            if solver == "box":
                make_plot_step(
                    runs, topnfilter, self.name + "_filter", self.working_folder
                )

                # In case of using "filter" option (remove structure that has already appeared in n-k run, k > 0
                if combine_filter:
                    trials = combine_pos(top_patterns=topnfilter)
                else:
                    trials = combine_pos(top_patterns=topn)
                    # If dont use filter, attempt to clean up repeated runs and runs with same structure.
                    aux_df = pd.DataFrame(trials)
                    # Get ids from each combination (id is a unique identifier of structure of each database)
                    aux_df["ids"] = aux_df.comb.apply(
                        lambda comb: [dic["id"] for dic in comb]
                    )
                    # Remove ids where there is more than the same id in the list
                    aux_clean = aux_df[
                        aux_df.ids.apply(
                            lambda lst: all(
                                [lst.count(element) == 1 for element in lst]
                            )
                        )
                    ].copy()
                    # Clean up now runs that are identical
                    aux_clean["ids_str"] = aux_clean["ids"].apply(
                        lambda lst: ",".join(lst)
                    )
                    aux_clean.drop_duplicates(subset="ids_str", inplace=True)
                    aux_clean.reset_index(drop=True, inplace=True)
                    print(
                        f"Removed {len(aux_df) - len(aux_clean)} repeated combinations."
                    )
                    # Remove it
                    aux_clean.drop(["ids", "ids_str"], axis=1, inplace=True)

            if solver == "box":
                # Set up the dataframe to refine the patterns in parallel using modin
                if combine_filter:
                    to_run = pd.DataFrame(trials)
                else:
                    to_run = aux_clean

                # Set up arguments for multiprocessing
                data = self.exp_data_file
                wf = self.working_folder
                comb_run = to_run.to_dict(orient="records")
                args = [[f] + [data, wf] for f in comb_run]

                # Refine all combinations in parallel
                with Pool(processes=n_jobs) as p:
                    result = p.starmap(refine_comb, args)
                    p.close()
                    p.join()

                # Parse results and report.
                to_run["rwp"] = [r["rwp"] for r in result]
                to_run["wt"] = [r["wt"] for r in result]
                outpath = os.path.join(self.working_folder, "trials_result.json")
                to_run.sort_values(by="rwp", inplace=True)
                trials = to_run.to_dict(orient="records")

                # Save raw results in json
                save_json(trials, outpath)

                # Load first run from correction analysis and concatenate all
                first_trial = load_json(
                    os.path.join(self.working_folder, "first_run.json")
                )
                first_df = pd.DataFrame(first_trial)
                first_df["nruns"] = 1
                first_df["pos"] = 1
                trials_df = pd.DataFrame(trials)
                df_list = []
                # aggregation of data for final report
                for adict, rwp, wt in zip(trials_df.comb, trials_df.rwp, trials_df.wt):
                    newdict = {}
                    for k, v in adict[0].items():
                        newdict[k] = [[x[k] for x in adict]]
                    df_ = pd.DataFrame(newdict)
                    df_["rwp"] = rwp
                    df_["wt"] = [wt]
                    df_["nruns"] = len(wt)
                    df_list.append(df_)

                # Setup final dataframe.
                df_final = pd.concat([first_df] + df_list)
                df_final["gpx_path"] = df_final.apply(
                    lambda row: make_gpx(row["name"], row.spacegroup, row.full_path),
                    axis=1,
                )
                df_final.sort_values(by="rwp", inplace=True)
                df_final.reset_index(drop=True, inplace=True)
                df_final.to_csv(
                    os.path.join(self.working_folder, "Refinement Results.csv"),
                    index=False,
                )
                print("Analysis complete")
                best_rwp = df_final["rwp"].iat[0]
                best_gpx = df_final["gpx_path"].iat[0]
                best_cifs = df_final["filename"].iat[0]
                print(
                    "Best combination: {}, Rwp: {}%".format(
                        "_".join(best_cifs), best_rwp
                    )
                )

            elif solver == "rietveld":
                todrop = ["simulated_files", "simulated_reflects", "gpx"]
                columns_stack = [
                    "provider",
                    "id",
                    "name",
                    "filename",
                    "spacegroup",
                    "spacegroup_number",
                    "crystal_system",
                    "system_type",
                    "full_path",
                    "Cij",
                ]
                for df in topn:
                    df.drop(todrop, axis=1, inplace=True)
                    df.reset_index(inplace=True, drop=True)
                for i in range(1, len(topn)):
                    best_last = (
                        topn[i - 1]
                        .loc[0, columns_stack]
                        .apply(lambda e: e if type(e) == list else [str(e)])
                    )
                    for k in range(len(topn[i])):
                        topn[i].loc[k, columns_stack] = best_last + topn[i].loc[
                            k, columns_stack
                        ].apply(lambda e: [str(e)])
                df_final = pd.concat(topn)
                df_final.sort_values(by="rwp", inplace=True)
                df_final.reset_index(drop=True, inplace=True)
                df_final.to_csv(
                    os.path.join(self.working_folder, "Refinement Results.csv"),
                    index=False,
                )
                print("Analysis complete")
                best_rwp = df_final["rwp"].iat[0]
                best_gpx = df_final["gpx_path"].iat[0]
                best_cifs = df_final["filename"].iat[0]
                print(
                    "Best combination: {}, Rwp: {}%".format(
                        "_".join(best_cifs), best_rwp
                    )
                )

        self.rwp = best_rwp
        self.gpx_best = best_gpx
        self.cifs_best = best_cifs
        self.results = df_final
        if self.use_preprocessed:
            print(f"Returning data to original format {self._old_fmt}")
            self.exp_data_file = self._old_data
            self.datafmt = self._old_fmt

        _outf = os.path.join(self.working_folder, "fig")
        create_folder(_outf)
        gpx_files = os.path.join(self.working_folder, "gsas2_files")

        # To plot all refinements or not. Attempt to do in parallel by Ray as it can be timeconsuming (parsing gpxs..)
        if plot_all:
            plot_all_gpx(gpx_files + os.sep)
        return df_final

    def plot_highest_correlated(
        self,
        top: int = 30,
        offset: float = 0.3,
        imult: float = 1.0,
        width: int = 1280,
        height: int = 968,
        **kwargs,
    ) -> XRay:
        """
        Plot the top N highest correlated patterns with given experimental data using Plotly.
        Parameters
        ----------
        top : Number of patterns to plot (int), defaults to 30
        offset : Spacing between the plots, defaults to 0.3 (yoffset = y - (ki + 1) * offset) where ki is the kth pattern
        imult : Intensity multiplier. (yintensity = yoffset * imult)
        width : Figure width
        height : Figure height
        kwargs : Keyword arguments to be passed to px.line

        Returns
        -------
        self
        """
        # Set up the "simulated" dataframe in a easy way for plotting
        self._read_simulations()
        dummy_data = [
            ("name", "Exp.Data"),
            ("provider", "Xray machine"),
            ("mid", 0),
            ("spacegroup", "?"),
            ("spacegroup_number", "?"),
            ("crystal_system", "?"),
            ("Cij", 1),
        ]
        col_order = self.simulated_gsas2[0][0].columns
        exp_data = self.exp_data.copy()
        for col, val in dummy_data:
            exp_data[col] = val
        patterns = [t[0] for t in self.simulated_gsas2]
        temp = self.cif_info.copy()
        temp.sort_values(by="Cij", inplace=True)
        temp.reset_index(drop=True, inplace=True)
        # j = self.cif_info.loc[0, c_to_grab]
        exp_data = exp_data[col_order]
        patterns_sorted = sorted(patterns, key=lambda c: c.Cij.iat[0], reverse=True)
        top_patterns = patterns_sorted[:top]
        make_offset(top_patterns, offset=offset, imult=imult)

        # Patterns scaling based on exp data
        SCALER_MAX = exp_data.int.max()
        exp_data.int = exp_data.int / SCALER_MAX
        plotly_df = pd.concat([exp_data] + top_patterns)
        hover_data = ["spacegroup", "crystal_system", "spacegroup_number", "Cij"]
        figure = px.line(
            plotly_df,
            x="theta",
            y="int",
            color="filename",
            width=width,
            height=height,
            hover_data=hover_data,
            **kwargs,
        )
        config = dict({"scrollZoom": True})
        ## some default plotly options ##
        hoverlabel = dict(font_size=18, font_family="Arial")
        legend = {"font": dict(size=18)}

        # Edit First and Second Plot
        figure.data[0]["mode"] = "markers"
        figure.data[0].marker.color = "black"
        figure.data[0].marker.symbol = "circle-open"
        figure.data[0].name = "Exp. Data"
        figure.update_layout(hoverlabel=hoverlabel, legend=legend, modebar_add=to_add)
        name = os.path.basename(self.filename)
        name = name.replace(".ras", "")
        outfolder = self.working_folder
        outname = name + "_Top" + str(top) + ".html"
        outpath = os.path.join(outfolder, outname)
        print(outname)
        figure.write_html(outpath, config=config)
        # figure.show()
        return self

    def plot_result(
        self,
        index: int,
        engine: str = "plotly",
        mode: str = "simul",
        save: bool = False,
        **kwargs
    ) -> Any[None, object]:
        """
        Plot a result of a specific index of the results dataframe
        Parameters
        ----------
        index : Index to plot (int)
        engine : Which plotting engine to use, "matplotlib" for matplotlib or "plotly" for Plotly. Defaults to "matplotlib".
        type :  "bragg" for plotting bragg positions as ticks or "simul" for full simulation patterns

        Returns
        -------
        None or a matplotlib figure
        """
        if engine == "matplotlib":
            plot_gpx(
                path=self.results.gpx_path.iat[index],
                save=save,
                engine=engine,
                outfolder=self.working_folder,
            )
        else:
            if mode == "bragg":
                return plot_gpx(
                    path=self.results.gpx_path.iat[index],
                    save=save,
                    engine=engine,
                    outfolder=self.working_folder,
                )
            else:
                simul_folder = os.path.join(self.working_folder, "Simul")
                names = self.results["name"].iat[index]
                spacegroups = self.results["spacegroup"].iat[index]
                cif_names = self.results["filename"].iat[index]

                # If its a string, we setup as a list (required for plotting func.) and make the new name.
                if type(cif_names) == str:
                    cif_names = [cif_names]
                    phasenames = [f"{names}-{spacegroups.replace('/', '-')}"]
                else:
                    phasenames = [
                        n + "-" + spg.replace("/", "-")
                        for n, spg in zip(names, spacegroups)
                    ]
                # now we make simul names
                # now we make simul names
                simul = [
                    pd.read_csv(
                        os.path.join(simul_folder, file.replace(".cif", "_Simul.csv"))
                    )
                    for file in cif_names
                ]
                for i, df in enumerate(simul):
                    df["phase_name"] = phasenames[i]
                return plot_gpx(
                    path=self.results.gpx_path.iat[index],
                    save=save,
                    engine=engine,
                    s_info=simul,
                    type=mode,
                    outfolder=self.working_folder,
                    **kwargs
                )

    def initialize_optimizer(self, index: Any[List[int], int]) -> XRay:
        """
        Initialize an optmizer for a index or a list of indexes of the result dataframe.
        This function provides a option for easily joining different set of results for a optimizaiton of phases gaven
        by different combinations. (See a list of indexes below)

        For only one index: (No union of phases to refine)
            - Specify index = index number

        A list of indexes:
            - This option is for when a union of different phases is required. For example, lets suppose that
            for a two phase search, the top 2 lowest rwp results (indexes 0 and 1) HoB2, Ho and HoB2, Ho2O3. When investigating the results,
            that is, the quickly refined patterns, it was evident that that are actually three phases, consisting of
            HoB2 + Ho and Ho2O3.
            For this particular case, if it wished to run the optmizer for a three-phase sample, one can pass to this
            function a list of indexes, which then the structures to be used will simple be a set union of both indexes,
            leading to HoB2,Ho and Ho2O3.
            Thus, in this case, index can receive a list of indexes, in this way, index=[0,1].
        Parameters
        ----------
        index : One index or a list of indexes.

        Returns
        -------
        Self.
        """
        # first lets check if index is a list or
        if type(index) == list:
            l = []
            for i in index:
                l += self.results.full_path.iat[i]
            cifs = list(set(l))
        else:
            cifs = self.results.full_path.iat[index]
            if type(cifs) != list:
                cifs = [cifs]

        self.optimizer = GSASRefiner(
            expdata=self.exp_data_file, working_folder=self.working_folder, cifs=cifs
        )
        self.optimizer.quick_refinement()  # show up default

        return self

    def run_optimizer(
        self,
        plot_best: bool = True,
        show_lattice: bool = True,
        random_state: int = 42,
        n_jobs: int = -1,
        n_trials: int = 200,
        n_startup: int = 20,
        allow_pref_orient: bool = True,
        allow_atomic_params: bool = False,
        allow_broad: bool = False,
        allow_strain: bool = False,
        allow_angle: bool = False,
        force_ori: bool = False,
        verbose: str = "silent",
        param: str = "rwp",
    ) -> XRay:
        """
        Start the Blackbox Optimization as suggested by NPJ Comp. Mat. XX xXX (2020)
        Parameters
        ----------
        plot_best : Plot the best trial result. Defaults to True.
        show_lattice : Show the best lattice params result. Defaults to True
        random_state : Optuna random seed for optmization.
        n_jobs : Number of jobs for optimization. Defaults to -1 (all cores)
        n_trials : Number of trials to search for best refinement. Defaults to 200
        n_startup : Number of startup trials of Optuna. Defaults to 20
        allow_pref_orient : Allow Preferred Orientation to be considered for optmization. Defaults to True
        allow_atomic_params : Allow refinement of atomic params (U, XU, XUF) to be considered for optmizaiton. Defaults to False
        allow_broad : Allow broadening terms to be considered for optmization. Defaults to False.
        allow_strain : Allow strain terms to be considered for optmization. Defaults to False
        allow_angle : Allow acute angle to be optmized (Small angle refinement before all other refinemenets followed by full angle). Defaults to False
        force_ori : To force prefered orientation to be always considered. Defaults to False.
        verbose : Suppress optuna messages. "silent" for yes, None to no
        param : Optmization objective param. "rwp" for Rwp or "gof" for GoF. Defaults to "rwp".

        Returns
        -------
        Self.
        """

        if self.optimizer is not None:
            self.optimizer.fit_optimize(
                n_jobs=n_jobs,
                n_startup=n_startup,
                n_trials=n_trials,
                random_state=random_state,
                allow_angle=allow_angle,
                allow_broad=allow_broad,
                allow_pref_orient=allow_pref_orient,
                allow_atomic_params=allow_atomic_params,
                allow_strain=allow_strain,
                force_ori=force_ori,
                verbose=verbose,
                param=param,
            )

        if plot_best:
            self.optimizer.plot_result()
        if show_lattice:
            print(self.optimizer.lattice_best)

        self.optgpx = os.path.join(self.working_folder, "optuna", self.optimizer.best)
        self.optlat = self.optimizer.lattice_best
        self.optrwp = self.optimizer.optim.rwp_best
        return self

    def export_results(self, outfolder="export") -> XRay:
        """
        Export the obtained results AFTER an optimization attempt has been made.
        Exports:
         - Final crystal structure params (lattice constants, volume, density, phase fractions and rwp values)
         - Yobs, ydiff and ycalc (normalized by yobs.max())
         - Bragg positions of each phase.
         - GSASII final .lst file for the best trial.
        Parameters
        ----------
        outfolder : Foldername to save

        Returns
        -------
        Self.

        """
        if self.optgpx == None:
            return "Please run optimizer first."
        else:
            self.optimizer.export_results(outfolder=outfolder)


    def load_results(self, filename: str) -> None:
        """Load past run results.
           This allows to use Xerus plotting functions to investigate the results, or re-run the analysis
           with different conditions.
           Also allows to easily start optimization if necessary

        Parameters
        ----------
        filename : str
            Path to a 'Refinement Results.csv' of an old run
        """
        def convert_data(data: str) -> Union[str, List]:
            """Helper function to convert a string into the correct data structure..

            Parameter
            ----------
            data : str
                A string representation of the data.

            Returns
            -------
            Any[str, List]
                Returns the correct representation of a string into the correct structure.
                For example "['HoB2', 'HoB4']" (str) -> ['HoB2', 'HoB4'] (List[str])
            """
            if not isinstance(data, str):
                return data
            try:
                return eval(data)
            except (NameError, SyntaxError, TypeError) as e:
                return data

        # Read the data and convert
        data = pd.read_csv(filename)
        self.results = data.applymap(convert_data)
        
    
