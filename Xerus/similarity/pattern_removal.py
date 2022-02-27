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
import pandas as pd
import numpy as np
from Xerus.engine.gsas2riet import quick_gsas
from Xerus.readers.datareader import DataReader
from Xerus.settings.settings import INSTR_PARAMS
from Xerus.utils.tools import save_json, sort_and_reset
from Xerus.utils.preprocessing import remove_baseline
from itertools import product, combinations
from typing import List, Tuple, Any, Dict
from pandas import DataFrame

instr_params = INSTR_PARAMS
to_export = ['provider',
             'id',
             'name',
             'filename',
             'spacegroup',
             'spacegroup_number',
             'crystal_system',
             'system_type',
             'full_path',
             'rwp']


def make_pos_array_and_diff(trials: List[dict]) -> List[dict]:
    """

    Parameters
    ----------
    trials : Non-filtered refinement trials

    Returns
    -------
    A list of dictionaries with updated position and difference of neighbours
    """

    _trials = trials[:]
    for i, c in enumerate(trials):
        x_array = sorted([d['pos'] for d in trials[i]['comb']])
        _trials[i]['pos_array'] = x_array
        _trials[i]['diff_array'] = [x_array[i] - x_array[i + 1] for i in range(len(x_array) - 1)]
        # print(x_array)
        # print(len(trials))
        # trials[i]['comb']
    return _trials


def filter_comb(all_comb: List[dict]) -> List[dict]:
    """
    Filter combinations based on their positional arrays.
    If there is any value <-1 in the position array and the first element is etiher differen than 0 or 1, it will be skipped.
    Parameters
    ----------
    all_comb : A list of dictionarise

    Returns
    -------
    A filtered list of dictionaries based on their distance.

    """
    return_list = []
    total_combs = len(all_comb)
    print("Filtering combinations..")
    for comb in all_comb:
        p_array = comb['pos_array']
        p_diff = comb['diff_array']
        if p_array[0] != 1:
            if p_array[0] != 2:
                # print("Pos diff then 0 or 1!")
                continue
        if any([val < -1 for val in p_diff]):
            # print("There is a diff < -1!")
            continue

        return_list.append(comb)
    final_len = len(return_list)
    rem = total_combs - final_len
    print("Initial # of combs: {}".format(total_combs))
    print("Final # of combs: {}, # removed: {}".format(final_len, rem))
    return return_list


def combine_pos(top_patterns: List[pd.DataFrame]) -> List[dict]:
    """
    Combine found structures for refinement based on the distance-array approach
    Parameters
    ----------
    top_patterns : Highest correlated patterns found in each run.

    Returns
    -------
    A list of dictionaries with the structures to be refined.
    """


    # _trials_name = "trials.json"
    # _trials_path = os.path.join(working_folder, _trials_name)
    runs = [df.drop(['bragg',
                     'simulated_files', 'simulated_reflects',
                     'Cij', 'rwp'], axis=1).to_dict(orient="records") for df in top_patterns]

    def transform(tup):
        """
        transform the [({dic1}, {dic2}), {dic3}}] -> [(dic1, dic,di3)]
        """
        res = tuple()
        for item in tup:
            if type(item) == dict:
                item = (item,)
            else:
                item = tuple(item)
            res += item
        return res

    # setup position array
    for pos in range(0, len(runs)):
        for phase in runs[pos]:
            # print(phase)
            phase['pos'] = pos + 1
    # this might be unnecessary.
    runs_copy = runs[:]
    final_ = []

    for i in range(1, len(runs_copy)):
        runs = runs_copy[:i + 1]
        winner = [runs[0][0]]  # get winner
        rest = runs[0][1:]  # get leftover
        others = runs[1:]  # runs N>1
        final = []
        a0 = list(product(winner, rest))
        if len(runs) == 2:
            others = list(product(winner, *runs[1:]))
            final = a0 + others
        # Messy logic below
        elif len(runs) > 2:
            runs = runs_copy[:i + 1]
            top1 = [runs[0][0]]  # grab the top1
            rest_first = runs[0][1:]  # grab the rest ofthe first run n=1
            other_phases = runs[1:]  # grab the remaining runs information
            unpack = [item for sublist in other_phases for item in sublist]  # flatten the other runs
            to_combine = rest_first + unpack  # put together the other runs
            combs = list(combinations(to_combine, len(runs) - 1))  # Combine the remaining patterns
            combs = [{"comb": comb} for comb in combs]
            combs = make_pos_array_and_diff(combs)
            removed = [dic['comb'] for dic in filter_comb(combs)]
            final = [transform(comb) for comb in
                     list(product(top1, removed))]  # combine the top1 with the other combina= final
        final_ += final
    trials = [{"comb": comb} for comb in final_]
    return trials


def remove_bragg(exp_data: pd.DataFrame, bragg_pos: set, delta: float) -> List[int]:
    """
    Return the indices corresponding to the bragg positions to be removed

    Parameters
    ----------
    exp_data : A pandas dataframe containing the experimental data
    bragg_pos : A set of consiting of bragg positions to be removed
    delta :  Width of peak to be removed

    Returns
    -------
    A list of indexes to be removed from the experimental data pattern.
    """
    modify = exp_data.copy()  # make a copy so we dont change the raw exp. data
    ranges = [(angle - delta / 2, angle + delta / 2) for angle in
              bragg_pos]  # make the lsit of angles to remove
    for rmin, rmax in ranges:
        modify = modify[
            (modify.theta <= rmin) | (modify.theta >= rmax)]  # Get the data that is not inside a bragg pos +- delta/2
    return modify.index.tolist()


def calculate_Cij(patterns: List[pd.DataFrame], exp_data: pd.DataFrame) -> Tuple[np.array, pd.DataFrame]:
    """
    Calculate the Pearson`s correlation (Cij) between the exp data and the simulated patterns.
    Parameters
    ----------
    patterns : A list of dataframes containing each simulated pattern
    exp_data : A pandas dataframe containing the experimental pattern.

    Returns
    -------
    A Tuple with numpy array with the Pearson`s correlation of each pattern in relation to the experimental pattern and
    the full correlation dataframe
    """

    # Stack patterns in columns
    correlation_df = pd.concat([exp_data] + patterns, axis=1)
    # Remove first entry as it represents the exp. data
    Cij = correlation_df.corr().iloc[1:, 0]
    return np.array(Cij), correlation_df


def remove_first_Cij(df: pd.DataFrame, exp_data: pd.DataFrame, delta: float, bragg_old: set) -> Tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Receives a pandas DataFrame of size N where at pos 0 lies the "Top1" pattern,
    the exp data dataframe, the delta (width) to remove  and the already removed bragg positions
    Remove the following bragg position of current top1 pattern (exp data and simul) and then recalculates Cij
    Then returns a df of size N-1 with the old "top1" pattern already removed. (not sorted for new Cij)

    Parameters
    ----------
    df : A pandas dataframe of size N, where the current top1 patterns sits at index 0
    exp_data : The experimental data DataFrame.
    delta : The width to remove of each peak corresponding to the bragg position
    bragg_old : The already removed bragg positions.

    Returns
    -------
    A tuple of 3 dfs (new data [removed], modified exp data [removed] and the correlation dataframe)
    """
    bragg_old = set(bragg_old)  # incase
    bragg_new = df.bragg.iat[0]
    bragg = bragg_new.union(bragg_old)
    index_new = remove_bragg(exp_data, bragg, delta)
    df_new = df.iloc[1:, :].copy().reset_index(drop=True)

    def clean_up(df_new: pd.DataFrame, index_new: List[float]) -> pd.DataFrame:
        """
        Attempts to clean up patterns with almost no high intesnity peaks left.
        Parameters
        ----------
        df_new : pd.DataFrame
            The pandas dataframe containing information regarding to the current run
            Assumes that the pd.DataFrame starts from 0 -> 1000. (reseted index.)
        index_new : List[int]
            A list of integers to be read.
        Returns
        -------
        """
        idx_drop = []
        df_mod = df_new.copy()
        for i in range(len(df_mod)):

            # Read up patterns and get maximum intensity.
            pattern_path = df_new.simulated_files.iat[i]
            pattern_data = pd.read_csv(pattern_path)
            INT_MAX = pattern_data.int.max()
            # Remove braggs and calculate newest maximum.
            pattern_new = pattern_data.loc[index_new, :].copy()
            INT_MAX_NEW = pattern_new.int.max()
            # Check if drops below 10%.
            pct = 100 * (INT_MAX_NEW) / INT_MAX
            if pct <= 10:
                # print(f"Intensity dropped to {pct}%, below 10%. Removing this pattern.")
                idx_drop.append(i)

        # Edit dataframe and return new.
        print (f"Removing {len(idx_drop)} patterns that maximum intensity dropped below 10%.")
        df_mod.drop(idx_drop, inplace=True)
        df_mod.reset_index(drop=True, inplace=True)

        return df_mod

    df_new = clean_up(df_new, index_new)
    patterns = [pd.read_csv(path).loc[index_new, :].drop('theta', axis=1) for path in
                df_new.simulated_files]

    exp_mani = exp_data.copy()
    exp_mani = exp_mani.loc[index_new, :]
    exp_mani.drop('theta', axis=1, inplace=True)
    Cij, cdf = calculate_Cij(patterns, exp_mani)  # First position is always the experimental data.
    exp_data_mod = exp_data.copy().loc[index_new, :]
    df_new['Cij'] = np.array(Cij)  # update Cij
    return df_new, exp_data_mod, cdf


def filter_runs(topn: List[pd.DataFrame], n: int = 3) -> List[list]:
    """
    Filter repeated cifs in the runs and allows the next g+1 structure to be considered for refinement.
    Parameters
    ----------
    topn : Topn patterns
    n : This refers to g (grabtop)

    Returns
    -------
    List a of lists.
    """

    entries = [df.to_dict(orient='records') for df in topn]  # unpack this boys

    # loop over the entries. always comparing n->(n+1).
    # In this way we *should* be sure that the found patterns of n+1 run will not be in n and so forth
    # Wrong! If a pattern show in n, n+1 and n+2, if we do n->n+1, meaning we remove from n+1, then it will be able to repear at n+2 since n+1->n+2 means n+1 didnt have the pattern.
    saw_ids = []
    new_top = [entries[0][:n]]  # first run should not be modified
    for i in range(len(entries) - 1):
        current_tops = entries[i]  # get current top
        next_tops = entries[i + 1]  # get next one

        top_curr = current_tops[:n]  # grab current top3
        top_curr_ids = [entry['id'] for entry in top_curr]  # get the unique id of each cif
        saw_ids += top_curr_ids
        not_in_next = [entry for entry in next_tops if entry['id'] not in saw_ids]  # remove them from n+1 patterns
        entries[i + 1] = not_in_next  # save it so in next loop we have the updated version
        new_top.append(not_in_next[:n])
        # return not_in_next

    return new_top


def run_correlation_analysis(experimental_data: str, patterns_data: pd.DataFrame, delta: float, number_of_runs: int,
                             working_folder, grabtop: int = 3,
                             datafmt: str = "ras",
                             remove_background: bool = False,
                             poly_degree: int = 8,
                             auto: bool = False,
                             auto_threshold: int = 10) -> Tuple[List[Tuple[DataFrame, list, float]], List[DataFrame], List[List[Dict]]]:
    """
    Run the correlation analysis.
    First the main phase is defined by (n=0):
        - Get top 3 Cij
        - refine top 3
        - order by 'rwp'
        - top is lowest rwp
    Then the subsequent phases will be decided by the newest highest Cij after the latest "best" pattern bragg positions
    has been removed from all the patterns.


    Parameters
    ----------
    experimental_data : Path to experimental data file
    patterns_data : Patterns data.
    delta : Delta to be removed
    number_of_runs : Number of phases to search (including main phase)
    working_folder : Working folder
    grabtop : How many patterns to look at each run
    datafmt : Format of experimental datafile.
    remove_background: To remove background or not. Defaults to False
    poly_degree: If remove_backgroud is True, which degree of polynomial to remove.
    auto: If auto sets to True, it will stop when experimental data drops intensity below `auto_threshold` (%)
    auto_threshold: Percentage of maximum intensity drop to automatically stop analysis.

    Returns
    -------
    Alot of information about the patterns, correlated patterns etc.
    """
    data = patterns_data
    data['bragg'] = data.simulated_reflects.apply(lambda path: set(pd.read_csv(path)['2-theta']))  # Get bragg positions
    exp_data = DataReader().read_data(experimental_data, datafmt)[0]  # Read experimental data
    exp_data.drop('filename', axis=1, inplace=True)  # Drop the useless filename
    max_intensity = exp_data.int.max()
    # Remove background
    if remove_background:
        exp_data = remove_baseline(exp_data, poly_degree=poly_degree, plot=False)
    removed_braggs = []  # Setup the removed bragg
    bragg_list = set([])  # Setup
    top_n = []  # Setup
    topz = []  # for filter
    exp_mani = exp_data.copy()  # Copy so we dont modify.
    exp_data_mod = exp_mani.copy()
    exp_mani.drop('theta', inplace=True, axis=1)  # Remove the theta
    for n in range(number_of_runs):
        # First find main phase.
        if n == 0:
            if grabtop > len(data):
                print(f"Grabtop ({grabtop}) is higher than data length ({len(data)})")
                grabtop = len(data)
            print(f"Refining {grabtop} highest correlated patterns and determining main phase (n=0)...")
            data = sort_and_reset(data)  # Sort & Reset for descending for Cij
            data.loc[:grabtop - 1, 'rwp'] = data.loc[:grabtop - 1].apply(
                lambda row: quick_gsas(powder_data=experimental_data,
                                       instr_params=instr_params,
                                       cif_name=row.full_path,
                                       phasename=row['name'] + "-" + row.spacegroup.replace("/", "-"),
                                       outfolder=working_folder)[0], axis=1)  # OK SAVING.
            data.sort_values(by='rwp', inplace=True)  # Resort by rwp
            rwpmain = data.rwp.iat[0]  # Get which sample has lowest rwp.
            # data.drop('rwp', axis=1, inplace=True)
            first_main = data.filename.iat[0]
            current_top = data.copy()  # Get top
            top_n.append(current_top.iloc[:grabtop, :])  #
            print("Phase detected {} with rwp {}".format(first_main, rwpmain))
            topz.append(current_top)
            # Export first run info.
            out = current_top.iloc[:grabtop].loc[:, to_export].copy()
            out['wt'] = [100] * grabtop
            out = out.to_dict(orient='records')
            save_json(out, os.path.join(working_folder, "first_run.json"))
            # In case of single phase search, return only the refined.
            if number_of_runs == 1:
                return [], top_n, [] # For consistency with n>1

        if n > 0:
            # Accumulate data that is going to be removed  (n-1) data
            bragg = data.bragg.iat[0]
            removed_braggs.append((exp_data_mod, list(bragg), delta))
            # Set up the theta to remove #
            bragg_list = bragg_list.union(bragg)  # <- union means we dont have 'repeated' considered.
            # Remove and go to next n
            data, exp_data_mod, cdf = remove_first_Cij(df=data, exp_data=exp_data, bragg_old=bragg_list, delta=delta)
            # Check auto
            if auto:
                # Check if intensity dropped below 10%
                if (exp_data_mod.int.max() / max_intensity) * 100 < auto_threshold:
                    print(f"Intensity lower than {auto_threshold}%. Finishing run")
                    print(f"Number total of possible phases is {n}")
                    topn_filter = filter_runs(topz, grabtop)
                    return removed_braggs, top_n, topn_filter
            data = sort_and_reset(data)
            next_phase = data.filename.iat[0]
            print("Phase detected {}".format(next_phase))
            current_top = data.copy()
            topz.append(current_top)
            top_n.append(current_top.iloc[:grabtop, :])

    # Add last run information for visualization purposes.
    bragg = data.bragg.iat[0]
    removed_braggs.append((exp_data_mod, list(bragg), delta))
    topn_filter = filter_runs(topz, grabtop)
    return removed_braggs, top_n, topn_filter


def run_correlation_analysis_riet(experimental_data: str, patterns_data: pd.DataFrame, number_of_runs: int,
                             working_folder, grabtop: int = 3,
                             datafmt: str = "ras", allow_ori: bool = False,
                             ) -> Tuple[List[Tuple[DataFrame, list, float]], List[DataFrame], List[List[Dict]]]:
    """
    This implements the fast rietveld method.
    This method seems to be really dependent on the quality of the starting data and of the CIFs tried.
    If there is no issues (mostly during refining well and not having alot of "wrong" structures), this method can
    quickly find the phases of a pattern.

    It works as follow
    First the main phase is defined by (n=0) (same as box)
        - Get top 3 Cij
        - refine top 3
        - order by 'rwp'
        - top is lowest rwp

    Then, for the next step the second phase will be determined by:
        - Get the top 3 Cij with the Residual (Residual = Yobs - Ycalc, ycalc is from last step)
        - Multiphase refinement of the top3. Lowest rwp -> Next phase.
        - Repeat, where new residual is now Yobs - ycalc from the multiphase run

    Parameters
    ----------
    experimental_data : Path to experimental data file
    patterns_data : Patterns data.
    delta : Delta to be removed
    number_of_runs : Number of phases to search (including main phase)
    working_folder : Working folder
    grabtop : How many patterns to look at each run
    datafmt : Format of experimental datafile.
    allow_ori: To try to fit pref. orientation while doing rietveld search. Defaults to False.

    Returns
    -------
    Alot of information about the patterns, correlated patterns etc.
    """
    data = patterns_data
    exp_data = DataReader().read_data(experimental_data, datafmt)[0]  # Read experimental data
    exp_data.drop('filename', axis=1, inplace=True)  # Drop the useless filename
    removed_braggs = []  # Setup the removed bragg
    top_n = []  # Setup
    exp_mani = exp_data.copy()  # Copy so we dont modify.
    exp_data_mod = exp_mani.copy()
    exp_mani.drop('theta', inplace=True, axis=1)  # Remove the theta
    cifs = []
    pnames = []
    for n in range(number_of_runs):
        if n == 0:  # First case.
            # Safety in case we dont have enough candidates structures.
            if grabtop > len(data):
                grabtop = len(data)
            print(f"Refining {grabtop} highest correlated patterns and determining main phase (n=0)...")
            from Xerus.engine.gsas2riet import quick_gsas
            data = sort_and_reset(data)  # Sort & Reset for descending for CijSSS
            results = data.loc[:grabtop - 1].apply(
            lambda row: quick_gsas(powder_data=experimental_data,
                                       instr_params=instr_params,
                                       cif_name=row.full_path,
                                       phasename=row['name'] + "-" + row.spacegroup.replace("/", "-"),
                                       outfolder=working_folder, ori=allow_ori), axis=1)  # OK SAVING.
            rwp = [res[0] for res in results]            ## Accumulate data that is going to be removed ## (n-1) data
            gpx = [res[1] for res in results]
            gpx_path = [res[1].filename for res in results]
            # Save the information
            data.loc[:grabtop-1, 'rwp'] = rwp
            data.loc[:grabtop-1, 'gpx'] = gpx
            data.loc[:grabtop-1, 'wt'] = [100] * grabtop
            data.loc[:grabtop-1, 'gpx_path'] = gpx_path
            data.sort_values(by='rwp', inplace=True)  # Resort by rwp
            rwpmain = data.rwp.iat[0]  # Get which sample has lowest rwp.
            first_main = data.filename.iat[0]
            current_top = data.copy()  # Get top
            top_n.append(current_top.iloc[:grabtop, :])  #
            print("Phase detected {} with rwp {:.3f}".format(first_main, rwpmain))
            # Export first run info.
            out = current_top.iloc[:grabtop].loc[:, to_export].copy()
            out['wt'] = [100] * grabtop
            out = out.to_dict(orient='records')
            save_json(out, os.path.join(working_folder, "first_run.json"))
            if number_of_runs == 1:
            	# In case one run just return these.
                # What we will show is top_n, since there is no need for filter etc.
                return [], top_n, [] # For consistency with n>1

        if n > 0:
            from Xerus.engine.gsas2riet import run_gsas_mp2
            print(f'Calculating correlations against residual and refining top {grabtop}. n = {n}')
            cifs.append(data.full_path.iat[0]) # Append last run best cif to add to refinement.
            pnames.append(data['name'].iat[0] + "-" + data.spacegroup.iat[0].replace("/", "-"))

            removed_braggs.append((exp_data_mod, [], 1))
            residual = data.gpx.iat[0].histograms()[0].getdata("residual")
            res_clean = [r if r > 0 else 0 for r in residual]
            data, exp_data_mod, cdf = remove_first_Cij_riet(df=data, exp_data=exp_data, residual=res_clean)
            data = sort_and_reset(data)
            results = data.loc[:grabtop - 1].apply(
                lambda row: run_gsas_mp2(powder_data=experimental_data,
                                       cifs=[row.full_path] + cifs,
                                       phasenames=[row['name'] + "-" + row.spacegroup.replace("/", "-")] + pnames,
                                       outfolder=working_folder, return_gpx=True, ori=allow_ori), axis=1)  # OK SAVING.

            rwp = [res[0] for res in results]
            gpx = [res[1] for res in results]
            wt = [str(res[2]) for res in results]
            gpx_path = [res[1].filename for res in results]
            data.loc[:grabtop-1, 'rwp'] = rwp
            data.loc[:grabtop-1, 'gpx'] = gpx
            data.loc[:grabtop-1, 'wt'] = wt
            data.loc[:grabtop-1, 'wt'] = data.loc[:grabtop-1, 'wt'].apply(eval)
            data.loc[:grabtop-1, 'gpx_path'] = gpx_path
            data.sort_values(by='rwp', inplace=True)  # Resort by rwp
            next_phase = data.filename.iat[0]
            # Ok lets add refinement now
            print("Phase detected {}".format(next_phase))
            current_top = data.copy()  # Get current top3
            top_n.append(current_top.iloc[:grabtop, :])


    removed_braggs.append((exp_data_mod, [], 1))
    return removed_braggs, top_n, []


def remove_first_Cij_riet(df: pd.DataFrame, exp_data: pd.DataFrame, residual: np.array) -> Tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Receives a pandas DataFrame of size N where at pos 0 lies the "Top1" pattern,
    the exp data dataframe, the delta (width) to remove  and the already removed bragg positions
    Remove the following bragg position of current top1 pattern (exp data and simul) and then recalculates Cij
    Then returns a df of size N-1 with the old "top1" pattern already removed. (not sorted for new Cij)

    Parameters
    ----------
    df : A pandas dataframe of size N, where the current top1 patterns sits at index 0
    exp_data : The experimental data DataFrame.
    delta : The width to remove of each peak corresponding to the bragg position
    bragg_old : The already removed bragg positions.

    Returns
    -------
    A tuple of 3 dfs (new data [removed], modified exp data [removed] and the correlation dataframe)
    """
    df_new = df.iloc[1:, :].copy().reset_index(drop=True)
    patterns = [pd.read_csv(path).drop('theta', axis=1) for path in
                df_new.simulated_files]


    exp_mani = exp_data.copy()
    exp_mani['int'] = residual
    exp_mani.drop('theta', axis=1, inplace=True)
    Cij, cdf = calculate_Cij(patterns, exp_mani)  # First position is always the experimental data.
    exp_data_mod = exp_data.copy()
    exp_data_mod.loc[:, 'int'] = residual
    df_new['Cij'] = np.array(Cij)  # update Cij
    return df_new, exp_data_mod, cdf
