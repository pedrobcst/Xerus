import codecs
import io
import os
from typing import List, Union

import pandas as pd

from conf import AppSettings


def read_input(file: io.BytesIO) -> str:
    """Helper function to save the input file into a temporary folder.

    Parameters
    ----------
    file : io.BytesIO
        An UploadedFile from Streamlit st.file_uploader.

    Returns
    -------
    _type_
        Nothing. Save files into the temporary folder.
    """
    stringio = io.StringIO(file.getvalue().decode("shift-jis"))

    # Save
    with codecs.open(os.path.join(AppSettings.TMP_FOLDER, file.name), "w", "shift-jis") as f:
        f.write(stringio.read())

    # Return path
    return os.path.join(AppSettings.TMP_FOLDER, file.name)



def process_input(input: str, return_int: bool = False) -> Union[None, List]:
    """Helper function to process streamlit inputs.

    Parameters
    ----------
    input : str
        A streamlit input that is seperated by comma.

    Returns
    -------
    Union[None, List]
        Returns a list if the input is not empty. Otherwise, returns None.
    """
    if input:
        if return_int:
            return_list = [int(element) for element in input.split(",")]
            if len(return_list) > 1:
                return return_list
            else:
                return return_list[0]
        return input.split(",")
    else:
        return None

    
def make_simulation_df(cif_info: pd.DataFrame) -> pd.DataFrame:
    """_summary_

    Parameters
    ----------
    cif_info : pd.DataFrame
        _description_

    Returns
    -------
    pd.DataFrame
        _description_
    """
    to_keep = [
        "filename",
        "spacegroup",
        "crystal_system",
        "spacegroup_number",
        "Cij"
    ]
    df = cif_info.copy()
    dfs = []
    for i, row in df.iterrows():
        simul_data = pd.read_csv(row['simulated_files'])
        for tk in to_keep:
            simul_data[tk] = row[tk]
        dfs.append(simul_data)
    return dfs



def make_selection_label(idx: int, df: pd.DataFrame) -> str:
    """Make a label for the selection.
    """
    name = ""
    if isinstance(df.loc[idx,"name"], list):
        for phase, spacegroup in zip(df.loc[idx, "name"], df.loc[idx, "spacegroup"]):
            name += f"{phase} ({spacegroup}) / "
        name += f" Rwp: {df.loc[idx, 'rwp']:.3f}% (id:{idx})"
    else:
        name = f"{df.loc[idx, 'name']} ({df.loc[idx, 'spacegroup']}) / Rwp: {df.loc[idx, 'rwp']:.3f}% (id:{idx})"
    return name
