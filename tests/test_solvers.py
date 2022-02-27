import sys, os
from pathlib import Path
from Xerus.settings.settings import  GSAS2_BIN
from . import INSTALL_PATH
sys.path.append(GSAS2_BIN)
from Xerus import XRay
import pytest
import shutil


@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_boxauto():
    """
    Test boxauto
    Returns
    -------
    1 if works.
    """

    r = XRay(
        name="HoB2",
        elements=["Ho", "B"],
        exp_data_file=INSTALL_PATH / "data/HoB21.ras",
        data_fmt="ras",
        working_folder="tests_box_auto/",
        max_oxy=2,
        standarize_int=True,
        use_preprocessed=True,
        remove_background=True
    )
    r.analyze(n_runs="auto", solver="box", ignore_provider=["AFLOW", "COD"])
    assert r.results.name.iat[0] == ['HoB2', 'HoB4']
    shutil.rmtree(r.working_folder)

@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_boxmethod():
    """
    Test box method with n_runs=2
    Returns
    -------
    1 if works.
    """

    r = XRay(
        name="HoB2",
        elements=["Ho", "B"],
        exp_data_file=INSTALL_PATH / "data/HoB21.ras",
        data_fmt="ras",
        working_folder="tests_box_method/",
        max_oxy=2,
        standarize_int=True,
        use_preprocessed=True,
        remove_background=True
    )
    r.analyze(n_runs=2, solver="box", ignore_provider=["AFLOW", "COD"])
    assert r.results.name.iat[0] == ['HoB2', 'HoB4']
    shutil.rmtree(r.working_folder)

@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_residualmethod():
    """
    Test rietveld method with n_runs=2
    Returns
    -------
    1 if works.
    """

    r = XRay(
        name="HoB2",
        elements=["Ho", "B"],
        exp_data_file=INSTALL_PATH / "data/HoB21.ras",
        data_fmt="ras",
        working_folder="tests_residual_method/",
        max_oxy=2,
        standarize_int=True,
        use_preprocessed=True,
        remove_background=True
    )
    r.analyze(n_runs=2, solver="rietveld", ignore_provider=["AFLOW", "COD"])
    assert r.results.name.iat[0] == ['HoB2', 'HoB4']
    shutil.rmtree(r.working_folder)

@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_onephaserun():
    """
    Test n_runs=1
    Returns
    -------
    1 if works.
    """

    r = XRay(
        name="HoB2",
        elements=["Ho", "B"],
        exp_data_file=INSTALL_PATH / "data/HoB21.ras",
        data_fmt="ras",
        working_folder="tests_box_method/",
        max_oxy=2,
        standarize_int=True,
        use_preprocessed=True,
        remove_background=True
    )
    r.analyze(n_runs=1, ignore_provider=["AFLOW", "COD"])
    assert r.results.name.iat[0] == 'HoB2'
    shutil.rmtree(r.working_folder)







