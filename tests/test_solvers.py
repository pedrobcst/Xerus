import pytest
import shutil
import unittest.mock
import pymatgen
import monty.serialization
import mongomock
from . import INSTALL_PATH
from Xerus import XRay


@pytest.fixture
def working_folder(request):
    """A fixture to set up working folders named according to the test name."""
    folder_name = request.node.originalname
    shutil.rmtree(folder_name, ignore_errors=True)
    yield folder_name
    shutil.rmtree(folder_name, ignore_errors=True)


@pytest.fixture(scope="function")
def mock_mp_rester():
    """A mock fixture for querying the Ho-B chemical system from the MP API."""
    with unittest.mock.patch("pymatgen.MPRester") as mock:
        instance = mock.return_value
        instance.query.side_effect = [
            monty.serialization.loadfn(INSTALL_PATH / "data/Ho-test.json"),
            monty.serialization.loadfn(INSTALL_PATH / "data/B-test.json"),
            monty.serialization.loadfn(INSTALL_PATH / "data/Ho-B-test.json"),
        ]
        yield


@pytest.fixture(scope="module")
def mock_mongo():
    """A fixture for using module-scoped mongomock collections."""
    with mongomock.patch("mongodb://localhost:27017", on_new="create"):
        yield


@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_boxauto(working_folder, mock_mp_rester, mock_mongo):
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
        working_folder=working_folder,
        max_oxy=2,
        standarize_int=True,
        use_preprocessed=True,
        remove_background=True
    )
    r.analyze(n_runs="auto", solver="box", ignore_provider=["AFLOW", "COD"])
    assert r.results.name.iat[0] == ['HoB2', 'HoB4']


@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_boxmethod(working_folder, mock_mp_rester, mock_mongo):
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
        working_folder=working_folder,
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
def test_residualmethod(working_folder, mock_mp_rester, mock_mongo):
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
        working_folder=working_folder,
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
def test_onephaserun(working_folder, mock_mp_rester, mock_mongo):
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
        working_folder=working_folder,
        max_oxy=2,
        standarize_int=True,
        use_preprocessed=True,
        remove_background=True
    )
    r.analyze(n_runs=1, ignore_provider=["AFLOW", "COD"])
    assert r.results.name.iat[0] == 'HoB2'
    shutil.rmtree(r.working_folder)
