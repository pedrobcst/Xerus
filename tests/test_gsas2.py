import sys, os
from pathlib import Path
projectpath = Path(os.getcwd()).parent.as_posix() + os.sep
sys.path.append(projectpath)
from Xerus.settings.settings import TEST_XRD, INSTR_PARAMS, GSAS2_BIN
sys.path.append(GSAS2_BIN)
from Xerus.engine.gsas2riet import simulate_pattern, quick_gsas
import pytest
import shutil


@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_simulation():
    """
    Test GSAS II simulation
    Returns
    -------
    1 if works.
    """

    status = simulate_pattern(filename="cif/HoB2_MP_mp-2267.cif")
    files = [
        "dummy.bak0.gpx",
        "dummy.csv",
        "dummy.gpx",
        "dummy.lst"
    ]
    for file in files:
        os.remove(file)
    assert status == 1

@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_refinement():
    """
    Test GSAS II Quick Refinement
    Returns
    -------

    """
    rwp, gpx = quick_gsas(cif_name="cif/HoB2_MP_mp-2267.cif", powder_data=TEST_XRD, phasename="Any", outfolder=".")
    shutil.rmtree("gsas2_files")
    assert rwp <= 30