import os
import sys
from pathlib import Path

from Xerus.settings.settings import GSAS2_BIN, INSTR_PARAMS, TEST_XRD
from . import INSTALL_PATH

sys.path.append(GSAS2_BIN)
import shutil

import pytest
from Xerus.engine.gsas2riet import quick_gsas, simulate_pattern


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

    status = simulate_pattern(filename=INSTALL_PATH / "cif/HoB2_MP_mp-2267.cif")
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
    rwp, gpx = quick_gsas(cif_name=INSTALL_PATH / "cif/HoB2_MP_mp-2267.cif", powder_data=TEST_XRD, phasename="Any", outfolder=".")
    shutil.rmtree("gsas2_files")
    assert rwp <= 35
