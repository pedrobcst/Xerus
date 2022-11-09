# Test for configuration & GSAS Binaries
import sys, os
from pathlib import Path
from . import INSTALL_PATH
from Xerus.settings.settings import TEST_XRD, INSTR_PARAMS, GSAS2_BIN
sys.path.append(GSAS2_BIN)
import GSASIIscriptable as G2sc
from Xerus.db.localdb import LocalDB
from Xerus.settings.settings import MP_API_KEY
# from pymatgen.ext.matproj import MPRestError
# from pymatgen.ext.matproj import MPRester
from mp_api.client import MPRestError
from mp_api.client import MPRester
MP_API_KEY_WRONG = "bQEFQ!"
import pytest

@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.user
def test_dbconn():
    """
    Tests connection to localdatabase
    The purpose is to check if local database is correctly configurated
    Returns
    -------

    """
    def connect():
        client = LocalDB().client
        try:
            info = client.server_info()
            return True
        except:
            return False
    assert connect() == True, "CIF Database connection failed"

@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.user
def test_mpconn():
    """
    Tests Materials Project Connection
    Purpose is to test mainly API Key is correctly set.
    Returns
    -------
    True if OK, False or Not
    """
    try:
        # MPRester(MP_API_KEY).get_data("HoB2")
        MPRester(MP_API_KEY).summary.search(formula = "HoB2", fields=["material_id", "band_gap"])
        assert True
    except MPRestError:
        assert False, "Failed to connect to MP Project"


@pytest.mark.filterwarnings('ignore::RuntimeWarning')
@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_gsas2settings():
    """
    Test configuration for GSASII
    Returns
    -------
    True or False
    """
    try:
        import sys, os

        # Create a gpx project
        filename="test.gpx"
        gpx = G2sc.G2Project(filename=filename)
        # Add histogram
        histogram = gpx.add_powder_histogram(datafile=TEST_XRD, iparams=INSTR_PARAMS)
        # Add test cif
        gpx.add_phase(phasefile=INSTALL_PATH / "cif/HoB2_MP_mp-2267.cif", phasename="HoB2 test", histograms=[histogram])

        # Refine:
        refdict0 = {"set": {"Background": {"no. coeffs": 6, "refine": True},
                        "Cell": True,
                        "Instrument Parameters": ["Zero"],
                        "Scale": True}}
        refdict4a = {"set": {'Sample Parameters': ['Shift', 'DisplaceX', 'DisplaceY', 'Scale']}}
        refdict5c = {"set": {'Instrument Parameters': ['X', 'Y']}}
        dictList = [refdict0, refdict4a, refdict5c]
        gpx.do_refinements(dictList)
        # Check rwp
        rwp = gpx.histograms()[0].get_wR()
        os.remove(filename)
        os.remove(filename.replace(".gpx", ".lst"))
        os.remove(filename.replace(".gpx", ".bak0.gpx"))

        if rwp <= 50:
            assert True, "RWP lower than 50%"
        else:
            assert False, "Large rwp."
    except:
        assert False, "GSAS II test failed."


