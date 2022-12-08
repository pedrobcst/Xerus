import os
import sys
from pathlib import Path

projectpath = Path(os.getcwd()).parent.as_posix() + os.sep
sys.path.append(projectpath)
import shutil

import pytest
import requests
from Xerus.queriers.optimade import OptimadeQuery

elements_list = ['Ho', 'B']
oqmd_option_1 = {
    "_oqmd_stability": "<0.01"
}
oqmd_option_2 = {
    "_oqmd_stability": "<0.01",
    "nsites": "<=4"
}
cod_option = {
    "_cod_volume": "<10001"
}
COD_URL = "https://www.crystallography.net/cod/optimade/v1"
OQMD_URL = "https://oqmd.org/optimade/v1"


def test_oqmd_onefilter():
    oqmd_query = OptimadeQuery(
        base_url=OQMD_URL,
        elements=elements_list,
        folder_path="optimade_oqmd_1",
        extra_filters=oqmd_option_1
    )
    oqmd_query.query()
    # Should return 4 cifs.
    assert len(os.listdir("optimade_oqmd_1")) == 4
    shutil.rmtree("optimade_oqmd_1")

def test_oqmd_twofilter():
    oqmd_query = OptimadeQuery(
        base_url=OQMD_URL,
        elements=elements_list,
        folder_path="optimade_oqmd_2",
        extra_filters=oqmd_option_2
    )
    oqmd_query.query()
    # Should return 1 cif
    assert len(os.listdir("optimade_oqmd_2")) == 1
    shutil.rmtree("optimade_oqmd_2")

def test_cod_query():
    cod_query = OptimadeQuery(
        base_url=COD_URL,
        elements=elements_list,
        folder_path="optimade_cod",
        extra_filters=cod_option
    )
    cod_query.query()
    assert len(os.listdir("optimade_cod")) == 1
    shutil.rmtree("optimade_cod")


@pytest.mark.skip(reason="Might be time consuming.")
def test_cod_query_large_element_list():
    cod_query = OptimadeQuery(
        base_url=COD_URL,
        elements=['Si', 'O'],
        folder_path="optimade_cod",
        extra_filters=cod_option
    )

def test_oqmd_response_status():
    """
    Test OQMD using HAS ONLY filters
    """
    oqmd_query = OptimadeQuery(
        base_url=OQMD_URL,
        elements=elements_list,
        folder_path="optimade_oqmd_1",
        extra_filters=oqmd_option_1
    )
    headers = oqmd_query.headers
    query_url = f"{oqmd_query.base_url}/{oqmd_query.optimade_endpoint}?{oqmd_query.optimade_filter}&{oqmd_query.optimade_response_fields}&page_limit=10"
    assert requests.get(query_url, headers=headers).status_code == 200
    shutil.rmtree("optimade_oqmd_1")

def test_cod_response_status():
    """
    Test COD using HAS ALL filters
    """
    cod_query = OptimadeQuery(
        base_url=COD_URL,
        elements=elements_list,
        folder_path="optimade_cod",
        extra_filters=cod_option
    )
    headers = cod_query.headers
    query_url = f"{cod_query.base_url}/{cod_query.optimade_endpoint}?{cod_query.optimade_filter_hasall}&{cod_query.optimade_response_fields}&page_limit=10"
    assert requests.get(query_url, headers=headers).status_code == 200
    shutil.rmtree("optimade_cod")

