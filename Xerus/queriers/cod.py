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
from itertools import combinations
from pathlib import Path
from typing import List, Tuple

import requests

from Xerus.utils.cifutils import rename_multicif
from Xerus.settings.settings import REQUESTS_TIMEOUT


def make_name(url: str):
    """
    Auxiliary function to retrieve CIF name from an URL

    Parameters
    ----------
    url : An URL for a CIF in COD.

    Returns
    -------
    A CIFname string from an URL from COD
    """
    return url.split("/")[-1].split(".")[0] + ".cif"


class CODQuery(object):
    """
    Query COD
    """

    def __init__(self, elements: List = None, max_num_elements: int = None, min_num_elements: int = 1,
                 file_format: str = "urls", combination: bool = True,
                 folder_path: os.PathLike = ""):

        self.elements = elements
        self.maxn = max_num_elements
        self.minn = min_num_elements
        self.file_format = file_format
        self.combination = combination
        self.folder_path = folder_path
        self.base_url = "http://www.crystallography.net/cod/result?"
        if not os.path.isdir(self.folder_path):
            os.mkdir(self.folder_path)


    def query_one(self, element_list: Tuple[str], rename: bool = False) -> None:
        """
        Query one list of elements in the COD and download the associated CIFS to specified folder.
        Parameters
        ----------
        element_list : A list of elements to be queried, ie ['Ho', 'B'] -> Ho & B binary elements
        rename : To automatically rename the CIFs using the database covention (Formula_Provider_providerID.CIF)

        Returns
        -------

        """
        element_string = "&".join(['el' + str(i + 1) + '=' + element_list[i] for i in range(0, len(element_list))])
        strictmax = len(element_list)
        strictmin = self.minn
        # query_url = self.base_url + element_string + "&strictmin=" + str(strictmin) + "&strictmax=" + str(
            # strictmax) + "&format=" + str(self.file_format)
        query_url = f"{self.base_url}{element_string}&strictmin={strictmin}&strictmax={strictmax}&maxv=10000&format={self.file_format}"

        if element_string == "el1=O" or element_string == "el1=C&el2=O":
            print("Oxygen only or Carbon and Oxygen only found. Skipping this query")
            return -1
        else:
            print("I am querying the following combination from COD: {}".format("-".join(element_list)))
            print(query_url)
            cif_data = requests.get(query_url,  timeout=REQUESTS_TIMEOUT)
            if cif_data.status_code == 200:
                cif_urls = cif_data.text.split("\n")[:-1]
                for url in cif_urls:
                    cif = requests.get(url, timeout=REQUESTS_TIMEOUT)
                    if cif.status_code == 200:
                        name = make_name(url)
                        with open(os.path.join(self.folder_path, name), "wb") as f:
                            f.write(cif.content)
                    else:
                        print("Status {}".format(cif.status_code))
            else:
                print("Status: {}".format(cif_data.status_code))
        if rename:
            rename_multicif(self.folder_path, provider="COD")

    def query(self) -> None:
        """
        Query the COD in the following fashion for each combination of elements ranging from min to max:
        "http://www.crystallography.net/cod/result?el1=Cu&el2=B&strictmin=1&strictmax=3&format=urls"

        Returns
        -------
        None
        """
        if self.combination:
            comb = [list(combinations(self.elements, x)) for x in range(1, self.maxn + 1)]
            combinations_flat = [item for sublist in comb for item in sublist]
            for tup in combinations_flat:
                self.query_one(tup, rename=False)
            rename_multicif(self.folder_path, provider="COD")
        else:
            self.query_one(self.elements, rename=True)
        return self
