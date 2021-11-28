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
import requests
from typing import Tuple
from pymatgen import Structure
import os


class OQMDQuery:
    """
    This class provides an easy way to query and access crystal structures from the OQMD (Open Quantum Materials Database)

    Saal, J. E., Kirklin, S., Aykol, M., Meredig, B., and Wolverton, C.
    "Materials Design and Discovery with High-Throughput Density Functional Theory:
    The Open Quantum Materials Database (OQMD)", JOM 65, 1501-1509 (2013). doi:10.1007/s11837-013-0755-4

    It automatically obtains POSCAR and transforms into CIFs using pymatgen+spglib as backend.

    By default it only downloading structured that has a corresponding entry into the ICSD, to avoid obtaning alot of
    theoretical / prototypes structure.
    Parameters
    ----------
    element_list : Tuple[str]
        An array containing the desired elements to query.
    stability : float, default: 0.5
        Obtains materials with calculated stability (eV/Atom) below threshold. Defaults to 0.5
    limit : int, default: 1000
        Query limit per page. Its kept to a large number as to obtain everything at once.
    """
    def __init__(self, element_list: Tuple[str], dumpfolder: str, stability: float = 0.5, limit: int = 1000):


        self._elements = ",".join(element_list)
        self._element_restriction = len(element_list)
        self._stability = str(stability)
        self._limit = str(limit)
        self._query_url = f"http://oqmd.org/oqmdapi/formationenergy?icsd=True&noduplicate=False&desc=False&sort_offset=0&limit={self._limit}&offset=0&fields=name,entry_id,icsd_id,composition_generic,spacegroup,prototype,ntypes,natoms,volume,delta_e,band_gap,stability&filter=element_set={self._elements}%20AND%20ntypes={self._element_restriction}%20AND%20stability%3C{self._stability}&response_format=json"
        try:
            self._json_data = requests.get(self._query_url).json()
            self._rawdata = self._json_data['data']
        except:
            self._rawdata = None
        self._materialurl = "http://oqmd.org/materials/entry/"
        self._poscarurl = "http://oqmd.org/materials/export/conventional/poscar/"
        self._dumpfolder = dumpfolder

        if not os.path.isdir(self._dumpfolder):
            os.mkdir(self._dumpfolder)

    def get_POSCAR_url(self, entry_id: str) -> str:
        """
        Auxiliary function to obtain the POSCAR from a OQMD id number.
        Scraps the material entry_url to obtain the link to the POSCAR (not in OQMD API)
        Parameters
        ----------
        entry_id : str
            A string representing the entry id of the material in OQMD.

        Returns
        -------
        str
            Returns the link to download the conventional POSCAR from OQMD.
        """

        page_content = requests.get(f"{self._materialurl}{entry_id}").content.decode()
        look_for = page_content.find("/materials/structure/")
        structure_id = page_content[look_for:].split('"')[0].split("/")[-1]
        poscar_link = f"{self._poscarurl}{structure_id}"
        return poscar_link

    def get_structure(self, entry: dict) -> None:
        """
        Get structure from a parsed .JSON dictionary from OQMD.
        Parameters
        ----------
        entry : dict-like
            A pair-key value containing an entry information obtained from the API response of OQMD.

        Returns
        -------
        None
            Returns nothing, save obtained struvures into folder

        """
        entry_id = entry['entry_id']
        poscar_url = self.get_POSCAR_url(entry_id)
        poscar_data = requests.get(poscar_url).content.decode()
        try:
            structure = Structure.from_str(poscar_data, fmt="POSCAR")
            cifname = f"{entry['name']}_OQMD_{entry_id}.cif"
            structure.to(fmt="CIF", filename=os.path.join(self._dumpfolder, cifname), symprec=0.01)
        except:
            print(f"Pymatgen failed to convert OQMD POSCAR to a Structure: {entry['name']}-{entry_id}")


    def query(self) -> None:
        """
        Loop through all the entries obtained from the query and export the structures int othe folder.
        Returns
        -------
        None

        """
        # Loop thru every entry
        if self._rawdata is None:
            return "Failed to Query OQMD."
        for entry in self._rawdata:
            # Save structure into dump folder
            self.get_structure(entry)






