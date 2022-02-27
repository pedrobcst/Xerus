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

# Query AFLOW using AFLUX
import os
import sys
from pathlib import Path
import requests
from urllib.request import urlopen
from pymatgen.io.cif import CifParser, CifWriter
import json
import shutil
from typing import List, Tuple
from Xerus.db.localdb import LocalDB
from Xerus.utils.cifutils import get_ciflist

class AFLOWQuery:
    """
    Queryer for downloading CIFS from AFLOW based on the AFLUX API.
    Obtain "P1" CIFS and symmetrize using pymatgen (spglib backend)
    The P1 cifs are saved into a different database.
    Adapted from the AFLUX Workshop @ AFLOW

    """
    def __init__(self, element_list, outfolder="aflow_dump/"):
        self.element_list = element_list
        self.nspecies = len(element_list)
        self.SERVER="http://aflow.org"
        self.API="/API/aflux/v1.0/?"
        self.MATCHBOOK="species(" + ",".join(self.element_list) +"),nspecies(" + str(self.nspecies) + "),catalog(ICSD),files"
        self.DIRECTIVES="$paging(0)"
        self.outfolder = outfolder
        self.SUMMONS=self.MATCHBOOK+","+self.DIRECTIVES
        self.data = json.loads(urlopen(self.SERVER+self.API+self.SUMMONS).read().decode("utf-8"))
        if not os.path.isdir(outfolder):
            os.mkdir(self.outfolder)
        self.dbconn_p1 = LocalDB().cif_p1
        self.dbhandler = LocalDB()

    def save_and_process(self, url: str, filename: str, outfolder: os.PathLike, symprec: float = 0.01) -> None:
        """
        Process the urls and download the cifs and rename it

        Parameters
        ----------
        url : An URL for a CIF in AFLOW
        filename : Filename from AFLOW
        outfolder : Folder to save CIFS
        symprec : SPGlib symprec parameter for obtaining cifs symmetries.

        Returns
        -------
        None
        """
        tempf = 'p1/'
        path = os.path.join(outfolder, filename)
        tempf_path = os.path.join(outfolder, tempf)
        if not os.path.isdir(tempf_path):
            os.mkdir(tempf_path)
        data = requests.get(url)
        with open(path, "wb") as f:
            f.write(data.content)
        cif_file = CifParser(path)
        try:
            struct = cif_file.get_structures(primitive=False)[0]
        except:
            print("Problem with this CIF, removing..")
            os.remove(path)
            return None
        writer = CifWriter(struct, symprec=symprec)
        comp = struct.composition.reduced_formula
        filename = filename.replace(filename.split("_")[0], comp)
        path_old = path
        path = os.path.join(outfolder, filename)
        p1_path = os.path.join(tempf_path, filename)
        os.rename(path_old, path)
        shutil.move(path, tempf_path)
        p1_meta = get_ciflist(tempf_path + os.sep)[0]
        if not self.dbhandler.check_id_p1(p1_meta['id']):
            # print("Adding raw Aflow P1 entry: {} to P1 database".format(filename))
            self.dbconn_p1.insert(p1_meta)
        os.remove(p1_path)
        path = os.path.join(outfolder, filename)
        writer.write_file(path)

    def make_filename(self, url: str, auid: str) -> str:
        """
        Make a filename using the database standard name: (Material_Provider_ProviderID.cif)
        Parameters
        ----------
        url : An AFLOW URL
        auid : An AFLOW Unique Identifier (aUID)

        Returns
        -------
        A std. filename
        """
        raw_name = url.split("/")[-1]
        code = raw_name.split("_")[-1]
        raw_name = raw_name.replace("ICSD", "AFLOW")
        final_id = auid.split(":")[-1] + "-" + code
        raw_name = raw_name.replace(code, final_id)
        return raw_name

    def getcif(self, response_files: dict) -> List[str]:
        """
        Helper function for getting CIF URLS from the response object.

        Parameters
        ----------
        response_files : A reponse['files'] from a response from aflux api

        Returns
        -------
        A list of filenames
        """
        return [entry for entry in response_files if entry.endswith(".cif")][0]

    def getcif_url(self, data: dict) -> List[Tuple[str, str]]:
        """
        Get the CIF Urls and associated filenames from an AFLUX response dict.

        Parameters
        ----------
        data : A dict coming from AFLUX Response

        Returns
        -------
        A list of tuples contaning the (url, filename)
        """
        urls = ["http://" + entry['aurl'].replace(":", "/") + "/" + self.getcif(entry['files']) for entry in data]
        entry_urls = zip(urls, data)
        filenames = [self.make_filename(url, entry['auid']) for (url, entry) in entry_urls]
        return list(zip(urls, filenames))

    def process_response(self):
        """
        Process all the data, download, rename, get symmetry etc.

        Returns
        -------
        None
        """
        info = self.getcif_url(self.data)
        for url, filename in info:
            self.save_and_process(url, filename, self.outfolder)
        #os.rmdir(os.path.join(self.outfolder, 'p1')) #remove p1 folder

    def query(self):
        """
        For consistency with the COD Query.
        Runs process_response()
        """
        print("Querying AFLOW for following combination {}".format("-".join(self.element_list)))
        self.process_response()
        return self
