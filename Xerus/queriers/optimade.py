"""This submodule implements the `OptimadeQuery` class, which enables filtering 
on any crystal structure database that implements an [OPTIMADE API](https://optimade.org).

"""
import logging
import os
import sys
from pathlib import Path
from typing import List
project_path = str(Path(os.path.dirname(os.path.realpath(__file__))).parent) + os.sep # so convoluted..
if project_path not in sys.path:
    sys.path.append(project_path)

import requests
from optimade.adapters import Structure, get_cif


class OptimadeQuery(object):

    def __init__(
        self, 
        base_url: str, 
        elements: List[str] = None, 
        folder_path: os.PathLike = ""
    ):
        """Initialise the query objects for a given database.
        
        Parameters:
            elements: The list of element symbols that define the chemical space to query.
            base_url: The base URL of the OPTIMADE API for the database.

        """

        self.elements = elements
        self.folder_path = Path(folder_path)
        self.base_url = base_url
        self.optimade_endpoint = "structures"
        self.optimade_filter = "filter=elements HAS ONLY " + ",".join(f'"{e}"' for e in self.elements)
        os.mkdir(self.folder_path, exist_ok=True)

    def make_name(self, entry, meta) -> str:
        return f'{meta["provider"]["prefix"]}_{entry["id"]}.cif'


    def query(self, query_url: str | None) -> None:

        if not query_url:
            query_url = f"{self.base_url}/{self.optimade_endpoint}?{self.optimade_filter}"

        response = requests.get(query_url)
        logging.info("Query %s returned status code %s", query_url, response.status_code)
        next_query_url = None
        if response.status_code == 200:
            data = response.json()
            meta = data["meta"]
            if meta["more_data_available"]:
                next_query_url = data["links"]["next"]

            structures = [Structure(**entry) for entry in data["data"]]

            for structure in structures:
                cif_fname = self.make_name(structure.entry, structure)
                with open(self.folder_path / cif_fname, "w") as f:
                    f.write(get_cif(structure))

        if next_query_url:
            self.query(next_query_url)