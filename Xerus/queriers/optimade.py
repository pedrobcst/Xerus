"""This submodule implements the `OptimadeQuery` class, which enables filtering 
on any crystal structure database that implements an [OPTIMADE API](https://optimade.org).

"""
import logging
import os
import sys
from pathlib import Path
from typing import List, Union

project_path = str(Path(os.path.dirname(os.path.realpath(__file__))).parent) + os.sep # so convoluted..
if project_path not in sys.path:
    sys.path.append(project_path)

import requests
from Xerus.utils.tools import create_folder

from optimade.adapters import Structure


class OptimadeQuery:

    def __init__(
        self, 
        base_url: str, 
        elements: List[str] = None, 
        folder_path: os.PathLike = "",
        symprec: float = 0.01,
        
    ):
        """Initialise the query objects for a given database.
        
        Parameters:
            elements: The list of element symbols that define the chemical space to query.
            base_url: The base URL of the OPTIMADE API for the database.
            symprec: Symmetry tolerance to pass to spglib for symmetrization purposes (default = 0.01)

        """

        self.elements = elements
        self.folder_path = Path(folder_path)
        self.base_url = base_url
        self.optimade_endpoint = "structures"
        self.optimade_filter = "filter=elements HAS ONLY " + ",".join(f'"{e}"' for e in self.elements)
        self.symprec = symprec
        os.makedirs(self.folder_path, exist_ok=True)
        self.optimade_response_fields = "response_fields=cartesian_site_positions,species,elements,nelements,species_at_sites,lattice_vectors,last_modified,elements_ratios,chemical_formula_descriptive,chemical_formula_reduced,chemical_formula_anonymous,nperiodic_dimensions,nsites,structure_features"

    def make_suffix(self, entry: dict, meta: dict) -> str:
        """Makes CIF suffix name from an OPTIMADE entry dictionary and meta-data information

        Parameters
        ----------
        entry : dict
            OptimadeStructure.struc.dict()
        meta : dict
            Meta-data dictionary

        Returns
        -------
        str
            Returns the suffix of _Provider_ProviderID
        """
        return f'{meta["provider"]["prefix"].upper()}_{entry["id"]}.cif'

    def query(self, query_url: Union[str, None] = None) -> None:

        if not query_url:
            query_url = f"{self.base_url}/{self.optimade_endpoint}?{self.optimade_filter}&{self.optimade_response_fields}"

        response = requests.get(query_url)
        logging.info("Query %s returned status code %s", query_url, response.status_code)
        print(f"Query {query_url} returned status code {response.status_code}")
        next_query_url = None
        if response.status_code == 200:
            data = response.json()
            meta = data["meta"]
            if meta["more_data_available"]:
                next_query_url = data["links"]["next"]

            structures = [Structure(entry) for entry in data["data"]]

            for structure in structures:
                # Get the suffix from provider and provider id
                cif_suffix = self.make_suffix(entry=structure.entry.dict(), meta=meta)
                # Convert the optimade structure into pymatgen format
                pymatgen_structure = structure.convert("pymatgen")
                # Get the reduced formula
                formula = pymatgen_structure.composition.reduced_formula
                # Make the cifname
                cifname = f"{formula}_{cif_suffix}"
                print(f"Saving {cifname}... to {self.folder_path}/{cifname}...")
                # Save cif to path
                pymatgen_structure.to(fmt="cif", filename=self.folder_path.joinpath(cifname), symprec=self.symprec)

        if next_query_url:
            self.query(next_query_url)
