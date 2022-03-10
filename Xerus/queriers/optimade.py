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
        extra_filters: dict = None
        
    ):
        """Initialise the query objects for a given database.
        
        Parameters:
            elements: The list of element symbols that define the chemical space to query.
            base_url: The base URL of the OPTIMADE API for the database.
            symprec: Symmetry tolerance to pass to spglib for symmetrization purposes (default = 0.01)
            extra_filters: extra paramaters to pass to filters. For example, if you want to query for a specific stability, you can pass {"stability": "{condition}{value}"} ie: {"stability": ">=0.5"}.

        """
        # TODO: Add option to pass a dictionary of extra filters to the query that can be done based on provider.
        if extra_filters:
            self.extra_filter = " AND ".join([f"{key}{value}" for key, value in extra_filters.items()])
        else:
            self.extra_filter = None
        self.elements = elements
        self.folder_path = Path(folder_path)
        self.base_url = base_url
        self.optimade_endpoint = "structures"
        self.symprec = symprec
        os.makedirs(self.folder_path, exist_ok=True)
        self.optimade_response_fields = "response_fields=cartesian_site_positions,species,elements,nelements,species_at_sites,lattice_vectors,last_modified,elements_ratios,chemical_formula_descriptive,chemical_formula_reduced,chemical_formula_anonymous,nperiodic_dimensions,nsites,structure_features,dimension_types"

    @property
    def optimade_filter(self) -> str:
        """Optimade filter string for the query.

        Returns
        -------
        str
            Returns the string of the optimade filter and take care with there are any extra filter or not.
        """
        filter = "filter=elements HAS ONLY " + ",".join(f'"{e}"' for e in self.elements)
        if self.extra_filter is not None:
            filter += f" AND {self.extra_filter}"
        return filter

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
        print("Querying......")
        if not query_url:
            query_url = f"{self.base_url}/{self.optimade_endpoint}?{self.optimade_filter}&{self.optimade_response_fields}&page_limit=10"

        response = requests.get(query_url)
        logging.info("Query %s returned status code %s", query_url, response.status_code)
        print(f"Query {query_url} returned status code {response.status_code}")
        next_query_url = None
        if response.status_code == 404:
            raise ValueError("Query returned 404, check provider URL")
        if response.status_code == 200:
            data = response.json()
            meta = data["meta"]
            if meta["more_data_available"]:
                next_query_url = data["links"]["next"]


            # Lets move this to try catch block
            # structures = [Structure(entry) for entry in data["data"]]

            for entry in data["data"]:
                try:
                    structure = Structure(entry)
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
                except ValueError:
                        print(f'Failed to convert {entry["id"]} to pymatgen structure..')
                
        if next_query_url:
            self.query(next_query_url)
