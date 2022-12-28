"""This submodule implements the `OptimadeQuery` class, which enables filtering
on any crystal structure database that implements an [OPTIMADE API](https://optimade.org).

"""
import logging
import os
from pathlib import Path
from typing import List, Union

import requests
from Xerus.settings.settings import REQUESTS_TIMEOUT, REQUESTS_HEADER

from optimade.adapters import Structure
from pymatgen.core import Structure as PymatgenStructure


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
        # Create folder for saving structures.
        self.folder_path = Path(folder_path)
        os.makedirs(self.folder_path, exist_ok=True)
        # Set up extra fitlers if needed
        if extra_filters:
            self.extra_filter = " AND ".join([f"{key}{value}" for key, value in extra_filters.items()])
        else:
            self.extra_filter = None
        # Headers to pass through requuests
        self.headers = REQUESTS_HEADER

        self.elements = elements
        self.base_url = base_url
        self.optimade_endpoint = "structures"
        self.symprec = symprec
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

    @property
    def optimade_filter_hasall(self) -> str:
        """Explicit optimade filter string in case that HAS ONLY has not being implemented.

        Returns
        -------
        str
            Returns the optimade string filter with an explicit written elements in case of HAS ONLY it is not implemented.
        """
        filter = "filter=(elements HAS ALL " + ",".join(f'"{e}"' for e in self.elements) + f"AND nelements={len(self.elements)})"
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

        # Construct the query URL first
        if not query_url:
            query_url = f"{self.base_url}/{self.optimade_endpoint}?{self.optimade_filter}&{self.optimade_response_fields}&page_limit=10"

        # Make the request and get the response
        response = requests.get(query_url, headers=self.headers, timeout=REQUESTS_TIMEOUT)
        logging.info("Query %s returned status code %s", query_url, response.status_code)
        next_query_url = None
        
        # First handle the issues with the response
        if response.status_code == 404:
            raise ValueError("Query returned 404, check provider URL")

        
        if response.status_code == 403:
            raise ConnectionRefusedError("Query returned 403. Check requests settings.")
        
        if response.status_code == 501:
            # If the query returns 501 Not Implemented, assume it is the HAS ONLY which failed and try again with an explicit filter
            query_url = f"{self.base_url}/{self.optimade_endpoint}?{self.optimade_filter_hasall}&{self.optimade_response_fields}&page_limit=10"
            # print(f"Retrying query with {query_url} ....")
            response = requests.get(query_url, headers=self.headers, timeout=REQUESTS_TIMEOUT)

        if response.status_code == 504:
            raise ValueError("Query returned 504, try changing query filters..")
        
        if response.status_code == 200:
            data = response.json()
            meta = data["meta"]
            if meta["more_data_available"]:
                next_query_url = data["links"]["next"]
                # Check if what is returned is dict or straight the URL
                if isinstance(next_query_url, dict):
                    next_query_url = next_query_url['href']

        # Only if query was successful we parse the data.
            for entry in data["data"]:
                try:
                    structure = Structure(entry)
                    # Get the suffix from provider and provider id
                    cif_suffix = self.make_suffix(entry=structure.entry.dict(), meta=meta)
                    # Convert the optimade structure into pymatgen format
                    pymatgen_structure = PymatgenStructure.from_str(structure.convert("cif"), fmt="cif")
                    # Get the reduced formula
                    formula = pymatgen_structure.composition.reduced_formula
                    # Make the cifname
                    cifname = f"{formula}_{cif_suffix}"
                    print(f"Saving {cifname}... to {self.folder_path}/{cifname}...")
                    # Save cif to path
                    pymatgen_structure.to(fmt="cif", filename=self.folder_path.joinpath(cifname), symprec=self.symprec)
                except (ValueError, TypeError) as e:
                        print(f'Failed to convert {entry["id"]} to pymatgen structure..')
        else:
            print("Something else went wrong")
            raise ValueError(f'Query return {response.status_code}.....')

        if next_query_url:
            self.query(next_query_url)
