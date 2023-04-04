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
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import pymongo
from pymongo.errors import ConnectionFailure
from Xerus.settings.settings import DB_CONN
from Xerus.utils.cifutils import make_system_types, write_cif
from Xerus.utils.tools import create_folder, load_json


class LocalDB:
    """
    This is the main class for handling all database related operations and controlling the local MongoDB
    Its responsible for:
        -> Downloading crystal structures from outside databases and storing locally
        -> Querying and providing to main needed crystal structures
        -> Further operation for consistency such as forced resync (download all CIFs again and check if there is any missing structure)

    Parameters
    ----------
    DB_NAME : str, default: CIF
        MongoDB Database name
    COLLECTION_NAME : str, default: cifs
        Collection name inside the database
    """


    def __init__(self, DB_NAME="CIF", COLLECTION_NAME="cifs"):
        self.client = pymongo.MongoClient(DB_CONN)  # connect to local database
        self.DB_NAME = DB_NAME
        self.COLLECTION_NAME = COLLECTION_NAME
        try:
            self.client.server_info()
            print('Sucessfuly connected to the database')
        except:
            raise ConnectionFailure
        self.database = self.client[self.DB_NAME]
        # If collection does not exist yet, create it and create a unique ID for providers.
        if self.COLLECTION_NAME not in self.database.list_collection_names():
            self.database.create_collection(self.COLLECTION_NAME)
            self.database[self.COLLECTION_NAME].create_index("id", unique=True)
        self.cif_meta = self.database[self.COLLECTION_NAME]
        
        # This will be removed in the future.
        # self.cif_p1 = self.database['AFLOW-P1']

    def check_system(self, system_type: str) -> bool:
        """
        Check for a certain combination of elements in the database.

        Parameters
        ----------
        system_type : str
            `system_type` is a string representation of the chemical system of elements. For example, for Ho and B containing
            materials, `system_type` would be equal to: "Ho-B". For constructing the string in correct manner, please refer
            to Xerus.utils.tools make_system_type function.

        Returns
        -------
        bool
            Returns True if system already present in database, False otherwise.
        """
        query = {"system_type": system_type}
        amount = self.cif_meta.count_documents(query, limit=1)
        if amount > 0:
            return True
        else:
            print("No matching system type {} found in the database".format(system_type))
            return False

    def check_id_p1(self, uid: str) -> bool:
        """
        Check for if a "P1" CIF ID already exists in the P1 database

        Parameters
        ----------
        uid : str
            An AFLOW unique identifier.

        Returns
        -------
        bool
            Returns True if exists, False otherwise.
        """
        query = {"id": uid}
        amount = self.cif_p1.count_documents(query, limit=1)
        if amount > 0:
            return True
        else:
            return False

    def check_id(self, uid: str) -> bool:
        """
        Check for if a CIF ID already exists in the  database

        Parameters
        ----------
        uid : str
            A provider-dependent unique identifier. ie, "mp-xxx" for Materials Project"
        Returns
        -------
        bool
            Returns True if exists, False otherwise.
        """
        query = {"id": uid}
        amount = self.cif_meta.count_documents(query, limit=1)
        return amount > 0



    def upload_many(self, data: Tuple[dict]) -> None:
        """
        Performs a MongoDB `update_many` action to insert new data into the database

        Parameters
        ----------
        data : array-like
            An array-like containing dictionaries with data for upload into the database
        """
        # Ignore error due to duplicate keys (if we are ensuring uniqueness for provider ID).
        try:
            self.cif_meta.insert_many(data, ordered = False)
        except pymongo.errors.BulkWriteError:
            pass



    def check_and_download(self, system_type: str, name : str) -> LocalDB:
        """
        Check for system type in the database. If it is missing, query providers and download the respective CIFs
        pertaining to that system.


        Parameters
        ----------
        system_type : str
            `system_type` is a string representation of the chemical system of elements. For example, for Ho and B containing
            materials, `system_type` would be equal to: "Ho-B". For constructing the string in correct manner, please refer
            to Xerus.utils.tools make_system_type function.
        name: str
            Data set name used for folder creation.

        Returns
        -------
        self
        """
        from Xerus.queriers.multiquery import multiquery
        if not self.check_system(system_type):
            elements = system_type.split("-")
            multiquery(elements, name = name)
        return self

    def check_all(self, system_types: Tuple[str], name: str) -> LocalDB:
        """
        Check for a list of system types and download missing.

        Parameters
        ----------
        system_types : array-like
            Check for a list of `system_type` strings in the database and download the CIFs for respective types if
            not present in the database.


        Returns
        -------
        self
        """
        for combination in system_types:
            # Skip "oxygen" and "CO" cifs.
            if combination == "O" or combination == "C-O":
                pass
            else:
                print("Checking the following combination:{}".format(combination))
                self.check_and_download(combination, name = name)
        return self

    def get_cifs_and_write(self, element_list : List[str], name: str, outfolder: str, maxn: int, max_oxy: int = 2, oxide:bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """

        Automatically receives a list of elements and make `system_types`.
        Check the database for these systems. In case a missing system is present,
        automatically query CIF providers, download the CIFs and upload to the database.

        The CIFs and the metadata / crystal structure information will be automatically exported to `outfolder`,
        and the CIF files will also be written to folder.


        Parameters
        ----------
        element_list : array-like
            An array-like with elements to be queried for. Automatically make `system_types`.
        outfolder : str
            Path to save queried information (CIFs, metadata, etc)
        maxn : int
            System limit. For example if `maxn` = 3, it will grab up to ternary system
        max_oxy : int
            Maximum oxygen for `system_type`. For example `max_oxy` = 2, will allow up to binary oxides, however in case
            of len(element_list) > 2 and O is one the elements, the `system_type`: "A-B-O" will not be considered.
        oxide : bool
            if it is oxide, when ture, only oxide will be in query.

        Returns
        -------
        array-like
            Returns a tuple containing three dataframes
                - cif_meta: CIFs and information that passed through all checks (simulation and refinement).
                            These are the CIFS that will be used for correlation and refinement.
                - cif_notrun: CIFs that failed refinement testing
                - cif_notsim: CIFs that failed simulation testing.
        """
        if not os.path.isdir(outfolder):
            os.mkdir(outfolder)  # make folder if nto exist.
        folder_to_write = 'cifs/'
        final_path = os.path.join(outfolder, folder_to_write)
        queries = make_system_types(element_list, maxn, ceramic_mode=oxide)
        self.check_all(queries, name = name)

        # check oxygen limit
        if max_oxy is not None:
            to_query = []
            for q in queries:
                if q == 'O':
                    continue
                print(q)
                if 'O' in q:
                    if len(q.split("-")) <= max_oxy:
                        to_query.append(q)
                else:
                    to_query.append(q)
            print('Modified Query for given {} to : {}'.format(max_oxy, to_query))
            queries = to_query

        # queries
        query_all_notsimul = {"system_type": {"$in": queries}, "simul_status.status": -1, "num_elem": {"$lte": maxn}}
        query_all_notran = {"system_type": {"$in": queries}, "ran": False, "num_elem": {"$lte": maxn}}
        query_all_ran = {"system_type": {"$in": queries}, "ran": True, "simul_status.status": 1,
                         "gsas_status.status": 1, "num_elem": {"$lte": maxn}}

        # data and meta
        cif_data = pd.DataFrame(self.cif_meta.find(query_all_ran, {"CIF": 1, "filename": 1}))
        cif_meta = pd.DataFrame(self.cif_meta.find(query_all_ran, {"CIF": 0,
                                                                   "ran": 0,
                                                                   "simul_status": 0,
                                                                   "gsas_status": 0, "_id": 0,
                                                                   "pymatgen_status": 0,
                                                                   "num_elem": 0,
                                                                   "elements": 0}))
        cif_meta['full_path'] = cif_meta.filename.apply(lambda filename: os.path.join(final_path, filename))
        cif_not_sim = pd.DataFrame(self.cif_meta.find(query_all_notsimul, {"CIF": 0}))
        cif_not_ran = pd.DataFrame(self.cif_meta.find(query_all_notran, {"CIF": 0}))

        cif_data.apply(lambda row: write_cif(row.filename, row.CIF, final_path), axis=1)  # write cifs to outfolder.
        return (cif_meta, cif_not_ran, cif_not_sim)


    def find_duplicate(self, field):
        """
        Find fields in the database with duplicate value by doing a group aggregation
        :param field: database field.
        :return: pandas dataframe with _id as the field and count as the number of times this field is repeated.
        """
        field = '$' + field

        pipeline = [
            {"$group": {"_id": field, "count": {"$sum": 1}}},  # group by field and sum occurence of field
            {"$match": {"count": {"$gt": 1}}}  # match to values of count > 1.
        ]
        return pd.DataFrame(self.cif_meta.aggregate(pipeline))

    def get_existing_field(self, field: str) -> pd.DataFrame:
        """
        Auxiliary function for getting the existing values of a `field` in the local MongoDB.
        This is used for checking for existence of `system_types` to be provided for the resync_all option.

        Parameters
        ----------
        field : str
            A field name that is actually present in the database. ie: `id`.

        Returns
        -------
        pd.DataFrame
            Returns a pandas dataframe containing a column with all existing `field` values in the database.
        """
        field = '$' + field

        pipeline = [
            {"$group": {"_id": field}}
        ]
        return pd.DataFrame(self.cif_meta.aggregate(pipeline))

    def resync_one(self, system_type: str) -> None:
        """
        Tries to "resync" one `system_type`:
         -> Automatically query all CIF providers for that `system_type`
         -> Check if the CIFs are already present in the database.
            -> In case of no: Test and upload the missing cifs
            -> In case of yes: Does nothing.


        Parameters
        ----------
        system_type : str
            `system_type` is a string representation of the chemical system of elements. For example, for Ho and B containing
            materials, `system_type` would be equal to: "Ho-B". For constructing the string in correct manner, please refer
            to Xerus.utils.tools make_system_type function.

        Returns
        -------
        None
        """
        from Xerus.queriers.multiquery import multiquery
        element_list = system_type.split("-")
        multiquery(element_list , name="resync",resync=True)

    def resync_all(self):
        """
        Resync all `system_types` currently present in the database.
        """
        system_types = self.get_existing_field('system_type')
        system_types['_id'].apply(lambda row: self.resync_one(row))



