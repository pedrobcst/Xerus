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
import json
import os
import shutil
import sys
from pathlib import Path
from typing import List

from Xerus.db.localdb import LocalDB
from Xerus.queriers.cod import CODQuery
from Xerus.queriers.mp import querymp
from Xerus.queriers.optimade import OptimadeQuery
from Xerus.utils.cifutils import (get_provider, make_system, rename_multicif,
                                  standarize)

dbhandler = LocalDB()
abs_path = Path(__file__).parent
proj_path = os.sep.join(str(abs_path).split(os.sep)[1:-1])
dump_folders = ["mp_dump/", "cod_dump/", "aflow_dump/", "oqmd_dump/"]
dump_folders = [os.path.join(abs_path, f) for f in dump_folders]
test_folder = os.path.join(abs_path,"queried_cifs/")


def movecifs(dump_folders: List[str] = dump_folders, test_folder: os.PathLike = test_folder) -> None:
    """
    Helper function to move the queried cifs from each folder into the folder for GSASII testing
    Parameters
    ----------
    dump_folders : A list of each provider "dump" folders.
    test_folder : The folder that will be used to keep the cifs for testing

    Returns
    -------
    None

    """
    if not os.path.isdir(test_folder):
        os.mkdir(test_folder)

    for folder in dump_folders:
        if os.path.isdir(folder):
            files = os.listdir(folder)
            for file in files:
                if not os.path.exists(os.path.join(test_folder, file)):
                    shutil.move(os.path.join(folder, file), test_folder)
            shutil.rmtree(folder)


def load_json(path: os.PathLike) -> dict:
    """
    Helper function to load a JSON file
    Parameters
    ----------
    path : Path to a JSON file

    Returns
    -------
    A dictionary containing the JSON file information
    """

    with open(path, "r") as fp:
        return json.load(fp)


def multiquery(element_list: List[str], max_num_elem: int,  resync:bool = False) -> None:
    """
    Query multiple providers for a given element list and maximum number of elements


    Parameters
    ----------
    element_list : A list of elements to be queried from all providers
    max_num_elem : The maximum number of elements to consider in each combination
    resync : To use resync method or not. If its True, it will redownload all the cifs for those elements and check
    if they already exist in the database or not. Newest CIFS (not current in database) will be them added.
    This defaults to False

    Returns
    -------
    None
    """
    query_type = make_system(standarize(element_list))
    td=[]
    if dbhandler.check_system(query_type) and not resync:
        return "Asked to query {}, but already present in database".format(query_type)

    cod_path = str(os.path.join(abs_path, "cod_dump")) + os.sep
    mp_path = os.path.join(abs_path, "mp_dump")
    aflow_path = os.path.join(abs_path, "aflow_dump")
    oqmd_path = os.path.join(abs_path, "oqmd_dump")
    ## MP query
    querymp(inc_eles=element_list,
             max_num_elem=max_num_elem,
             combination=False, min_hull=0.1, folder_path=mp_path) # Moving hull to 0.1
    td.append(mp_path)

    ## COD query
    cod = CODQuery(elements=element_list,
                   max_num_elements=max_num_elem,
                   combination=False,
                   folder_path=cod_path + os.sep)
    td.append(cod_path)
    cod.query_one(element_list=element_list, rename=True)

    # AFLOW query
    # AFLOWQuery(element_list, outfolder=aflow_path).query()
    # td.append(aflow_path)

    # # OQMD Query
    # print("Querying OQMD...")
    # oqmd = OQMDQuery(element_list=element_list, dumpfolder=oqmd_path)
    # oqmd.query()
    # td.append(oqmd_path)
    

    # TODO: Evaluate on how to implement each querier through optimade via the generic OptimadeQuery interface
    # # Some example OPTIMADE queries
    # print("Querying some OPTIMADE APIs")
    # cod_optimade = OptimadeQuery(
    #     base_url="https://www.crystallography.net/tcod/optimade/v1/",
    #     elements=element_list,
    #     folder_path=Path("optimade")
    # )

    # Currently only MP works.
    # mp_optimade = OptimadeQuery(
    #     base_url="https://optimade.materialsproject.org/v1/",
    #     elements=element_list,
    #     folder_path=Path("optimade")
    # )
    print(f"Querying OQMD through OPTIMADE....")
    oqmd_optimade = OptimadeQuery(
        base_url="https://oqmd.org/optimade/v1",
        elements=element_list,
        folder_path=Path(oqmd_path),
        extra_filters={"_oqmd_stability": "<0.05"}
    )
    oqmd_optimade.query()
    td.append(oqmd_path)

    # td.append(Path("optimade"))

    movecifs(dump_folders=td)
    print("Finished downloading CIFs.")

    if resync:
        print('Multiquery has been ran with resync option.')
        print('Checking for existence of current existing cifs in database and removing...')
        cif_folder = os.path.join(abs_path,"queried_cifs/")
        for filename in os.listdir(cif_folder):
            if not filename.endswith("cif"):
                continue
            uid = get_provider(filename)[1] # get id from filename
            if dbhandler.check_id(uid): # check if id exist
                path = os.path.join(cif_folder, filename)
                print("Removing {} ...".format(path))
                os.remove(path) # remove file if already in db
        ## lets resync here the obtained files to check with the database..
    ## TEST CIFS ##
    print("Testing, uploading and deleting cifs...")
    #cmd = 'python ' +os.sep + str(os.path.join(proj_path,'test_cif.py'))
    cmd = "python " + str(os.path.join(abs_path, "tcif.py"))
    print(cmd)
    os.system(cmd)

    # ## UPDATE DB ##
    print("Uploading database with cifs..")
    data = load_json(os.path.join(abs_path,'queried_cifs','cif.json'))
    print(len(data))
    if len(data) == 0:
        if resync:
            add_dummy = False
            print("No new structure to add.")
        else:
            add_dummy = True
        print('It seems that there is no CIF for the following combination: {}'.format('-'.join(element_list)))
        if add_dummy:
            print('adding dummy entry, ran=False guarantees that there will be no trouble.')
            dummy_entry = {
                "id": "-1",
                "system_type": query_type,
                "gsas_status": {"tested": True, "status": -1},
                "ran": False,
                "type": "dummy"
            }
            dbhandler.upload_many([dummy_entry])
    else:
        dbhandler.upload_many(data)

    ## Remove rest ##
    print("Deleting..")
    shutil.rmtree(os.path.join(abs_path,'queried_cifs'))

