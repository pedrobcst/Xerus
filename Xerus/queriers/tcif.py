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
#     This file Tests CIFS using GSASII
#     First they are tested to see if GSASII can successfuly simulate patterns.
#     and them a "trial" refinement is ran in "default" experimental pattern. If the refinement is successful, regardless
#     of the Rwp, the CIF is then deemed its OK.
#     This is done to avoid the issue that we cannot control the quality of the CIFS from the databases. Sometimes, one or
#     another file will break GSASII: It will lead to a FloatingPointError that breaks all subsequesent refinements (all
#     subsequent refinements will also give FloatingPointErrors), and it this can only be avoided with a full restart of
#     the system. Thus, to minimize this "problematic" CIFs issues, this test was found to be necessary to maintain stability.
#     Note that this method cannot guarantee that all "bad" structures will be successfuly filtered and also it might lead to some "good"
#     structures to be set as "problematic" (altough I never observed such case), as we are using a "test" pattern that probably differs
#     strongly from most of the tested structures.
#
#     Current issues:
#     - This method is currently quite slow, as it does not run in parallel. To avoid the problem of "FloatingPointErrors",
#     when such Exception is raised (currently we catch all), the script will self restart. To keep a track of what has been tested,
#     the already ran structures information are kept in a json file, and when the problem occurs, this json will be saved with
#     information of all currently tested cifs and their statuses. Thus, when the script restarts, he will started from the
#     next non-tested pattern (that is broken pattern position +1)

import os
import sys
from pathlib import Path
from Xerus.engine.gsas2riet import run_gsas, simulate_pattern
import json
from Xerus.settings.settings import TEST_XRD, INSTR_PARAMS


instr_params = INSTR_PARAMS
powder_data = TEST_XRD
abs_path = Path(__file__).parent
cif_folder = os.path.join(abs_path, "queried_cifs")
cif_name = "cif.json"
datapath = os.path.join(cif_folder,cif_name)
from tqdm import tqdm
from Xerus.utils.cifutils import get_ciflist
if not os.path.exists(datapath):
    get_ciflist(str(cif_folder) + os.sep, "r", True)
with open(datapath, "r") as fp:
    data = json.load(fp)

def load_json(path):
    with open(path, "r") as fp:
        return json.load(fp)

def check_extension(file, allowed=(".json", ".cif", ".ipynb_checkpoints")):
    return any([file.endswith(extension) for extension in allowed])

def save_json(dict, filename=cif_name, folder=cif_folder):
    """

    :param dict:
    :param path:
    :return:
    """
    path = os.path.join(folder, filename)
    with open(path, "w") as fp:
        json.dump(dict, fp)

if __name__ == "__main__":

    for i, entry in tqdm(enumerate(data)):
        filename = entry['filename']
        cifpath = os.path.join(cif_folder, filename)
        ss = entry['simul_status']
        gs = entry['gsas_status']
        if not ss['tested']:
            print("Testing simulation of {}".format(cifpath))
            status = simulate_pattern(filename=cifpath,
                                                   instr_params=instr_params)
            data[i]['simul_status']['tested'] = True
            data[i]['simul_status']['status'] = status
        if not gs['tested']:
            try:
                print("We are going to try to run {} in GSASII".format(cifpath))
                run_gsas(powder_data=powder_data,
                         instr_params=instr_params,
                         phasename="Test.gpx",
                         cif_name=cifpath,
                         max_cyc=1)
            except:
                print("It broke!")
                data[i]['gsas_status']['status'] = -1
                data[i]['gsas_status']['tested'] = True
                data[i]['ran'] = False
                save_json(data)
                # Restart script?
                sys.stdout.flush()
                python = sys.executable
                os.execl(python, python, *sys.argv)
            else:
                print("PASS!")
                data[i]['gsas_status']['status'] = 1
                data[i]['gsas_status']['tested'] = True
                data[i]['ran'] = True
    save_json(data)
    ## clean up .lst, .gpx, .bak
    files = os.listdir(cif_folder)
    for file in files:
        if not check_extension(file):
            os.remove(os.path.join(cif_folder, file))
