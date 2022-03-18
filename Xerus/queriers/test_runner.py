# Pass few of the responsabilities from the tcif to here
# This file should be responsbile for the following:
# - parse the cifs
# - break the parsed data into smaller chunks
# - run in parallel for each chunk
# - aggregate the results into a final 'cif.json'

import json
import os
import subprocess
import sys
from pathlib import Path

from Xerus.utils.cifutils import get_ciflist
from Xerus.utils.tools import chunk_list

abs_path = Path(__file__).parent
cif_folder = os.path.join(abs_path, "queried_cifs")

# This should first break up the cifs into smaller chunks
data = get_ciflist(str(cif_folder) + os.sep)
chunks = chunk_list(data, os.cpu_count())

# Save chunks
process_chunks = []
for i, chunk in enumerate(chunks):
    filename = "cif_chunk_{}.json".format(i)
    path = os.path.join(cif_folder, filename)
    # with open(path, "w") as fp:
    #     json.dump(chunk, fp)
    process_chunks.append(filename)


if __name__ == "__main__":
    processes = []
    print('Starting tests...')
    # for chunk_name in process_chunks:
    #     command = ["python", os.path.join(abs_path,"tcif_parallel.py"), os.path.join("queried_cifs", chunk_name)]
    #     process = subprocess.Popen(command)
    #     processes.append(process)
    # for proc in processes:
    #     proc.communicate()
    print('Finished. Preparing to merge chunks into a single file.')
    data = []
    for chunk_name in process_chunks:
        with open(os.path.join(cif_folder, chunk_name), "r") as fp:
            chunk = json.load(fp)
            data.extend(chunk)
    # Save merged file
    with open(os.path.join(cif_folder, "cif.json"), "w") as fp:
        json.dump(data, fp)
    print(f'Merged {os.cpu_count()} chunks into a single file of {len(data)} entries.')




