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
import os
import urllib.parse
from configparser import ConfigParser
from pathlib import Path



abs_path = Path(__file__).parent
project_path = str(Path(os.path.dirname(os.path.realpath(__file__))).parent) + os.sep
__version__ = '1.1b'

config = ConfigParser()
config.read(abs_path.joinpath("config.conf"))

## Get the settings entrys ##
INSTR_PARAMS = os.path.join(project_path, config['gsas2']['instr_params'])
# GSAS2_BIN = config['gsas2']['binpath']
GSAS2_BIN = os.path.join(project_path, 'GSASII') # bundled GSAS II
TEST_XRD = os.path.join(project_path,config['gsas2']['testxrd'])
MP_API_KEY = config['mp']['apikey']

REQUESTS_TIMEOUT = 60.0
REQUESTS_HEADER =  {"User-Agent": f"Xerus/{__version__}"}


if config['mongodb']['host'] == 'localhost':
    DB_CONN = 'localhost'
else:
    host=config['mongodb']['host']
    user=urllib.parse.quote(config['mongodb']['user'])
    password=config['mongodb']['password']
    DB_CONN = "mongodb://%s:%s@%s" % (user, password, host)
