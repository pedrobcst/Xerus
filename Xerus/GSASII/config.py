# -*- coding: utf-8 -*-
'''
*config.py: Configuration options*
----------------------------------
This file created in SelectConfigSetting on 05 02 2021 14:41
'''

import os.path
import GSASIIpath

Main_Pos = (2132, 48)
'''Main window location - will be updated & saved when user moves
it. If position is outside screen then it will be repositioned to default
'''

Main_Size = (1138, 935)
'''Main window size (width, height) - initially uses wx.DefaultSize but will updated
 and saved as the user changes the window
'''

Plot_Pos = (95, 129)
'''Plot window location - will be updated & saved when user moves it
these widows. If position is outside screen then it will be repositioned to default
'''

Plot_Size = (778, 655)
'''Plot window size (width, height) - initially uses wx.DefaultSize but will updated
 and saved as the user changes the window
'''

previous_GPX_files = [
	  "/home/pedrobcst/Dropbox/PycharmProjects/SSXRA/tests/Gd3Ni6Si2_ann/optuna/20210205_Gd3Ni6Si2_Pedro_ann_850C_1week_seed42_trial_378.gpx",
	  "/home/pedrobcst",
	  "/home/pedrobcst/Dropbox/PycharmProjects/pythonProject/test1/hob2/HoB10-Atom.gpx",
	  "/home/pedrobcst/Dropbox/PycharmProjects/pythonProject/simul_test/PbSO4sim_ok.gpx",
   ]
'''A list of previously used .gpx files
'''

