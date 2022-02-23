# -*- coding: utf-8 -*-
"""
Created on Sat May 23 16:51:32 2020
This script read the HOS_Convert debug file and plot the results to validate the results.
To output the debug file, set iDebug=1 in HOS_Converter program, the output file will be:
    1. Debug_SF_??????.dat: free surface elevation data and velocity.
    2. Debug_LF_??????.dat: left BC velocity data.
    3. Debug_RF_??????.dat: right BC velocity data.

@author: XHH
"""
## In[2]:
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

data_folder = r'E:\P_FPSO_Reproducible CFD JIP\P_FPSO_Reproducible_CFD_JIP_NWT\HOS2UDW2P_Converter\HOS2UDW2P_Converter_V20200717\debug_data'

files = [img for img in os.listdir(data_folder) if "Debug_SF_" in img]  

nfiles =len(files) #get the number of files

#get the image modification time and sort the files based on the modification time.
#Note: the file creation time may chage due to the copy and paste
files_time =np.zeros(nfiles)

for i in range(nfiles):
    fname = os.path.join(data_folder, files[i])
    files_time[i]=os.path.getmtime(fname)
        
#get the index of file nmae
file_index = np.argsort(files_time)

#get the file file
fname = os.path.join(data_folder, files[file_index[0]])
#load the file 
fdata = pd.read_csv(fname, header=None, index_col=None,names=["x","z","u","v","w"], delimiter=r"\s+")

plt.plot(fdata["x"],fdata["z"],label="Elev");
#plt.plot(fdata["x"],fdata["u"],label="u");
#plt.plot(fdata["x"],fdata["w"],label="w");
