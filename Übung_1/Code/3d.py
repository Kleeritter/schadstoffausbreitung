#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# run: python3 xz_plot.py

import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from netCDF4 import Dataset
from mpl_toolkits.mplot3d import Axes3D
#import gc # Garbage Collector


filename = "ubung1.nc"
fileout = "output_file.png"
units = "mg/m^3"
variable = "c"
levels = [0.01,0.02,0.05, 0.1, 0.15, 0.20, 0.25, 0.27,0.28,0.285,0.288, 0.3,0.5]

def nc_read_from_file_2d_all(filename, varname):


   import numpy as np
   import sys
   
   try:
      f = open(filename)
      f.close()
#      print("Load: " + filename + ".")
   except FileNotFoundError:
      print("Error: " + filename + ". No such file. Aborting...")
      sys.exit(1)
   
   nc_file = Dataset(filename, "r", format="NETCDF4")
   tmp_array = np.array(nc_file.variables[varname][:,:,:], dtype=type(nc_file.variables[varname]))

   return tmp_array

def nc_read_from_file_1d_all(filename, varname):


   import numpy as np
   import sys
   
   try:
      f = open(filename)
      f.close()
      print("Load: " + filename + ".")
   except FileNotFoundError:
      print("Error: " + filename + ". No such file. Aborting...")
      sys.exit(1)
   
   nc_file = Dataset(filename, "r", format="NETCDF4")
   tmp_array = np.array(nc_file.variables[varname][:] , dtype=type(nc_file.variables[varname]))


   return tmp_array


#  read data from file
conc = nc_read_from_file_2d_all(filename, variable)

conc = np.where(conc == -9999.0,np.nan,conc)

x = nc_read_from_file_1d_all(filename, "x")
y = nc_read_from_file_1d_all(filename, "y")
z = nc_read_from_file_1d_all(filename, "z")
print(len(x))
x,y,z= np.meshgrid(x,y,z)
print("plotting....")
plt.figure()
fig, ax = plt.subplots()
#CS = ax.tricontour(x,y,conc, levels,colors='k')#ax.contour(x, y, conc, levels, colors='k' )
# alternatively plot 10 levels
CS = ax.scatter(x, y, conc)
#ax.clabel(CS, fontsize=9, inline=1)
plt.title('Concentration (' + units + ')')
plt.xlabel('y (m)')
plt.ylabel('z (m)')
#optional: set y-axis maximum value to 5
#plt.ylim(top=5) 
plt.show()

# print to file
plt.savefig(fileout, format='png', dpi=300)
