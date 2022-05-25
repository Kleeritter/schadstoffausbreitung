#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from netCDF4 import Dataset




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
   tmp_array = np.array(nc_file.variables[varname][:,:], dtype=type(nc_file.variables[varname]))

   return tmp_array

def nc_read_from_file_1d_all(filename, varname):


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
   tmp_array = np.array(nc_file.variables[varname][:] , dtype=type(nc_file.variables[varname]))


   return tmp_array



filename = "file1.nc"
filename2 = "file1.nc"

units = "g/m$^3$"

conc = nc_read_from_file_2d_all(filename, "conc")
conc2 = nc_read_from_file_2d_all(filename2, "conc")

conc = np.where(conc == -9999.0,np.nan,conc)
x = nc_read_from_file_1d_all(filename, "x")
z = nc_read_from_file_1d_all(filename, "z")

conc2 = np.where(conc2 == -9999.0,np.nan,conc2)
x2 = nc_read_from_file_1d_all(filename2, "x")
z2 = nc_read_from_file_1d_all(filename2, "z")


print("plotting....")
plt.figure()

fig, ax = plt.subplots()

levels = [0.05, 0.1, 0.15, 0.20, 0.25]
CS = ax.contour(x, z, conc, levels, colors='black' )
CS2 = ax.contour(x2, z2, conc2, levels, colors='red' )

ax.clabel(CS, fontsize=9, inline=1)

plt.title('Concentration (' + units + ')')
plt.xlabel('x (m)')
plt.ylabel('z (m)')
plt.show()



#cbar = plt.colorbar(CS2,ax=axs[1], shrink=0.6)

#cbar.set_clim(12.0, 18.0)
#cbar.ax.set_ylabel('Lufttemperatur (°C)')


plt.savefig(fileout, format='png', dpi=300)
