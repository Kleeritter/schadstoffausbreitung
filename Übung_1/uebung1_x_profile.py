#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from netCDF4 import Dataset




def nc_read_from_file_3d(filename, varname, ind_z, ind_y):


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
   tmp_array = np.array(nc_file.variables[varname][ind_z,ind_y,:], dtype=type(nc_file.variables[varname]))

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



filename = "uebung1_result.nc"
fileout = "uebung1.png"
units = "mg/m$^3$"

conc = nc_read_from_file_3d(filename, "conc", 0, 100)

conc = np.where(conc == -9999.0,np.nan,conc)
x = nc_read_from_file_1d_all(filename, "x")


print("plotting....")
fig = plt.figure()
ax = plt.axes()

ax.plot(x, conc)

plt.ylabel('c (' + units + ')')
plt.xlabel('x (m)')
plt.ylim(top=0.5)
plt.show()



#cbar = plt.colorbar(CS2,ax=axs[1], shrink=0.6)

#cbar.set_clim(12.0, 18.0)
#cbar.ax.set_ylabel('Lufttemperatur (°C)')


plt.savefig(fileout, format='png', dpi=300)
