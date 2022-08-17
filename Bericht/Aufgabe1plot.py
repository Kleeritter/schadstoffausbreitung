import math as m
import random
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset



filename = "Bericht/b1.nc"
units = "g/m^3"

def nc_read_from_file_2d_all(filename, varname):
   nc_file = Dataset(filename, "r", format="NETCDF4")
   tmp_array = np.array(nc_file.variables[varname][:,:],
   dtype=type(nc_file.variables[varname]))
   return tmp_array

def nc_read_from_file_1d_all(filename, varname):
   nc_file = Dataset(filename, "r", format="NETCDF4")
   tmp_array = np.array(nc_file.variables[varname][:],
   dtype=type(nc_file.variables[varname]))
   return tmp_array

def nc_read_from_file_3d(filename, varname, z):
   nc_file = Dataset(filename, "r", format="NETCDF4")
   tmp_array = np.array(nc_file.variables[varname][:,0],
   dtype=type(nc_file.variables[varname]))
   return tmp_array

conc = nc_read_from_file_2d_all(filename, "c")
x = nc_read_from_file_1d_all(filename, "x")
z = nc_read_from_file_1d_all(filename, "z")



plt.figure()
fig, ax1 = plt.subplots()
levels = [0.1, 0.20, 0.3, 0.4, 0.5]
CS = ax1.contour(x, z, conc, levels, colors='black' )
ax1.clabel(CS, fontsize=9, inline=1)
plt.title('Konzentration (' + units + ')')
plt.xlabel('x (m)')
plt.ylabel('z (m)')
plt.savefig('Bericht/aufgabe_1a.png', format='png', dpi=300)

plt.figure()
fig, ax2 = plt.subplots()
concconc = nc_read_from_file_3d(filename, "c", z)
ax2.plot(x, concconc)
plt.title('Konzentration (' + units + ') am Erdboden')
plt.ylabel('c (' + units + ')')
plt.xlabel('x (m)')
plt.savefig('Bericht/aufgabe_1a_plot2.png', format='png', dpi=300)
