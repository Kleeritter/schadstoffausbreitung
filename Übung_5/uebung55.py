import math
import random
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
from tqdm import tqdm
import sympy as sy
from sympy import symbols

##### Einlesen #####
def nc_read_from_file_2d_all(filename, varname):
   f = open(filename)
   f.close()
   nc_file = Dataset(filename, "r", format="NETCDF4")
   tmp_array = np.array(nc_file.variables[varname][:,:], dtype=type(nc_file.variables[varname]))
   return tmp_array

filename = "input_uebung5.nc"

u_wind = np.transpose(nc_read_from_file_2d_all(filename, "u"))
u_wind = np.where(u_wind == -9999.0, np.nan, u_wind)

su = np.transpose(nc_read_from_file_2d_all(filename, "u2"))
su = np.where(su == 0, np.nan, su)

w_wind = np.transpose(nc_read_from_file_2d_all(filename, "w"))
w_wind = np.where(w_wind == -9999.0, np.nan, w_wind)

sw = np.transpose(nc_read_from_file_2d_all(filename, "w2"))
sw = np.where(sw == 0, np.nan, sw)

######## Eingabeparameter ########
n = 10
xq = 61.5
zq = 1
xg = 120
zg = 120
dx = 1
dz = 1
q = 1
lam = 5
k = 0.38
mw = 0.23

gitter_monte = np.zeros((xg, zg))

#####Definition der Reflexion#####
def reflexion(xi, zi, ui, wi, xa, za):

   a = (xi <= 30 and zi <= 60)

   b = (xi >= 90 and zi <= 60)

   c = (30 <= xi <= 90 and zi <= zq)

   while a or b or c:

      r = symbols("r")

      if (xi > xa and za < 60): #rechte Wand
         xi = xi - 2 * (abs(90 - xi)) 
         ui = -ui
         wi = -wi
         #print("rechte Wand")

      elif (xi < xa and za < 60): #linke Wand
         xi = xi + 2 * (abs(30 - xi)) 
         ui = -ui
         wi = -wi
         #print("linke Wand")

      elif (zi < za and za >= 60): #Dach
         zi = zi + 2 * (abs(60 - zi)) 
         ui = -ui
         wi = -wi
         #print("Dach")

      elif (zi <= zq and 30 < xi < 90): #Boden
         zi = zi + 2 * (abs(0 - zi)) 
         ui = -ui
         wi = -wi
         #print("Boden")

      """
      elif sy.solve(r * [1, 1] + [30, 60] - [xa, za], r):
         if xi < 30 and sy.solve(r + [30, 60] - [xi, zi], r):
            xi = xa
            zi = za
            ui = -ui
            wi = -wi
      
      elif sy.solve( r * [1, 1] + [30, 0] - [xa, za], r):
         if xi < 30 and sy.solve(r + [30, 0] - [xi, zi], r):
            xi = xa
            zi = za
            ui = -ui
            wi = -wi
      
      elif sy.solve( -r * [1, 1] + [90, 60] - [xa, za], r):
         if xi > 90 and sy.solve(r + [90, 60] - [xi, zi], r):
            xi = xa
            zi = za
            ui = -ui
            wi = -wi
      
      elif sy.solve( -r * [1, 1] + [90, 0] - [xa, za], r):
         if xi > 90 and sy.solve(r + [90, 0] - [xi, zi], r):
            xi = xa
            zi = za
            ui = -ui
            wi = -wi
      """
   return xi, zi, ui, wi

def positionen(xi, xa, zi, za, ui, wi, su, sw, u_wind, w_wind):
   
   tl = (0.05 * (k * zi) / (1 + k * (zi / lam))) / (mw * math.sqrt(su[int(xi), int(zi)] + sw[int(xi), int(zi)]))
   #print(tl)
   if (0.1 * tl) > (0.05 * (k * 2) / (1 + k * (2 / lam))) / (mw * math.sqrt(su[int(xi), int(zi)] + sw[int(xi), int(zi)])):
      dt = 0.1 * tl
   else:
      dt = (0.05 * (k * 2) / (1 + k * (2 / lam))) / (mw * math.sqrt(su[int(xi), int(zi)] + sw[int(xi), int(zi)]))
      print("alla")
   rl = math.exp(- dt / tl)
   
   if xi == xa:
      ut = rl * ui + math.sqrt(1 - rl **2) * math.sqrt(su[int(xi), int(zi)]) * random.gauss(0, 1) 
   else: 
      ut = rl * ui + math.sqrt(1 - rl **2) * math.sqrt(su[int(xi), int(zi)]) * random.gauss(0, 1)  + ( 1 - rl) * tl * ((su[int(xi), int(zi)] - su[int(xa), int(za)]) / (xi - xa))

   if zi == za:
      wt = rl * wi + math.sqrt(1 - rl **2) * math.sqrt(sw[int(xi), int(zi)]) * random.gauss(0, 1) 
   else:
      wt = rl * wi + math.sqrt(1 - rl **2) * math.sqrt(sw[int(xi), int(zi)]) * random.gauss(0, 1)  + ( 1 - rl) * tl * ((sw[int(xi), int(zi)] - sw[int(xa), int(za)]) / (zi - za))

   ui = u_wind[int(xi), int(zi)] + ut 
   
   wi = w_wind[int(xi), int(zi)] + wt
   
   xi = xi + ui * dt

   zi = zi + wi * dt
   
   xi, zi, ui, wi = reflexion(xi, zi, ui, wi, xa, za)
   #print(xi, zi)
   
   return xi, zi, wi, ui, dt

def gitter(xi, xa, zi, za, dt):

   max_x = int(xi - xa)
   max_z = int(zi - za)
   
   for i in range(0, max_x):
      x_si = math.ceil(xa) + i
      t_k.append((x_si - xa) / (xi - xa))

   for j in range(0, max_z):
      z_si = math.ceil(za) + j
      t_k.append((z_si - za) / (zi - za))
   
   t_ks = sorted(t_k)

   for i in range(1, len(t_ks)):
      
      t_m = np.mean(np.array([t_ks[i-1], t_ks[i]]))
      
      x_gg = xa + t_m * (xi - xa)

      z_gg = za + t_m * (zi - za)
      
      if x_gg > xg or z_gg > zg:
         break
      else:
         
         gitter_monte[int(x_gg), int(z_gg)] += (t_ks[i] - t_ks[i-1]) * dt *((q * dt) / (n * dx * dz))
   return 

for i in tqdm(range(0, n)):
   x = []
   z = [] 
   t_k = []
   xi = xq
   zi = zq
   ui = 0
   wi = 0
   dt = 0

   while ((xi + ui * dt )) < xg:
      print(xi)
      xa = xi
      za = zi

      xi, zi, wi, ui, dt = positionen(xi, xa, zi, za, ui, wi, su, sw, u_wind, w_wind)
      
      if xi > xg:
         xi = xg

      #gitter(xi, xa, zi, za, dt)
      
      x.append(xi)
      z.append(zi)

   plt.plot(x, z)

plt.show()

"""
print('...Erstellen von NC-Datei für Monte-Modell')
fn = 'C:/Users/lisad/Downloads/uebung5.nc'
monte = Dataset(fn, 'w', format='NETCDF4')

monte.createDimension('x', xg)
monte.createDimension('z', zg)

c=monte.createVariable('c', 'f4', ('z', 'x'))
x_vars = monte.createVariable('x', 'f4', 'x')
z_vars =monte.createVariable('z', 'f4', 'z')

c[:]=np.transpose(gitter_monte)
x_vars[:] = np.arange(0.0, xg*dx, dx, dtype=float)
z_vars[:] = np.arange(0.0, zg*dz, dz, dtype=float)

c.units='s/m^2'
x_vars.units='m'
z_vars.units='m'
monte.close()

filename2 = "C:/Users/lisad/Downloads/uebung5.nc"
fileout = "C:/Users/lisad/Downloads/uebung5.png"

def nc_read_from_file_2d_all(filename2, varname):
   f = open(filename2)
   f.close()
   nc_file = Dataset(filename2, "r", format="NETCDF4")
   tmp_array = np.array(nc_file.variables[varname][:,:], dtype=type(nc_file.variables[varname]))
   return tmp_array

def nc_read_from_file_1d_all(filename2, varname):
   f = open(filename2)
   f.close()
   nc_file = Dataset(filename2, "r", format="NETCDF4")
   tmp_array = np.array(nc_file.variables[varname][:] , dtype=type(nc_file.variables[varname]))
   return tmp_array

conc = nc_read_from_file_2d_all(filename2, "c")
conc = np.where(conc == -9999.0,np.nan,conc)

x = nc_read_from_file_1d_all(filename2, "x")
z = nc_read_from_file_1d_all(filename2, "z")

plt.figure()
fig, ax = plt.subplots()

ax.contour(x, z, conc)
units = "mg/$m^3$"
'Concentration (' + units + ')'

plt.legend()
plt.xlabel('z(m)')
plt.ylabel('x(m)')
plt.show()  
"""