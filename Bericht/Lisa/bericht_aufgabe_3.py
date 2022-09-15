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

filename = "Bericht/input_uebung5.nc"

u = np.transpose(nc_read_from_file_2d_all(filename, "u"))
u = np.where(u == -9999.0, np.nan, u)

su = np.transpose(nc_read_from_file_2d_all(filename, "u2"))

w = np.transpose(nc_read_from_file_2d_all(filename, "w"))
w = np.where(w == -9999.0, np.nan, w)

sw = np.transpose(nc_read_from_file_2d_all(filename, "w2"))

######## Eingabeparameter ########
# Anzahl Partikel
n = 10^2
# Quellort (m)
xq = 60.5
zq = 0.5
# Modellgrenzen (m)
nx = 120
nz = 120
# Gitterweite (m)
dx = 1
dz = 1
# Quellstärke (g/m^3)
q = 1
# Karman-Konstante
kappa = 0.38

#Gitter-Array
gitter_monte = np.zeros((nx, nz))

#####Definition der Reflexion#####
def reflexion(xi, zi, ui, wi, xa, za):

   a = (xi <= 31 and zi <= 61) #rechtes Gebäude

   b = (xi >= 90 and zi <= 61) #linkes Gebäude

   c = (31 <= xi  and xi<= 91 and zi < 1)

   while  (xi <= 31 and zi <= 61) or  (xi >= 90 and zi <= 61) or (31 <= xi  and xi<= 91 and zi < 1):
     # print(31 <= xi  and xi<= 91 and zi < 1, c)
      if (xi >= 90 and zq < zi <= 61): #rechte Wand
         xi = xi - 2 * (abs(90 - xi)) 
         ui = -ui
         print("rechte Wand")

      elif (xi <= 31 and zq < zi <= 61): #linke Wand
         xi = xi + 2 * (abs(31 - xi)) 
         ui = -ui
         print("linke Wand")

      elif (zi <= 61 and xi < 31 or xi > 90): #Dach
         zi = zi + 2 * (abs(61 - zi)) 
         wi = -wi
         print("Dach")

      elif (31 <= xi and xi <= 90 and zi < 1): #Boden
         zi = zi + 2 * (abs(1 - zi)) 

         wi = -wi
         print("Boden")
      else :
         xi=xi
         zi=zi
         wi=wi
         ui=ui
         #print("wels ", c)

   return xi, zi, ui, wi

def positionen(xi, xa, zi, za, ui, wi, su, sw, u, w):
   zii= int(zi)
   zaa=int(za)
   zbb=int(zb)
   if zii <1:
      zii=1
   if zaa <1:
      zaa=1
   if zbb <1:
      zbb=1
   if su[int(xi), int(zii)] == 0 or sw[int(xi), int(zii)] == 0: # Wenn Standardabweichung Null, dann wird auch tl Null!
      print("alla ", xi , zi)
      tl = 0.0001
      dt = 0.1 * tl
   else:
      tl = (0.05 * (kappa * zi) / ((1 + kappa * (zi / 5))) / (0.23 * math.sqrt(su[int(xi), int(zii)] + sw[int(xi), int(zii)]))) 

      if (0.1 * tl) > (0.05 * (kappa * 2) / (1 + kappa * (2 / 5))) / (0.23 * math.sqrt(su[int(xi), int(zii)] + sw[int(xi), int(zii)])):
         dt = 0.1 * tl
      else:
         dt = (0.05 * (kappa * 2) / (1 + kappa * (2 / 5))) / (0.23 * math.sqrt(su[int(xi), int(zii)] + sw[int(xi), int(zii)]))
   
   rl = math.exp(- dt / tl)

   ut = rl * ui + math.sqrt(1 - rl **2) * math.sqrt(su[int(xi), int(zii)]) * random.gauss(0, 1) + ( 1 - rl) * tl * ((su[int(xa), int(zaa)] - su[int(xb), int(zbb)]) / (xa - xb))

   wt = rl * wi + math.sqrt(1 - rl **2) * math.sqrt(sw[int(xi), int(zii)]) * random.gauss(0, 1) + ( 1 - rl) * tl * ((sw[int(xa), int(zaa)] - sw[int(xb), int(zbb)]) / (za - zb))
       
   ui = u[int(xi), int(zi)] + ut 
   
   wi = w[int(xi), int(zi)] + wt
   
   xi = xi + ui * dt

   zi = zi + wi * dt
   
   xi, zi, ui, wi = reflexion(xi, zi, ui, wi, xa, za)

   print(xi, zi)
   
   return xi, zi, wi, ui, dt

def gitter(xi, xa, zi, za, dt):

   maxx = int(xi - xa)
   maxz = int(zi - za)
   
   for i in range(0, maxx):
      xsi = math.ceil(xa) + i
      tk.append((xsi - xa) / (xi - xa))

   for j in range(0, maxz):
      zsi = math.ceil(za) + j
      tk.append((zsi - za) / (zi - za))
   
   tks = sorted(tk)

   for i in range(1, len(tks)):
      
      tm = np.mean(np.array([tks[i-1], tks[i]]))

      xgg = xa + tm * (xi - xa)
      zgg = za + tm * (zi - za)
      
      if xgg > nx or zgg > nz:
         break
      else:
         
         gitter_monte[int(xgg), int(zgg)] += (tks[i] - tks[i-1]) * dt *((q * dt) / (n * dx * dz))
   return 

for i in tqdm(range(0, n)):
   x = []
   z = [] 
   tk = []
   xi = math.ceil(xq)
   zi = math.ceil(zq)
   ui = 0
   wi = 0
   dt = 0
   xb = 0
   zb = 0
   xa = 0
   za = 0

   while ((xi + ui * dt )) < nx:
      xb = xa
      zb = za

      xa = xi
      za = zi

      xi, zi, wi, ui, dt = positionen(xi, xa, zi, za, ui, wi, su, sw, u, w)

      gitter(xi, xa, zi, za, dt)

      x.append(xi)
      z.append(zi)

plt.figure()
fig, ax = plt.subplots()

xla=np.arange(0,120)
zla=np.arange(0,120)
ax.contour(xla, zla, gitter_monte)

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