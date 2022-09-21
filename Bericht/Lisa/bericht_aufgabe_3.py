import math as m
import random as r
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from tqdm import tqdm
from numba import jit
##### Einlesen #####
@jit
def main():
   def nc_2d(filename, varname):
      f = open(filename)
      f.close()
      nc_file = Dataset(filename, "r", format = "NETCDF4")
      tmp_array = np.array(nc_file.variables[varname][:, :], 
      dtype = type(nc_file.variables[varname]))
      return tmp_array

   filename = "/Users/alex/Documents/Dokumente - Mac mini von Alex/GitHub/schadstoffausbreitung/Bericht/Lisa/input_uebung5.nc"

   #Horizontalgeschwindigkeit
   u = np.transpose(nc_2d(filename, "u"))
   u = np.where(u == -9999.0, np.nan, u)
   #Vertikalgeschwindigkeit
   w = np.transpose(nc_2d(filename, "w"))
   w = np.where(w == -9999.0, np.nan, w)
   #Standardabweichung der Horizontalgeschwindigkeit
   su = np.transpose(nc_2d(filename, "u2"))
   #Standardabweichung der Vertikalgeschwindigkeit
   sw = np.transpose(nc_2d(filename, "w2"))
   # Anzahl Partikel
   n = 1000
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
   # lambda
   l = 5
   #Gitter-Array
   gitter_mc = np.zeros((nx, nz))

   def reflexion(xi, zi, ui, wi):
      while(xi<=32 and zi<=61)or(xi>=89 and zi<=61)or(32<=xi<=89 and zi<= 1):
         #linke Wand
         if (xi <= 32 and zi <= 61):
            xi = xi + 2 * (abs(32 - xi)) 
            ui = -ui
         #rechte Wand
         elif (xi >= 89 and zi <= 61):
            xi = xi - 2 * (abs(89 - xi)) 
            ui = -ui
         #Boden
         elif (32 <= xi <= 89 and zi <= 1):
            zi = zi + 2 * (abs(1 - zi)) 
            wi = -wi
         else:
            xi = xi
            zi = zi
            wi = wi
            ui = ui
         
      return xi, zi, ui, wi

   def pos(xi, xa, xb, zi, za, zb, ui, wi, dt, su, sw, u, w):

      if su[int(xi),int(zi)] == 0 or sw[int(xi),int(zi)] == 0:
         tl = 0.0001 #sehr klein wählen, da sonst dt = 0
      else: 
         tl_one = 0.05*(kappa*zi)/((1+kappa*(zi/l))) 
         tl_two = (0.23*m.sqrt(su[int(xi),int(zi)]+sw[int(xi),int(zi)]))
         tl = tl_one / tl_two #T_L(z)

         tl_2m_one = 0.05*(kappa*2)/((1+kappa*(2/l)))
         tl_2m_two = (0.23*m.sqrt(su[int(xi),int(zi)]+sw[int(xi),int(zi)]))
         tl_2m = tl_2m_one / tl_2m_two  #T_L(2m)

         if (0.1 * tl) > (tl_2m):
               dt = 0.1 * tl
         else:
               dt = tl_2m
      
   #lagrangsche Autokorrelationsfunktion
      rl = m.exp(- dt / tl)
      
      ut_alt = rl*ui+m.sqrt(1-rl**2)*m.sqrt(su[int(xi),int(zi)])*r.gauss(0, 1)
      ut_neu = (1-rl)*tl*((su[int(xa),int(za)]-su[int(xb),int(zb)])/(xa-xb))
      ut =  ut_alt + ut_neu #Abweichung vom Grundzustand
      ui = u[int(xi), int(zi)] + ut #Grundzustand und Abweichung

      wt_alt = rl*wi+m.sqrt(1-rl**2)*m.sqrt(sw[int(xi),int(zi)])*r.gauss(0, 1)
      wt_neu = (1-rl)*tl*((sw[int(xa),int(za)]-sw[int(xb),int(zb)])/(za-zb))
      wt = wt_alt + wt_neu #Abweichung vom Grundzustand
      wi = w[int(xi), int(zi)] + wt #Grundzustand und Abweichung

      xi = xi + ui * dt
      zi = zi + wi * dt

      xi, zi, ui, wi = reflexion(xi, zi, ui, wi)

      return xi, zi, wi, ui, dt

   def gitter(xi, xa, zi, za, dt):

      #Anzahl an Schnittpunkten
      maxx = abs(int(xi) - int(xa))
      maxz = abs(int(zi) - int(za))
   
      #Berechnung von ti, tj
      for i in range(0, maxx):
         if xi < xa:
            xsi = xa - i
         else:
            xsi = xa + i
         tk.append((xsi - xa) / abs(xi - xa))
      
      for j in range(0, maxz):
         if zi < za:
            zsi = za - j
         else:
            zsi = za + j
         tk.append((zsi - za) / abs(zi - za))

      #Sortieren
      tks = sorted(tk)
      
      if maxx >= 1 and maxz >= 1:
         #Einsetzen in die Geradengleichung

         for i in range(1, len(tks)):
            tm = np.mean(np.array([tks[i-1], tks[i]]))
            xgg = xa + tm * abs(xi - xa)
            zgg = za + tm * abs(zi - za)
            if xgg > nx:
               xgg = nx - 1

            if zgg > nz:
               zgg = nz - 1
            #Berechnen der Konzentration
            tkk = (tks[i] - tks[i-1]) * dt
            c = ((q * dt) / (n * dx * dz))
            gitter_mc[int(xgg), int(zgg)] += tkk * c
      else:
         if xi > nx:
            xi = nx - 1
         if zi > nz:
            zi = nz - 1
         c = ((q * dt) / (n * dx * dz))

         gitter_mc[int(xi), int(zi)] += c
            
      return gitter_mc

   for i in tqdm(range(0, n)):
      x = []
      z = [] 
      tk = []
      xi = m.ceil(xq)
      zi = m.ceil(zq)
      ui = 0
      wi = 0
      dt = 0
      xb = 0
      zb = 0
      xa = 0
      za = 0

      while ((xi + ui * dt ) < nx) and ((zi + wi * dt ) < nz) :
         xb = xa
         zb = za
         
         xa = xi
         za = zi
         
         xi,zi,wi,ui,dt = pos(xi, xa, xb, zi, za, zb, ui, wi, dt, su, sw, u, w)

         gitter(xi, xa, zi, za, dt)
         x.append(xi)
         z.append(zi)
      #plt.plot(x, z)
   #plt.show()

   conc = np.transpose(gitter_mc)

   x = np.arange(0, 120)
   z = np.arange(0, 120)

   fig, ax = plt.subplots()
   levels = [0.01, 0.025, 0.05, 0.5, 0.75, 1.0, 1.25, 1.5,  2]
   cs = ax.contour(x, z, conc, levels, colors='black')
   ax.clabel(cs, fontsize=9, inline=1)
   units = r"$\mathrm{\frac{kg}{s}}$"
   plt.title('Konzentration einer Straßenschlucht ('+units+')')
   plt.grid()
   plt.xlabel('x(m)')
   plt.ylabel('z(m)')
   plt.savefig('bericht_aufgabe_3a.png', format = 'png', dpi = 300)

main()