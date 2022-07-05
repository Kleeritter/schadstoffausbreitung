from collections import UserString
import math
from matplotlib.pyplot import tick_params
from tqdm import tqdm
#import pandas as pd
#mport matplotlib.pyplot as plt
import numpy as np
import random
import math
#import plotly.graph_objects as go
from netCDF4 import Dataset
import netCDF4 as nc4
#from plotly.subplots import make_subplots
#import numba # CUDA zur Auslagerung auf die GPU
#from numba import cuda
import vergleich_gauss_mc4
from vergleich_gauss_mc4 import vergleich #Visualisierung der NC Datei im Vergleich

def ubu3 ():
    n= 2000#!Anzahl Partikel
    ubalken = np.float64(5) #!m/s
    wbalken =np.float64( 0) #!m/s
    xq =np.float64( 0.5) #!m
    zq =np.float64( 0.008) #!m
    counter=np.float64(0)
    xgrenz= 110# !m
    zgrenz=25
    ges=[]
    tl = np.float64(100)  #s Zeit
    dt = np.float64(0.4) # Zeitschritt

    
    nx = 1000
    ny = 50
    nz = 2000
    dx = 2
    dy = 2
    dz = 2
    ui=5
    ustern=0.35
    k=0.38
    znull=0.008
    sigu= 2.5*ustern #m/s
    sigw= 1.3*ustern #m/s
    #sigu= np.float64(0) #m/s
    #sigw= np.float64(0.39) #m/s



    gitter = np.zeros((int(xgrenz/2 +1),int(zgrenz/2 + 1)))
    positionsliste=[]
    xvalues=[]
    zvalues=[]

    nxx=int(xgrenz/dx +1)
    nzz=int(zgrenz/dz +1)

    cd= np.zeros((nxx,nzz))



    q= 150
    dt= 0.4
    def zuordnen(vi):
        vr = round(vi)      # Runden der Partikelposition, um anschließend zu überprüfen, zu welchem Gittermittelpunkt zugeordnet werden muss
        if vr > vi:
            if (vr % 2) == 0:   # Fall: Aufgerundet und gerade
                vm = vr - 1
            else:
                vm = vr         # Fall: Aufgerundet und ungerade
            
        elif vr < vi:
            if (vr % 2) == 0:
                vm = vr + 1     # Fall: Abgerundet und gerade
            else:
                vm = vr         # Fall: Abgerundet und ungerade
        else:
            if (vr % 2) == 0:
                vm = vr + 1     # Fall: Bereit gerade
            else:
                vm = vr         # Fall: Bereits ungerade

        return vm
    
    def gita(vi):
        vr = round(vi)      # Runden der Partikelposition, um anschließend zu überprüfen, zu welchem Gittermittelpunkt zugeordnet werden muss
        if vr > vi:
            if (vr % 2) == 0:   # Fall: Aufgerundet und gerade
                vm = vr - 1 
            else:
                vm = vr         # Fall: Aufgerundet und ungerade
            
        elif vr < vi:
            if (vr % 2) == 0:
                vm = vr + 1     # Fall: Abgerundet und gerade
            else:
                vm = vr         # Fall: Abgerundet und ungerade
        else:
            if (vr % 2) == 0:
                vm = vr + 1     # Fall: Bereit gerade
            else:
                vm = vr         # Fall: Bereits ungerade

        return vm

    def prandl(zi):
        if zi ==0:
            ubalken= 0
        else:
            ubalken= (ustern/k)*math.log(abs(zi)/znull)
            
        return ubalken

    def prandltl(zi):
        tl= ((k*ustern)/sigw**2) *abs(zi)
        return tl

    def prairie():

        return
    def exaktgitter(xi,xold,zi,zold):
        maxstep= math.ceil(xi)-math.floor()

        xsi=math.ceil(xi)
        zsi=math.ceil(zi)

        ti=(xsi -xold)/(xi-xold)
        tj=(zsi -zold)/(zi-zold)
        tku=sorted(ti,tj)
        return
    def positionen(xi,wi,zi,tl,ui):
        rl= math.exp(- dt/tl)
        rr=np.float64( random.gauss(0,1))
        xi= xi + ui*dt
        wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
        zi= zi + wi*dt
        #print(xi,zi)
        return xi,wi,zi

    def konki(xi,zi,xold,zold):
        #for k in range():
        t=1/xi * dt
        #poser=(xold,zold)+ t*(xi -xold,zi-zold)
        #langposer= np.dot(poser)
        #print(langposer)
    def gitweis(xi,zi):
        xm = int((zuordnen(xi) - 1) / dx)
        zm = int((zuordnen(zi) - 1) / dz)
        gitter[xm,zm] = gitter[xm,zm] + 1
        return 




    def main():
        for i in tqdm(range(n)):
            
            xi=xq
            zi=zq
            posi=[]
            ui=ubalken
            wi=wbalken

            while (math.ceil(xi)< xgrenz):
                #print(zi)
                xold=xi
                zold=zi
                
                #ui=uberechnung(zi)
                if (zi<znull):
                    zi=-zi
                    wi= -wi
                    tl=prandltl(zi)
                    #tl=
                    ui=prandl(zi)
                    xi,wi,zi =positionen(xi,wi,zi,tl,ui)
                    #xm = int((zuordnen(xi) - 1) / dx)
                    #zm = int((zuordnen(zi) - 1) / dz)
                    #gitter[xm,zm] = gitter[xm,zm] + 1
                    gitweis(xi,zi)
                    
                    
                else:
                    tl=prandltl(zi)
                    ui=prandl(zi)
                    xi,wi,zi =positionen(xi,wi,zi,tl,ui)
                    xm = int((zuordnen(xi) - 1) / dx)
                    zm = int((zuordnen(zi) - 1) / dz)
                    gitter[xm,zm] = gitter[xm,zm] + 1
                   # konki(xi,zi,xold,zold)
                    gitweis(xi,zi)
                    
           
            #print(i)
    main()
    return [i * ((q * dt)/(n * dx * dz)) for i in gitter]

konzentrationen=ubu3()
xgrenz= 110# !m
zgrenz=25
dx=2
dz=2
nxx=int(xgrenz/dx +1)
nzz=int(zgrenz/dz +1)
d = nc4.Dataset('prada.nc','w', format='NETCDF4') #'w' stands for write
d.createDimension('x', nzz)
d.createDimension('z', nxx)
cnet= d.createVariable('c', 'f8', ('x','z'))
xnet= d.createVariable('x', 'f8', 'x')
znet= d.createVariable('z', 'f8', 'z')

cnet[:]=np.transpose(konzentrationen)
znet[:]=np.arange(0.,nxx*dx,dx,dtype=float)
xnet[:]=np.arange(0.,nzz*dz,dz,dtype=float)
d.close()

vergleich_gauss_mc4.vergleich()