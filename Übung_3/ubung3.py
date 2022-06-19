import math
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import math
import plotly.graph_objects as go
from netCDF4 import Dataset
import netCDF4 as nc4
from plotly.subplots import make_subplots

n= 20000#!Anzahl Partikel
ubalken = np.float64(5) #!m/s
wbalken =np.float64( 0) #!m/s
zq =np.float64( 45) #!m
xq =np.float64( 51) #!m
counter=np.float64(0)
xgrenz= np.float64(2000)# !m
zgrenz=400
ges=[]
tl = np.float64(100)  #s Zeit
dt = np.float64(0.4) # Zeitschritt
sigu= np.float64(0) #m/s
sigw= np.float64(0.39) #m/s
rl= math.exp(- dt/tl)
nx = 1000
ny = 50
nz = 2000
dx = 2
dy = 2
dz = 2
ui=5


gitter = np.zeros((int(xgrenz/2 + 1),int(zgrenz/2 + 1)))
bins=np.arange(0, 2000,2)
positionsliste=[]
xvalues=[]
zvalues=[]

cdc=[]
nxx=int(xgrenz/dx +1)
nzz=int(zgrenz/dz +1)
cdkack=[]
cd=np.zeros((nxx, nzz))
def gauss2d():
    Q= 150 #!540 kg/h also 5.4e+8mg/h und so 150000 
    ubalken= 5
    dx=2
    dz=2
    for i,indexei in zip(range(1 , nzz),(range(1,nxx,2))):
        for j,indexej in zip(range(1,nxx),(range(1,nxx,2))):
            #dez = gg*((dx*i)**gk)
       	    #dez = 2*sigw**2*tl*((i/ubalken)- ubalken +ubalken*math.exp(- (i/(ubalken*tl))))
            dez= 2*sigw**2*tl*((i/ubalken)-tl +tl*math.exp(-i/(ubalken*tl)))
            cd[j,i]=(Q/(math.sqrt(2* math.pi)*math.sqrt(dez)*ubalken) *(math.exp((-((j*dx)-zq)**2)/(2*dez)) +math.exp((-((j*dz)+zq)**2)/(2*dez))))
            #cdc.append(co)

    return


q= 150
dt= 0.4
def concentration(nj):
    dx=2
    dz=2
    q= 150
    dt= 0.4
    c= q*(nj*dt)/(n*dx*dz)

    return c
def positionen(xi,wi,zi):
    
    rr=np.float64( random.gauss(0,1))
    xi= xi + ui*dt
    wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
    zi= zi + wi*dt
    return xi,wi,zi


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
for i in tqdm(range(n)):
    xi=xq
    zi=zq
    ui=ubalken
    wi=wbalken
    posi=[]


    while (xi<= xgrenz):
        if (zi<0):
            zi=-zi
            wi= -wi
            xi,wi,zi =positionen(xi,wi,zi)
            posi.append([xi,zi])
            xvalues.append(xi)
            zvalues.append(zi)
            positionsliste.append([xi,zi])
            xm = int((zuordnen(xi) - 1) / dx)
            zm = int((zuordnen(zi) - 1) / dz)
            gitter[xm,zm] = gitter[xm,zm] + 1


        else:
            xi,wi,zi =positionen(xi,wi,zi)
            posi.append([xi,zi])
            xvalues.append(xi)
            zvalues.append(zi)
            positionsliste.append([xi,zi])
            xm = int((zuordnen(xi) - 1) / dx)
            zm = int((zuordnen(zi) - 1) / dz)
            gitter[xm,zm] = gitter[xm,zm] + 1

    ges.append(posi)

konzentrationen = [i * ((q * dt)/(n * dx * dz)) for i in gitter]
 

gauss2d()

f = nc4.Dataset('koks.nc','w', format='NETCDF4') #'w' stands for write
#tempgrp = f.createGroup('Temp_data')
#tempgrp.createDimension('c', (n))
f.createDimension('x', nxx)
f.createDimension('z', nzz)
cnet= f.createVariable('c', 'f8', ('x','z'))
xnet= f.createVariable('x', 'f8', 'x')
znet= f.createVariable('z', 'f8', 'z')

cnet[:,:]=cd
xnet[:]=[x for x in range(0,nxx)]
znet[:]=[x for x in range(0,nzz)]
f.close()


d = nc4.Dataset('mc.nc','w', format='NETCDF4') #'w' stands for write
d.createDimension('x', nzz)
d.createDimension('z', nxx)
cnet= d.createVariable('c', 'f8', ('x','z'))
xnet= d.createVariable('x', 'f8', 'x')
znet= d.createVariable('z', 'f8', 'z')

cnet[:,:]=np.transpose(konzentrationen)
znet[:]=np.arange(0.,nxx*dx,dx,dtype=float)
xnet[:]=np.arange(0.,nzz*dz,dz,dtype=float)
d.close()

"""
lig = make_subplots(    rows=2,
    specs=[[{"type": "contour"}],
            [{"type": "contour"}]],
)

lig.add_trace(
    go.Contour(
        z=cplot,
        x=xplot,#.insert(, # horizontal axis
        y=zplot, # vertical axis
        contours_coloring='lines',
        contours=dict(
            start=0.1,
            end=0.5,
            size=0.1,
        ),
    ), row=1, col=1)


#lig.show()
#lig.write_image("Übung_3/conturallagaus.png")
"""