from collections import UserString
import math
from matplotlib.pyplot import tick_params
from tqdm import tqdm

import numpy as np
import random
import math

from netCDF4 import Dataset
import netCDF4 as nc4

import matplotlib.pyplot as plt

ges=[]

n= 1000#!Anzahl Partikel
ubalken = np.float64(5) #!m/s
wbalken =np.float64( 0) #!m/s
xq =np.float64( 0) #!m
zq =np.float64( 0.5) #!m
counter=np.float64(0)
xgrenz= 110# !m
zgrenz=25
ges=[]
tl = np.float64(100)  #s Zeit
dt = np.float64(0.4) # Zeitschritt


nx = 1000
ny = 50
nz = 2000
dx = 1
dy = 1
dz = 1
ui=5
ustern=0.35
k=0.38
znull=0.008
sigu= 2.5*ustern #m/s
sigw= 1.3*ustern #m/s

gitter = np.zeros((int(xgrenz*2),int(zgrenz*2)))
positionsliste=[]
xvalues=[]
zvalues=[]
ges=[]

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
          if vr >= xgrenz or zgrenz:
             vm = vr    # Fall: Abgerundet und gerade
          else:
             vm = vr + 1
        else:
            vm = vr         # Fall: Abgerundet und ungerade
    else:
        if (vr % 2) == 0:
          if vr >= xgrenz or zgrenz:
             vm = vr    # Fall: Abgerundet und gerade
          else:
             vm = vr + 1   # Fall: Bereit gerade
        else:
            vm = vr         # Fall: Bereits ungerade

    return vm

def gg(xold,zold,xi,zi,t):
    #print(t)
    xg=xold +t*(xi-xold)
    zg=zold +t*(zi-zold)
    return int(xg),int(zg)

def prandl(zi):
    if zi ==0:
        ubalken= 0
    else:
        ubalken= (ustern/k)*math.log(abs(zi)/znull)
        
    return ubalken

def prandltl(zi):
    tl= ((k*ustern)/sigw**2) *abs(zi)
    #dt=0.1*tl
    return tl#,dt

def prairie():

    return
def rangecheck(xi,xold,zi,zold):
    rangex=int(xi-xold)
    rangez=int(zi-zold)
    #print(rangez,rangex)
    if (rangex + rangez) <2:
        #print("oof")
        gitweis(xi,zi)
    else: 
        exaktgitter(xi,xold,zi,zold)
def exaktgitter(xi,xold,zi,zold):
    #print(xi,xold)
    #print(zi,zold)
    ti=[]
    tj=[]
    toks=[]  
    rangex=int(xi-xold)
    rangez=int(zi-zold)
    
    for i in range(0,rangex):
        if i==0:
            xsi=math.ceil(xold)
            toks.append((xsi -xold)/(xi-xold))
        else:
            xsi+=1
            toks.append((xsi -xold)/(xi-xold))
    #toks.append(ti)

    for i in range(0,rangez):
        #print("geht iwas")
        if i==0:
            zsi=math.ceil(zold)
            #print("alla")
            #print(zsi)
            toks.append((zsi -zold)/(zi-zold))
        else:
            zsi+=1
            toks.append((zsi -zold)/(zi-zold))
    #toks.append(tj)
    tku=sorted(toks)
    #print(tku)
    for i in range(1,len(tku)):
        ti=tku[i]
        told=tku[i-1]
        t=np.mean(np.array([told, ti]))
        posx,posz=gg(xold,xi,zold,zi,t)
        gitter[posz,posx]+=(tku[i]-tku[i-1]*dt)  
    return
def positionen(xi,wi,zi,tl,ui,):
    rl= math.exp(- dt/tl)
    rr=np.float64( random.gauss(0,1))
    xi= xi + ui*dt
    wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
    zi= zi + wi*dt
    #print(xi,zi)
    return xi,wi,zi


def gitweis(xi,zi):
    xm = int((zuordnen(xi) - 1) / dx)
    zm = int((zuordnen(zi) - 1) / dz)
    if xm >int(xgrenz/2):
      xm=int(xgrenz/2)
    if zm >int(zgrenz/2):
      zm=int(zgrenz/2)
    gitter[xm,zm] = gitter[xm,zm] + 1
    return 

for i in tqdm(range(n)):
    
    xi=xq
    zi=zq
    posi=[]
    wi=wbalken

    while (math.ceil(xi)< xgrenz):
        xold=xi
        zold=zi
        
        if (zi<0):
            zi=-zi
            wi= -wi
            tl=prandltl(zi)
            ui=prandl(zi)
            xi,wi,zi =positionen(xi,wi,zi,tl,ui)#,dt)
            rangecheck(xold,zold,xi,zi)
            posi.append([xi,zi])           
            
        else:
            tl=prandltl(zi)
            ui=prandl(zi)
            xi,wi,zi =positionen(xi,wi,zi,tl,ui)#,dt)
            #gitweis(xi,zi)
            rangecheck(xold,zold,xi,zi)
            posi.append([xi,zi])
            
    ges.append(posi)


#konzentrationen=[i * ((q * dt)/(n * dx * dz)) for i in gitter]
d = nc4.Dataset('prada.nc','w', format='NETCDF4') #'w' stands for write
d.createDimension('x', 2*zgrenz)
d.createDimension('z', 2*xgrenz)
cnet= d.createVariable('c', 'f8', ('x','z'))
xnet= d.createVariable('x', 'f8', 'x')
znet= d.createVariable('z', 'f8', 'z')

cnet[:]=np.transpose(gitter)
znet[:]=np.arange(0.,2*xgrenz,dx,dtype=float)
xnet[:]=np.arange(0.,2*zgrenz,dz,dtype=float)
d.close()

for i in tqdm(range(len(ges))):
    x=[]
    z=[]
    for j in range(len(ges[i])):
        x.append(ges[i][j][0])
        z.append(ges[i][j][1])
    plt.plot(x,z)

plt.title("Partikeltrajektorien")
plt.xlabel("Distanz  X in m")
plt.ylabel("Höhe Z in m")
plt.show()
#plt.savefig("Test")

