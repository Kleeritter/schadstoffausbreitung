import math
from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
import random
import math


n= 1000#!Anzahl Partikel
ubalken = np.float64(5) #!m/s
wbalken =np.float64( 0) #!m/s
zq =np.float64( 45) #!m
xq =np.float64( 50) #!m
counter=np.float64(0)
xgrenz= np.float64(2000)# !m


ges=[]
tl = np.float64(100)  #s Zeit
dt = np.float64(4) # Zeitschritt
sigu= np.float64(0) #m/s
sigw= np.float64(0.39) #m/s
rl= math.exp(- dt/tl)
nx = 1000
ny = 50
nz = 20
dx = 10
dy = 10
dz = 10


c=np.zeros((nx, ny, nz))
def concentration():
    Q= 1.5e+5 #!540 kg/h also 5.4e+8mg/h und so 150000 
    ubalken= 5
    h= 100
    zq= 100
    xq=0
    yq=250
    fg=0.504
    fk=0.818
    gg=0.265
    gk=0.818
    for i in range( 1, nx+1):
        for j in range(1,ny+1):
            for k in range(1, nz+1):
       	        dey = fg*((dx*i)**fk)
       	        dez = gg*((dx*i)**gk)
                c[i-1,j-1,k-1] =(Q/(2* math.pi*dey*dez*ubalken))*math.exp((-((j*dy)-yq)**2)/(2*dey**2)) *(math.exp((-((k*dz)-h)**2)/(2*dez**2)) 
           +math.exp((-((k*dz)+h)**2)/(2*dez**2)))
    return

concentration()
print(c)
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
            rr=np.float64( random.gauss(0,1))
            ui=5
            xi= xi + ui*dt
            wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
            zi= zi + wi*dt
            posi.append([xi,zi])


        else:
            rr= np.float64(random.gauss(0,1))
            ui=5
            xi= xi + ui*dt
            wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
            zi= zi + wi*dt
            posi.append([xi,zi])
    ges.append(posi)
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
plt.savefig("Partikeltrajektorien.png", dpi=150)

