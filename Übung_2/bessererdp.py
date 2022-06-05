import math
from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
import random


n= 1000#!Anzahl Partikel
ubalken = np.float64(5) #!m/s
wbalken =np.float64( 0) #!m/s
zq =np.float64( 45) #!m
xq =np.float64( 50) #!m
counter=np.float64(0)
xgrenz= np.float64(2000)# !m


ges=[]
tl = np.float64(100)  #s Zeit
dt = np.float64(0.4) # Zeitschritt
sigu= np.float64(0) #m/s
sigw= np.float64(0.39) #m/s
rl= math.exp(- dt/tl)


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

