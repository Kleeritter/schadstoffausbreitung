import math
from tqdm import tqdm

import matplotlib as plt
import numpy as np
import random
print(random.gauss(0,1))


n= 1000 #!Anzahl Partikel
ubalken = 5 #!m/s
wbalken = 0 #!m/s
zq = 45 #!m
xq = 2000 #!m
counter=0
xgrenz= 3000# !m
posi=np.array([])
ges=np.array([])
def berechnen(xi,ui,wi,zi,posi):
    tl = 100  #s Zeit
    dt = 4 # Zeitschritt
    sigu= 0 #m/s
    sigw= 0.39 #m/s

    rl= math.exp(- dt/tl)
    rr= random.gauss(0,1)
    ui= rl*ui + math.sqrt((1 - rl**2))*sigu* rr
    xi= xi + ui*dt
    wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
    zi= zi + wi*dt
    transport=np.array([xi,zi])
    return np.append(posi, transport), xi, zi, ui,wi
for i in tqdm(range(n)):
    xi=xq
    zi=zq
    ui=ubalken
    wi=wbalken

    while (xi<= xgrenz):
        if (zi<0):
            zi=-zi
            wi= -wi
            berechnen(xi,zi,posi,ui,wi)
            print("goose")
            ges=np.append(ges,posi)

        else:
            berechnen(xi,zi,posi,ui,wi)
            print(xi)
            ges=np.append(ges,posi)
    #counter=counter +1

for i in range(len(ges)):
    x=np.arrange(0,len(ges[i]))
    plt.plot(x,ges[i])

plt.show()

