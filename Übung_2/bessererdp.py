import math
from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
import random
print(random.gauss(0,1))


n= 1000 #!Anzahl Partikel
ubalken = 5 #!m/s
wbalken = 0 #!m/s
zq = 45 #!m
xq = 2000 #!m
counter=0
xgrenz= 2492# !m
#posi=np.array([])
#ges=np.array([])

ges=[]
tl = 100  #s Zeit
dt = 4 # Zeitschritt
sigu= 0 #m/s
sigw= 0.39 #m/s
rl= math.exp(- dt/tl)
print(rl)
xliste=[]
zliste =[]

for i in tqdm(range(n)):
    xi=xq
    zi=zq
    ui=ubalken
    wi=wbalken
    posi=[]

    while (xi<= xgrenz):
        if (zi<0):
            #print("alarm")
            zi=-zi
            wi= -wi
            rr= random.gauss(0,1)
            ui= rl*ui + math.sqrt((1 - rl**2))*sigu* rr
            xi= xi + ui*dt
            wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
            zi= zi + wi*dt
            #transport=np.array([xi,zi])
            #posi=np.append(posi, (xi,zi))
            #ges=np.append(ges,posi)
            posi.append([xi,ui])
            xliste.append(ui)
            zliste.append(wi)

        else:
            rr= random.gauss(0,1)
            ui= rl*ui + math.sqrt((1 - rl**2))*sigu* rr
            xi= xi + ui*dt
            wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
            zi= zi + wi*dt
            #transport=np.array([xi,zi])
            #posi=np.append(posi, (xi,zi))
            #ges=np.append(ges,posi)
            posi.append([xi,ui])
            xliste.append(ui)
            zliste.append(wi)
    #counter=counter +1
    ges.append(posi)
print(ges[0])
for i in tqdm(range(len(ges))):
    x=[]
    z=[]
    for j in range(len(ges[i])):
        x.append(ges[i][j][0])
        z.append(ges[i][j][1])
    plt.plot(x,z)
    #x=np.arange(0,len(ges[i]))
    #plt.plot(x,ges[i])

plt.show()

