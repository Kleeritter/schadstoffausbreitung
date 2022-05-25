import math
from tqdm import tqdm

import matplotlib.pyplot as plt
import numpy as np
import random
print(random.gauss(0,1))


n= 1000#!Anzahl Partikel
ubalken = np.float64(5) #!m/s
wbalken =np.float64( 0) #!m/s
ubalkenzw = np.float32(5) #!m/s
wbalkenzw =np.float32( 0) #!m/s
zqzw =np.float32( 45) #!m
xqzw =np.float32( 50) #!m
xgrenzzw= np.float32(2000)# !m
zq =np.float64( 45) #!m
xq =np.float64( 50) #!m
xgrenz= np.float64(2000)# !m


ges=[]
gesz=[]
tl = np.float64(100)  #s Zeit
dt = np.float64(4) # Zeitschritt
sigu= np.float64(0) #m/s
sigw= np.float64(0.39) #m/s
rl= math.exp(- dt/tl)


for i in tqdm(range(n)):
    xi=xq
    zi=zq
    xizw=xqzw
    zizw=zqzw
    ui=ubalken
    wi=wbalken
    uizw=ubalkenzw
    wizw=wbalkenzw
    posi=[]
    posizw=[]
    diff=[]
    diffx=[]
    diffz=[]
    while (xi<= xgrenz):
        if (zi<0):
            zi=-zi
            wi= -wi
            rr=np.float64( random.gauss(0,1))
            rrzw=np.float32(rr)
            ui=5
            xi= xi + ui*dt
            wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
            zi= zi + wi*dt
            posi.append([xi,zi])
            xizw= xizw + uizw*dt
            wizw= rl*wizw + math.sqrt((1 - rl**2))*sigw* rrzw
            zizw= zizw + wizw*dt
            posizw.append([xizw,zizw])
            difx= xi/xizw
            difz= zi/zizw
            diffx.append([xi,difx])
            diffz.append([zi,difz])
        else:
            rr= np.float64(random.gauss(0,1))
            rrzw=np.float32(rr)
            ui=5
            xi= xi + ui*dt
            wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
            zi= zi + wi*dt
            posi.append([xi,zi])
            xizw= xizw + uizw*dt
            wizw= rl*wizw + math.sqrt((1 - rl**2))*sigw* rrzw
            zizw= zizw + wizw*dt
            posizw.append([xizw,zizw])
            difx= xi/xizw
            difz= zi/zizw
            diffx.append([xi,difx])
            diffz.append([zi,difz])
    ges.append(diffx)
    gesz.append(diffz)
for i in tqdm(range(len(ges))):
    x=[]
    z=[]
    for j in range(len(ges[i])):
        x.append(ges[i][j][0])
        z.append(ges[i][j][1])
    plt.plot(x,z)

plt.title("Quotient DP/ SP für X Koordinate")
plt.xlabel("X in m")
plt.ylabel("Quotient")
plt.show()
plt.savefig("qX.png", dpi=150)

for i in tqdm(range(len(gesz))):
    x=[]
    z=[]
    for j in range(len(gesz[i])):
        x.append(gesz[i][j][0])
        z.append(gesz[i][j][1])
    plt.plot(x,z)

plt.title("Quotient DP/ SP für Z Koordinate")
plt.xlabel("Z in m")
plt.ylabel("Quotient")
plt.show()
plt.savefig("qZ.png", dpi=150)
