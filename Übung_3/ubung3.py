import math
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import math


n= 1000 #!Anzahl Partikel
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
ui=5

bins=np.arange(0, 2002,2)
positionsliste=[]
xvalues=[]
zvalues=[]
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

def positionen(xi,wi,zi):
    
    rr=np.float64( random.gauss(0,1))
    xi= xi + ui*dt
    wi= rl*wi + math.sqrt((1 - rl**2))*sigw* rr
    zi= zi + wi*dt
    return xi,wi,zi

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
            xi,zi,wi =positionen(xi,wi,zi)
            posi.append([xi,zi])
            xvalues.append(xi)
            zvalues.append(zi)
            positionsliste.append([xi,zi])


        else:
            #print(xi)
            xi,wi,zi =positionen(xi,wi,zi)
            #print(xi)
            posi.append([xi,zi])
            xvalues.append(xi)
            zvalues.append(zi)
            positionsliste.append([xi,zi])

    ges.append(posi)
    """
    df=pd.DataFrame.from_dict({
        'xvalues': xvalues,
        'zvalues': zvalues
    })
    df['xvaluesort']=pd.cut(df['xvalues'],bins, precision=6)
    df['zvaluesort']=pd.cut(df['zvalues'],bins, precision=6)
    """
    """
    df=pd.DataFrame.from_dict({
        'values': posi,
    })
df['valuesort']=pd.qcut(df['values'],100, precision=6)
"""
xvalues.pop()
zvalues.pop()
df=pd.DataFrame.from_dict({
        'xvalues': xvalues,
        'zvalues': zvalues
    })
df['xvaluesort']=pd.cut(df['xvalues'],bins, precision=6)
df['zvaluesort']=pd.cut(df['zvalues'],bins, precision=6)
#df['values'] = df['xvaluesort'] + df['zvaluesort']
sx= pd.Series(xvalues)
sx.value_counts(bins=bins)#.sort_index()
print(df)

#alla=pd.Series.to_frame(df.value_counts(["xvaluesort", "zvaluesort"]))#.rename(columns={4:"times"})
alla=df.value_counts(["xvaluesort", "zvaluesort"])
#alla.columns=["xvaluesort"]
#alla.rename(columns={alla.columns[0]:'times'})#,inplace=True)
#print(alla["(68,70] (44,46]"])
print(alla)
klar=[]
for i in range(len(alla)):
    klar.append([alla.index.values[i][0].mid,alla.index.values[i][1].mid,alla.values[i]])
print(klar)
print(alla.index.values[0][1].mid)
print(alla.index.values[0][0].mid)
for i in tqdm(range(len(ges))):#
    x=[]
    z=[]
    for j in range(len(ges[i])):
        x.append(ges[i][j][0])
        z.append(ges[i][j][1])
    plt.plot(x,z)

plt.title("Partikeltrajektorien")
plt.xlabel("Distanz  X in m")
plt.ylabel("Höhe Z in m")
#plt.show()
plt.savefig("Partikeltrajektorien.png", dpi=150)

