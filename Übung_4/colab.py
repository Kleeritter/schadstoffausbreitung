from tqdm import tqdm
import numpy as np
import random
import math
import netCDF4 as nc4
import matplotlib.pyplot as plt
from grassmodul import grass
ges = [] # was soll das sein?
n = 10**5 # !Anzahl Partikel
ubalken = np.float64(5)  # !m/s
wbalken = np.float64(0)  # !m/s
xq = np.float64(0)  # !m
zq = np.float64(0.5)  # !m
counter = np.float64(0)
xgrenz = 110  # !m
zgrenz = 25

# mach die mal weg zur übersicht
nx = 1000
ny = 50
nz = 2000
dx = 1
dy = 1
dz = 1

ustern = 0.35
k = 0.38
znull = 0.008
sigu = 2.5 * ustern  # m/s
sigw = 1.3 * ustern  # m/s

gitter = np.zeros((int(xgrenz), int(zgrenz)))
konk = np.zeros((int(xgrenz), int(zgrenz)))
q = 150

def gg(xold, zold, xi, zi, t):

    xg = xold + t * (xi - xold)
    zg = zold + t * (zi - zold)
    #print("der wert",xg, t, xold,xi)
    if xg>xgrenz:
        xg=xgrenz-1
    return int(xg), int(zg)

def prandl(zi):
    if zi < znull:
        ubalken = 0
    else:
        ubalken = (ustern / k) * math.log(abs(zi) / znull)

    return ubalken

def prandltl(zi):
    tl = ((k * ustern) / sigw ** 2) * abs(zi)
    if (0.1*tl)>((k * ustern) / sigw ** 2) * abs(2): #falls dt kleiner als tl in 2 m Höhe
        dt = 0.1*tl
    else:
        dt = 0.1*((k * ustern) / sigw ** 2) * abs(2)

    return tl,  dt

def rangecheck(xi, xold, zi, zold):
    rangex = int(xi - xold)
    rangez = int(zi - zold)
    # print(rangez,rangex)
    if (rangex + rangez) < 2:
        #print("oof")
        gitweis(xi, zi)
    else:
        exaktgitter(xi, xold, zi, zold)

def exaktgitter(xi, xold, zi, zold):
    ti = []
    tj = []
    toks = []
    rangex = int(xi - xold)
    rangez = int(zi - zold)
    #print(xi,xold,zi,zold,rangex,rangez)

    for i in range(0, rangex):
        if i == 0:
            xsi = math.ceil(xold)
            toks.append((xsi - xold) / (xi - xold))
        else:
            xsi += 1
            toks.append((xsi - xold) / (xi - xold))
    # toks.append(ti)

    for i in range(0, rangez):
        # print("geht iwas")
        if i == 0:
            zsi = math.ceil(zold)
            # print("alla")
            # print(zsi)
            toks.append((zsi - zold) / (zi - zold))
        else:
            zsi += 1
            toks.append((zsi - zold) / (zi - zold))
    tku = sorted(toks)
    #print(tku)
    for i in range(1, len(tku)):
        ti = tku[i]
        told = tku[i - 1]
        #print((tku[i] - tku[i - 1] * dt))
        t = np.mean(np.array([told, ti]))
        #print(ti, told,t)
        posx, posz = gg(xold, zold, xi, zi, t)
        if posz>zgrenz or posx>xgrenz:
            break
        gitter[posx, posz] += (tku[i] - tku[i-1] )* dt
        konk[posx, posz] += (tku[i] - tku[i-1])* dt*((q * dt)/(n * dx * dz))
    return


def positionen(xi, wi, zi, tl, ui, dt ):
    rl = math.exp(- dt / tl)
    rr = np.float64(random.gauss(0, 1))
    xi = xi + ui * dt
    wi = rl * wi + math.sqrt((1 - rl ** 2)) * sigw * rr
    zi = zi + wi * dt
    return xi, wi, zi


def gitweis(xi, zi):
    xm = int((xi)-1)
    zm = int((zi)-1)
    #if xm > int(xgrenz / 2):
     #   xm = int(xgrenz / 2)
    #if zm > int(zgrenz / 2):
     #   zm = int(zgrenz / 2)
    gitter[xm, zm] = gitter[xm, zm] + 1
    konk[xm, zm] += ((q * dt)/(n * dx * dz))
    return

for i in tqdm(range(n)):
    xi = xq
    zi = zq
    posi = []
    wi = wbalken
    dt=0
    while (math.ceil(xi + ubalken * dt) < xgrenz) and (math.ceil(zi) < zgrenz):
        xold = xi
        zold = zi

        if (zi < 0):
            zi = -zi
            wi = -wi
            tl, dt = prandltl(zi)
            ui = prandl(zi)
            xi, wi, zi = positionen(xi, wi, zi, tl, ui, dt )  # ,dt)
            rangecheck(xi, xold, zi, zold)


        else:
            tl, dt = prandltl(zi)
            ui = prandl(zi)
            xi, wi, zi = positionen(xi, wi, zi, tl, ui, dx)  # ,dt)
            # gitweis(xi,zi)
            rangecheck(xi, xold, zi, zold)

print(np.max(konk))
print(np.max(gitter))

#konzentrationen=[i * ((q * dt)/(n * dx * dz)) for i in gitter]
d = nc4.Dataset('alla.nc', 'w', format='NETCDF4')  # 'w' stands for write
d.createDimension('x', xgrenz)
d.createDimension('z',  zgrenz)
cnet = d.createVariable('c', 'f8', ('z', 'x'))
xnet = d.createVariable('x', 'f8', 'x')
znet = d.createVariable('z', 'f8', 'z')
cnet[:] = np.transpose(konk)
znet[:] = np.arange(0., zgrenz, dx, dtype=float)
xnet[:] = np.arange(0., xgrenz, dz, dtype=float)
d.close()
grass()