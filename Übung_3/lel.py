#

# Pakete

import random
import math as m
import matplotlib.pyplot as plt
import numpy as np
import os
from netCDF4 import Dataset

# Randbedingungen
Dx = 2
Dz = 2
tau_L = 100
Dt = 0.4
sig_u = 0
sig_w = 0.39
u_mittel = 5
w_mittel = 0
z_q = 45
x_q = 51
x_grenz = 2000
z_grenz = 400
Q = 150
N = 20000

# Lagrangsche Autokorrelationsfunktion (Näherung)

R_L = m.exp(-Dt/tau_L)

# Turbulenter Anteil der Windgeschwindigkeits z-Komponente

def w_turbo(R_L, w, sig_w):
    w_turb = R_L * w + m.sqrt(1 - R_L**2) * sig_w * random.gauss(0,1)
    return w_turb

# Gitter aufspannen

gitter = np.zeros((int(x_grenz/2 + 1),int(z_grenz/2 + 1)))

#print(gitter.shape)

# Partikel einem Gittervolumenmittelpunkt zuordnen

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
    
# Dies ist leider nicht allgemeingültig, sondern basiert auf der Gitterweite von 2. 
# Für eine allgemeingültige Lösung hat mir leider die Zeit gefehlt.
# Ziel ist es vor der nächsten Übung mit dem Panda-Paket eine richtige Lösung zu finden.

# N Partikeltrajektorien

for k in range(0,N):
    pos = [x_q,z_q]
    x = []
    z = []
    w_turb = 0
    u_turb = 0
    while pos[0] <= x_grenz:
        w = w_mittel + w_turb
        w_turb = w_turbo(R_L, w, sig_w)
        
        pos[0] = pos[0] + (u_mittel + u_turb)*Dt
        pos[1] = pos[1] + (w_mittel + w_turb)*Dt

        if pos[1] <= 0:
            pos[1] = -pos[1]
            w = -w
        
        xi = pos[0]
        zi = pos[1]
        
        xm = int((zuordnen(xi) - 1) / Dx)
        zm = int((zuordnen(zi) - 1) / Dz)
        
        gitter[xm,zm] = gitter[xm,zm] + 1
        
        x.append(pos[0])
        z.append(pos[1])
        
    if (k % 100) == 0:
        print(k)
        
    plt.plot(x,z, linewidth = 0.1)
    

konzentrationen = [i * ((Q * Dt)/(N * Dx * Dz)) for i in gitter]

#print(konzentrationen)

# Abspeichern als nc-Datei

if os.path.exists('uebung3.nc') == True:
    try:
        os.remove('uebung3.nc')
    except OSError as e:
        print('Error')
    else:
        print('Deleted')
         
u3 = Dataset('uebung3.nc', 'w', format='NETCDF4')
print(u3.data_model)

x_pkt =  int(x_grenz/Dx +1)
z_pkt =  int(z_grenz/Dz +1)

u3.createDimension('x', x_pkt)
u3.createDimension('z', z_pkt)

xu3 = u3.createVariable('x', 'f4', 'x')
zu3 = u3.createVariable('z', 'f4', 'z')
cu3 = u3.createVariable('c', 'f4', ('z', 'x'))

xu3[:]=np.arange(0.,x_pkt*Dx,Dx,dtype=float)
zu3[:]=np.arange(0.,z_pkt*Dz,Dz,dtype=float)
cu3[:]=np.transpose(konzentrationen)

xu3.units = 'meters'
zu3.units = 'meters'
cu3.units = 'g/m^3'

u3.close()


# Plotanpassung

plt.title('Partikeltrajektorien')
plt.xlabel('Strecke [m]')
plt.ylabel('Höhe [m]')
plt.axis([0,2000,-10,350])
plt.axhline(y=0, color='red', linestyle='-', linewidth=1)
plt.savefig('Trajektorien.jpg', format='jpg',dpi=2000)
plt.show()
plt.close()
