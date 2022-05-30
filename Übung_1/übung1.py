import netCDF4 as ncf
import math
dx=10
dy=10
dz=10
nx = 1000
ny = 50
nz = 20
#Randbedingungen
Q= 1.5e+5 #540 kg/h also 5.4e+8mg/h und so 150000 
ubalken= 5
h= 100
zq= 100
xq=0
yq=250

c=[]    
#Schichtung 
fg=0.504
fk=0.818
gg=0.265
gk=0.818
def berechneconc() :
    for i in range(0,nx):
        for j in range(0,ny):
            for k in range(0,nz):
                dey = fg*((dx*i)**fk)
                print(dey)
       	        dez =gg*((dx*i)**gk)
                c[i,j,k] =(Q/(2* math.pi *dey*dez*ubalken))* math.exp((-((j*dy)-yq)**2)/(2*dey**2)) *(math.exp((-((k*dz)-h)**2)/(2*dez**2)) + math.exp((-((k*dz)+h)**2)/(2*dez**2)))

    return c

berechneconc()
print(c)