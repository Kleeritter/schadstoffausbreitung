import netCDF4 as ncf
from netCDF4 import Dataset
import math
import numpy as np
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

   
#Schichtung 
fg=0.504
fk=0.818
gg=0.265
gk=0.818
c= np.full([nx,ny,nz],0.0)
def berechneconc() :
    
    for i in range(1,(nx +1)):
        for j in range(1,(ny +1)):
            for k in range(1,(nz+1)):
                dey = fg*((dx*i)**fk)
                dez =gg*((dx*i)**gk) 
                c[(i-1),(j-1),(k-1)] =(Q/(2* math.pi *dey*dez*ubalken))* math.exp((-((j*dy)-yq)**2)/(2*dey**2)) *(math.exp((-((k*dz)-h)**2)/(2*dez**2)) + math.exp((-((k*dz)+h)**2)/(2*dez**2)))

    print(np.max(c))
    return c


def netcdfcreawrite() :
    rootgrp = Dataset("testpyer.nc", "w", format="NETCDF4")
    rootgrp.description = "bogus example script"
    #rootgrp.Conventions='COARDS'
    rootgrp.title='Übung1'
    x=rootgrp.createDimension("x",nx)
    y=rootgrp.createDimension("y",ny)
    z=rootgrp.createDimension("z",nz)
    con=rootgrp.createDimension("c",len(c))
    xvalues = rootgrp.createVariable("x","f8",("x"))
    yvalues = rootgrp.createVariable("y","f8",("y"))
    zvalues = rootgrp.createVariable("z","f8",("z",))
    cvalues = rootgrp.createVariable("c","f8",("c","x","y","z"))
    xvalues.units="meters"
    yvalues.units="meters"
    zvalues.units="meters"
    cvalues.units="mgm^-3"
    
    #c=rootgrp.createDimension("c",nz)
    for dimobj in rootgrp.dimensions.values():
         print(dimobj)
    xvalues=np.arange(0,nx)
    yvalues=np.arange(0,ny)
    zvalues=np.arange(0,nz)
    cvalues=c
   
    rootgrp.close()
    return

def netcdfwrite() : 
 
    return

berechneconc()
netcdfcreawrite()
