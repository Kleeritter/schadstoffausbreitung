using NetCDF
using Random, Distributions
using ProgressBars
using Plots
using SymPy
n = 10^2 # !Anzahl Partikel
xq = 61.5  # !m
zq = 2  # !m
counter = 0
xgrenz = 120  # !m
zgrenz = 120

## Modellparameter
nx = 1000
ny = 50
nz = 20
dx = 1
dy = 1
dz = 1

##Randbedingungen
Q= 1.5e+5 #540 kg/h also 5.4e+8mg/h und so 150000 
ubalken= 5
h= 100
zq= 100
xq=0
yq=250

## Schichtung
fg=0.504
fk=0.818
gg=0.265
gk=0.818


## Array Initialisieren
c=zeros(nx,ny,nz)


function main()
    for i in ProgressBar(0:nx)
        for j in 0:ny
            for k in 0:nz
                dey = fg*((dx*i)^fk)
                dez =gg*((dx*i)^gk)
                c(i,j,k) =(Q/(2* pi*dey*dez*ubalken))*exp((-((j*dy)-yq)^2)/(2*dey^2)) *(exp((-((k*dz)-h)^2)/(2*dez^2)) +exp((-((k*dz)+h)^2)/(2*dez^2)))
            end
        end
    
    end   
end
main()

print(max(c))

if  isfile("u1.nc") == true
    rm("u1.nc",force=true)
end

nccreate("u1.nc", "c", "x", collect(0:nx),  "y", collect(0:ny),"z", collect(0:nz))
ncwrite(c, "u1.nc", "c")

