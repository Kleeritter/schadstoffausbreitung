using NetCDF
using Random, Distributions
using ProgressBars
using SymPy
using PlotlyJS
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
dz = 1

##Randbedingungen
Q= 1.5e+5 #540 kg/h also 5.4e+8mg/h und so 150000 
ubalken= 5
h= 100
zq= 45
xq=51

## Schichtung
tl = 100  #s Zeit
sigw= 0.39 #m/s

## Array Initialisieren

nxx=xgrenz +2
nzz=zgrenz +2
cd= zeros(nxx,nzz)
x = zeros(nx +1)
function gauss()
    for i in ProgressBar(1 : nxx)
        for j in 1:  nzz
            if x[i]-xq <= 0.0
                cd[i,j] = 0.0
            else
                dez= 2*sigw^2*tl*((i/ubalken)-tl +tl*exp(-i/(ubalken*tl)))
                cd[i,j]=(Q/(sqrt(2* pi)*sqrt(dez)*ubalken) *(exp((-((j*dx)-zq)^2)/(2*dez)) +exp((-((j*dz)+zq)^2)/(2*dez))))
            end
        end
    end
end

gauss()


if  isfile("b1.nc") == true
    rm("b1.nc",force=true)
end

nccreate("b1.nc", "c", "x", collect(0:nxx),  "z", collect(0:nzz))
ncwrite(cd, "b1.nc", "c")

plot(contour(cd))