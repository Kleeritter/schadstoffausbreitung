using NetCDF
using Random, Distributions
using ProgressBars
using PlotlyJS


n= 30000#!Anzahl Partikel
ubalken = 5 #!m/s
wbalken = 0 #!m/s
zq = 45 #!m
xq = 51 #!m
counter=0
xgrenz= 2000# !m
zgrenz=400
ges=[]
tl = 100  #s Zeit
dt = 0.4 # Zeitschritt
sigu= 0 #m/s
sigw= 0.39 #m/s
rl= exp(- dt/tl)
nx = 1000
ny = 50
nz = 2000
dx = 2
dy = 2
dz = 2
ui=5


## Modellparameter
nx = 2000
nz = 400
dx = 1
dz = 1
ges=[]
##Randbedingungen
Q= 150 #540 kg/h also 5.4e+8mg/h und so 150000 
ubalken= 5
zq= 45
xq=51
#sigu = 2.5 * ustern  # m/s
#sigw = 1.3 * ustern  # m/s
## Schichtung
tl = 100  #s Zeit
sigw= 0.39 #m/s

## Array Initialisieren

nxx=nx+1 
nzz=nz+1
units = "g/m^3"
gitter=zeros(nxx,nzz)
cd= zeros(nxx,nzz)
x = range(0,nxx)
positionsliste=[]
xvalues=[]
zvalues=[]

function gauss()
    for i in ProgressBar(1 : nxx)
        for j in 1:  nzz
            if x[i]-xq <= 0.0
                cd[i,j] = 0.0
            else
                dez= 2*sigw^2*tl*((i/ubalken)-tl +tl*exp(-(i)/(ubalken*tl)))
                cd[i,j]=(Q/(sqrt(2* pi)*sqrt(dez)*ubalken) *(exp((-((j*dx)-zq)^2)/(2*dez)) +exp((-((j*dz)+zq)^2)/(2*dez))))
            end
        end
    end
    cdground=[]
    for j in 1:nxx
        #print(cd[j,1])
        push!(cdground,cd[j,1])

    end
    cdground=convert(Array{Float64,1}, cdground)
end


function positionen(xi, wi, zi)

    rl = exp(- dt / tl)
    Random.seed!()
    d = Normal()
    rr = rand(d, 1)[1]

    xi = xi + ui * dt
    wi=rl*wi + sqrt((1 - rl^2))*sigw* rr
    zi = zi + wi * dt


return xi, wi, zi,ui

end




function monte()
    for i in ProgressBar(1:n+1)
        xi=xq
        zi=zq
        ui=ubalken
        wi=wbalken
        posi=[]
    
    
        while xi<= xgrenz
            if zi<0
                zi=-zi
                wi= -wi
                xi,wi,zi =positionen(xi,wi,zi)
                push!(posi,[xi,zi])
                push!( xvalues,xi)
                push!(zvalues, zi)
                push!(positionsliste,[xi,zi])
                xm = int(xi)
                zm = int(zi)
                gitter[xm,zm] = gitter[xm,zm] + 1
    
    
            else
                xi,wi,zi =positionen(xi,wi,zi)
                push!(posi,[xi,zi])
                push!( xvalues,xi)
                push!(zvalues, zi)
                push!(positionsliste,[xi,zi])
                xm = int(xi)
                zm = int(zi)
                gitter[xm,zm] = gitter[xm,zm] + 1
    
            end
        end

    end
    konzentrationen = [i * ((q * dt)/(n * dx * dz)) for i in gitter]
    print(konzentrationen)
    print("alla")
end

 
function main()
    gauss()
    monte()
    #expo()
    #grafen()
    
end
#print(argmax(cdground))
#println(cd[712,1])
#cdground=cd[cd[:,1] .== 0, :]
#print(cdground)

if  isfile("Bericht/b1.nc") == true
    rm("Bericht/b1.nc",force=true)
end

nccreate("Bericht/b1.nc", "c", "x", collect(0:nxx),  "z", collect(0:nzz))
ncwrite(cd, "Bericht/b1.nc", "c")



