using NetCDF
using Random, Distributions
using ProgressBars
using PlotlyJS


global n,ubalken,wbalken,zq,xq,xgrenz,zgrenz,tl,nx,ny,nz,dx,dy,dz,ui,q::Int
global dt,sigu,sigw::Float64
global units::String
global gitter,cd,x,cdground::Array
n= 10^3#!Anzahl Partikel
ubalken = 5 #!m/s
wbalken = 0 #!m/s
zq = 45 #!m
xq = 51 #!m
counter=0
xgrenz= 2000# !m
zgrenz=400
tl = 100  #s Zeit
dt = 0.4 # Zeitschritt
sigu= 0 #m/s
sigw= 0.39 #m/s

dx = 1
dy = 1
dz = 1
ui=5


## Modellparameter
nx = 2000
nz = 400
dx = 1
dz = 1
q= 150 #540 kg/h also 5.4e+8mg/h und so 150000 

#sigu = 2.5 * ustern  # m/s
#sigw = 1.3 * ustern  # m/s
## Schichtung

## Array Initialisieren

nxx=nx+1 
nzz=nz+1
units = "g/m^3"
gitter=zeros(nx,nz)
cd= zeros(nxx,nzz)
x = range(0,nxx)
rl= exp(- dt/tl)

function gauss()
    for i in ProgressBar(1 : nxx)
        for j in 1:  nzz
            if x[i]-xq <= 0.0
                cd[i,j] = 0.0
            else
                dez= 2*sigw^2*tl*((i/ubalken)-tl +tl*exp(-(i)/(ubalken*tl)))
                cd[i,j]=(q/(sqrt(2* pi)*sqrt(dez)*ubalken) *(exp((-((j*dx)-zq)^2)/(2*dez)) +exp((-((j*dz)+zq)^2)/(2*dez))))
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



return xi, wi, zi

end




function monte()
    for i in ProgressBar(1:n+1)
        xi=xq
        zi=zq
        ui=ubalken
        wi=wbalken
        posi=[]
    
    
        while (ceil(xi+ui*dt) < xgrenz)
            if zi<1
                zi=-zi
                wi= -wi
                xi,wi,zi =positionen(xi,wi,zi)
                xm = abs(convert(Int64,round(xi)))+1
                zm = abs(convert(Int64,round(zi)))+1
                gitter[xm,zm] = gitter[xm,zm] + 1
    
    
            else
                xi,wi,zi =positionen(xi,wi,zi)
                xm = abs(convert(Int64,round(xi)))+1
                zm = abs(convert(Int64,round(zi)))+1
                if zm> zgrenz
                    zm=zgrenz
                end
                gitter[xm,zm] = gitter[xm,zm] + 1
    
            end
        end

    end
    return konzentrationen = [i * ((q * dt)/(n * dx * dz)) for i in gitter]
end

 function grafen()
    print("grafen geht los")
    xm=ncread("sens.nc","z")
    ym=ncread("sens.nc","x")
    zm=transpose(ncread("sens.nc","c"))
    savefig(plot([contour(x=collect(0:nzz),y=collect(0:nxx),z=transpose(cd),
contours_coloring="lines",
line_width=2,
#colorscale="Hot",
contours_start=0.1,
contours_end=0.5,
contours_size=0.1,
name="Gauss",
showlegend=true ,),

contour(x=xm,y=ym,z=zm,
contours_coloring="lines",
line_width=2,
colorscale="electric",
contours_start=0.1,
contours_end=0.5,
contours_size=0.1,
name="Montecarlo",
showlegend=true ,)
],
Layout(
    title="Konzentration (" *units * ") im Vergleich für N=50000",
    xaxis_title="x (m)",
    yaxis_title="z (m)",
    legend=attr(
        x=1,
        y=1,
        yanchor="bottom",
        xanchor="right",
        orientation="h"
    )
)), "Bericht/Bilder/1b50k.png")#,width=1920, height=1080)
 end

 function expo(konzentrationen)
    if  isfile("Bericht/boob.nc") == true
        rm("Bericht/boob.nc",force=true)
    end
    cd=transpose(konzentrationen)
    nccreate("Bericht/boob.nc", "c", "x", collect(0:nxx-1),  "z", collect(0:nzz-1))
    ncwrite( konzentrationen,"Bericht/boob.nc", "c")
 end



function main()
    gauss()
    konzentrationen=monte()
    #expo(konzentrationen)
    grafen()
    
end
main()
#print(argmax(cdground))
#println(cd[712,1])
#cdground=cd[cd[:,1] .== 0, :]
#print(cdground)






