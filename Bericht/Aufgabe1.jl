using NetCDF
using Random, Distributions
using ProgressBars
using PlotlyJS


## Modellparameter
nx = 2000
nz = 400
dx = 1
dz = 1

##Randbedingungen
Q= 150 #540 kg/h also 5.4e+8mg/h und so 150000 
ubalken= 5
zq= 45
xq=51

## Schichtung
tl = 100  #s Zeit
sigw= 0.39 #m/s

## Array Initialisieren

nxx=nx+1 
nzz=nz+1
units = "g/m^3"
cd= zeros(nxx,nzz)
x = range(0,nxx)
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
end

gauss()
cdground=[]
for j in 1:nxx
    #print(cd[j,1])
    push!(cdground,cd[j,1])

end
cdground=convert(Array{Float64,1}, cdground)

#print(argmax(cdground))
println(cd[712,1])
#cdground=cd[cd[:,1] .== 0, :]
#print(cdground)

if  isfile("Bericht/b1.nc") == true
    rm("Bericht/b1.nc",force=true)
end

nccreate("Bericht/b1.nc", "c", "x", collect(0:nxx),  "z", collect(0:nzz))
ncwrite(cd, "Bericht/b1.nc", "c")

savefig(plot(contour(x=collect(0:nzz),y=collect(0:nxx),z=transpose(cd),
contours_coloring="lines",
line_width=2,
colorscale="electric",
contours_start=0.1,
contours_end=0.5,
contours_size=0.1,),
Layout(
    title="Konzentration (" *units * ")",
    xaxis_title="x (m)",
    yaxis_title="z (m)",
)), "Bericht/1a.png")#,width=1920, height=1080)




savefig(plot(scatter(x=collect(0:nxx),y=cdground),
Layout(
    title="Konzentration (" *units * ") am Erdboden",
    xaxis_title="x (m)",
    yaxis_title="Konzentration (" *units * ")",
)), "Bericht/1b.png")#,width=1920, height=1080)
