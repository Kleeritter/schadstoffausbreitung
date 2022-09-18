using NetCDF
using Random, Distributions
using ProgressBars
using LinearAlgebra
#using PlotlyJS
using CairoMakie
using FileIO
#using SymPy
#n = 10^3# !Anzahl Partikel

xgrenz = 120  # !m
zgrenz = 120
units = "g/m^3"     # Einheit fuer Graphen
#r = symbols("r")
xlist=[]
zlist=[]
dx = 1
dy = 1
dz = 1
ges=[]
k = 0.38
znull = 0.008
q = 0.2
gitter=zeros(xgrenz,zgrenz)
konk=zeros(xgrenz,zgrenz)

marongu = ncread("Bericht/input_uebung5.nc", "u")
marongw =ncread("Bericht/input_uebung5.nc","w")
marongus = ncread("Bericht/input_uebung5.nc","u2") 
marongws = ncread("Bericht/input_uebung5.nc","w2")
gurkenlist=findall(x->x==-9999.0,marongu)
gurkenlistw=findall(x->x==-9999.0,marongw)
for i in 1: length(gurkenlist)
marongu[gurkenlist[i]]=NaN
end 

for i in 1: length(gurkenlistw)
    marongw[gurkenlist[i]]=NaN
end

println(marongus[61,2])
function gg(xold, zold, xi, zi, t)

    xg = xold + t * (xi - xold)
    zg = zold + t * (zi - zold)
    if xg>xgrenz
        xg=xgrenz
    end
    if floor(zg) <1
        zg=1
    end
    return convert(Int64, floor(xg)), convert(Int64, floor(zg))
end

function prandltl(zi,xi)
    xii=convert(Int64,floor(xi))
    zii=convert(Int64,floor(zi))
    if floor(zi)==0
zii=1
    end

    tl= 0.05*((k*zii)/(1+k*(zii/5)))/(0.23*sqrt(marongus[xii+1,zii+1]+marongws[xii+1,zii+1]))
    if (0.1*tl)>0.05*((k*2)/(1+k*(2/5)))/(0.23*sqrt(marongus[xii+1,zii+1]+marongws[xii+1,zii+1])) #falls dt kleiner als tl in 2 m Hoehe
      dt = 0.1*tl
    else
        dt = 0.05*((k*2)/(1+k*(2/5)))/(0.23*sqrt(marongus[xii+1,zii+1]+marongws[xii+1,zii+1]))
    end
    return tl,  dt
end  

function rangecheck(xi, xold, zi, zold,dt,n)
    rangex = floor(xi - xold)
    rangez = floor(zi - zold)
    if (rangex + rangez) < 2
        gitweis(xi, zi,dt,n)
    else
        exaktgitter(xi, xold, zi, zold,dt)
    end
end

function exaktgitter(xi, xold, zi, zold,dt)
    ti = []
    toks = []
    rangex = convert(Int64,floor(xi - xold))
    rangez = convert(Int64,floor(zi - zold))
    xsi = ceil(xold)
    zsi = ceil(zold)
    for i in 0:rangex
        
        if i == 0
            xsi = ceil(xold)
            push!(toks,(xsi - xold) / (xi - xold))
        else
            xsi += 1
            push!(toks,(xsi - xold) / (xi - xold))
        end
    end
    for i in 0:rangez
        if i == 0
            zsi = ceil(zold)
            push!(toks,(zsi - zold) / (zi - zold))
        else
            zsi += 1

            push!(toks,(zsi - zold) / (zi - zold))
        end
    end
    tku = sort!(toks)
    for i in 2:length(tku)
        ti = tku[i]
        told = tku[i - 1]
        t = mean([told, ti])
        posx, posz = gg(xold, zold, xi, zi, t)
        gitter[posx, abs(posz)] += ((tku[i] - tku[i-1] )* dt)
        konk[abs(posx), abs(posz)] += ((tku[i] - tku[i-1])* dt*((q * dt)/(n * dx * dz)))
    return

    end
end

function positionen(xi, wi, zi, tl, ui, dt,xold,zold,xolder,zolder )

ixolder= floor(xolder )+1
izolder=floor(zolder )+1
    ixold= floor(xold )+1
    izold= floor(zold )+1
    izi= floor(zi )+1
    ixi= floor(xi)+1
    rl = exp(- dt / tl)
    Random.seed!()
    d = Normal()
    rr = rand(d, 1)[1]
    if izi == 1
        izi = 2
    end
    if ixi == 91
        ixi=90
    end
    if ixold== 91
        ixold=90
    end
    if izold== 1
        izold=2
    end

    difqu= abs(marongus[abs(convert(Int64,(ixolder))),abs(convert(Int64,(izolder)))]-marongus[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])/ 2* abs(xolder-xi)
    difqw=abs(marongws[abs(convert(Int64,(ixolder))),abs(convert(Int64,(izolder)))]-marongws[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])/ 2*abs(xolder-xi)

    
    ukack= rl *ui  + sqrt((1 - rl ^ 2)) * sqrt(marongus[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])*rr +(1-rl)*tl*difqu
    ui= marongu[abs(convert(Int64,(ixold))),abs(convert(Int64,(izold)))] +ukack
    xi = xi + ui * dt
    wkack = rl * wi + sqrt((1 - rl ^ 2)) *sqrt(marongws[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])*rr +(1-rl)*tl*difqw
    wi=marongw[abs(convert(Int64,(ixold))),abs(convert(Int64,(izold)))] +wkack
    zi = zi + wi * dt
    while ((xi<=31 && zi<=61)||(xi>=90 && zi <=61)||(30<=xi<=90 && zi<=1)||zi<1)
        if (xi>=90 && zi<=61) # rechte Wand
            br=1
            if br==1
                xi=xi-2*(abs(90-xi))
                ui=-ui
            end
        elseif (xi<=31 && zi<=61) && zi>1 #linke Wand
            br=1
            if br==1
            xi=xi+2*(abs(31-xi))
            ui=-ui
            end

        else  #Boden
            #println("Boden ist aus Lava")
            wi=-wi
            zi=zi+2*(abs(1-zi))
        end
    
end
return xi, wi, zi,ui

end

"""
function  eckendreck(xi,zi,xold,zold,ui,wi)
    if xi<=30&& typeof(solve(r*[1,1]+[30,60]-[xold,zold],r)) !=Vector{Any} && typeof(solve(r*[1,1]+[30,60]-[xi,zi],r)) !=Vector{Any}
        ui=-ui
        wi=-wi
        xi=xold
        zi=zold
        br=2

    elseif xi<=30&& typeof(solve(r*[1,1]+[30,0]-[xold,zold],r)) !=Vector{Any} && typeof(solve(r*[1,1]+[30,0]-[xi,zi],r)) !=Vector{Any}
            ui=-ui
            wi=-wi
            xi=xold
            zi=zold
            br=2
    elseif xi>=90 && typeof(solve(-r*[1,1]+[90,0]-[xold,zold],r))!=Vector{Any} && typeof(solve(-r*[1,1]+[90,0]-[xi,zi],r)) !=Vector{Any}
            ui=-ui
            wi=-wi
            xi=xold
            zi=zold
            br=2
    elseif xi>=90&& typeof(solve(-r*[1,1]+[90,60]-[xold,zold],r))!=Vector{Any} && typeof(solve(-r*[1,1]+[90,60]-[xi,zi],r)) !=Vector{Any}
        ui=-ui
        wi=-wi
        xi=xold
        zi=zold
        br=2
else
    ui=ui
    wi=wi
    br=1
end
return ui,wi, br
end 
"""

function gitweis(xi, zi,dt,n)
    if floor(zi) <0
        zi=0
    end
    xm = abs(convert(Int64,floor((xi +1))))
    zm = abs(convert(Int64,floor((zi +1))))
    gitter[xm, zm] = gitter[xm, zm] + 1
    konk[xm, zm] += 1*((q * dt)/(n * dx * dz))
    return
end

function monte(xq,zq,n)
for i in ProgressBar(1:n)
    xi = xq
    zi = zq
    dt=0
    ui=0
    wi=0
    xold=xq
    zold=zq
    while (ceil(xi+ui*dt) < xgrenz) 
        xolder=xold
        zolder=zold
        xold = xi
        zold = zi
            tl, dt = prandltl(zi,xi)
            xi, wi, zi,ui = positionen(xi, wi, zi, tl, ui, dt,xold,zold,xolder,zolder)

            rangecheck(xi, xold, zi, zold,dt,n)
            push!(xlist, xi)
            push!(zlist, zi)
    end

end
return  konk 
end



function vergleichmakie(xq,zq)
na= 10
nb= 100
nc =500
nd= 1000


 levels= [0,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2]#[0.001,0.005,0.01,0.05,0.1,0.5,0.75,1,2,5]#-1:0.1:1
 fig = Figure(resolution=(1080,1080))
 xs=LinRange(0, xgrenz,xgrenz)
 ys=LinRange(0, zgrenz, zgrenz)


contour(fig[1, 1],ys, xs, monte(xq,zq,na),levels=levels,label = "Red Dots")
contour(fig[1, 2],xs, ys, monte(xq,zq,nb),levels=levels,title= "N = " *string(nb))
contour(fig[2, 1],xs, ys, monte(xq,zq,nc),levels=levels,title= "N = " *string(nc))
contour(fig[2, 2],xs, ys, monte(xq,zq,nd),levels=levels,title= "N = " *string(nd))
Colorbar(fig[1:2,3], limits = (0.1, 2), colormap = :viridis,
    flipaxis = false)
axislegend()
Label(fig[0, :], text = "Partikelanzahlen im Vergleich fuer x = " *string(xq)*"m z = " *string(zq)*"m", textsize = 30)
save(
"Bericht/Bilder/3_vergleich_x = "*string(xq)*".png", fig)
end

function einzelmakie(xq,zq)
n=1000
levels= [0,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2]#-1:0.1:1
fig = Figure(resolution=(1080,1080))
xs=LinRange(0, xgrenz,xgrenz)
ys=LinRange(0, zgrenz, zgrenz)


contour(fig[1, 1],ys, xs, monte(xq,zq,n),levels=levels)
Colorbar(fig[1, 2], limits = (0.1, 2), colormap = :viridis,
    flipaxis = false)
save(
"Bericht/Bilder/3_single_x = "*string(xq)*".png", fig)
end 

    
    function main()
        ### Aufgbabe a
        xq = 60.5  # !m
        zq = 0.5 # !m
        vergleichmakie(xq,zq)
        einzelmakie(xq,zq)
        ### Aufgbabe b
        xq = 15.5  # !m
        zq =  65.5 # !m
        #vergleichmakie(xq,zq)
        #einzelmakie(xq,zq)
        
    end
    main()
    
