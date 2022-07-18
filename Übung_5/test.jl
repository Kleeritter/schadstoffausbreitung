using NetCDF
using Random, Distributions
using ProgressBars
using LinearAlgebra
using Plots
using SymPy
n = 0 # !Anzahl Partikel
ubalken = 5  # !m/s
wbalken =0  # !m/s
xq = 60.5  # !m
zq = 1  # !m
counter = 0
xgrenz = 120  # !m
zgrenz = 120

r = symbols("r")
xlist=[]
zlist=[]
dx = 1
dy = 1
dz = 1
ges=[]
ustern = 0.35
k = 0.38
znull = 1
sigu = 2.5 * ustern  # m/s
sigw = 1.3 * ustern  # m/s
q = 5
gitter=zeros(xgrenz+1,zgrenz+1)
konk=zeros(xgrenz+1,zgrenz+1)

marongu = ncread("Übung_5/input_uebung5.nc", "u")
marongw = ncread("Übung_5/input_uebung5.nc", "w")
marongus = ncread("Übung_5/input_uebung5.nc", "u2")
marongws = ncread("Übung_5/input_uebung5.nc", "w2")
gurkenlist=findall(x->x==-9999.0,marongu)
gurkenlistw=findall(x->x==-9999.0,marongw)
for i in 1: length(gurkenlist)
marongu[gurkenlist[i]]=NaN
end 

for i in 1: length(gurkenlistw)
    marongw[gurkenlist[i]]=NaN
end

function gg(xold, zold, xi, zi, t)

    xg = xold + t * (xi - xold)
    zg = zold + t * (zi - zold)
    #println(floor(zg))
    #prInt64("der wert",xg, t, xold,xi)
    if xg>xgrenz
        xg=xgrenz
    end
    if floor(zg) <0
        zg=0
    end
    #println(floor(zg))
    return convert(Int64, floor(xg+1)), convert(Int64, floor(zg+1))
end

function prandl(zi)
    if zi < znull
        ubalken = 0
    else
        ubalken = (ustern / k) * log(abs(zi) / znull)
    end
    return ubalken
end
function prandltl(zi,xi)
    xii=convert(Int64,floor(xi))+1
    zii=convert(Int64,floor(zi))+1
    #println("indexus",xii,zii)
    #println("so ein scheiss", sqrt(marongus[xii,zii]+marongws[xii,zii]))
    tl= 0.05*((k*zii)/(1+k*(zii/5)))/(0.23*sqrt(marongus[xii,zii]+marongws[xii,zii]))
    #print(tl)
    #tl = ((k * ustern) / sigw ^ 2) * abs(zi)
    if (0.1*tl)>0.05*((k*2)/(1+k*(2/5)))/(0.23*sqrt(marongus[xii,2]+marongws[xii,2])) #falls dt kleiner als tl in 2 m Hoehe
        dt = 0.1*tl
    else
        dt = 0.05*((k*2)/(1+k*(2/5)))/(0.23*sqrt(marongus[xii,2]+marongws[xii,2]))
    end    
    return tl,  dt
end  

function rangecheck(xi, xold, zi, zold,dt)
    rangex = floor(xi - xold)
    rangez = floor(zi - zold)
    if (rangex + rangez) < 2
        gitweis(xi, zi,dt)
    else
        exaktgitter(xi, xold, zi, zold,dt)
    end
end

function exaktgitter(xi, xold, zi, zold,dt)
    ti = []
    tj = []
    toks = []
    #print(xi,xold)
    #print(xi)
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
    for i in 2:length(tku)+1
        ti = tku[i]
        told = tku[i - 1]
        t = mean([told, ti])
        posx, posz = gg(xold, zold, xi, zi, t)
        gitter[posx, abs(posz)] += ((tku[i] - tku[i-1] )* dt)
        konk[abs(posx), abs(posz)] += ((tku[i] - tku[i-1])* dt*((q * dt)/(n * dx * dz)))
    return

    end
end

function positionen(xi, wi, zi, tl, ui, dt,xold,zold )

    rl = exp(- dt / tl)
    d = Normal()
    rr = rand(d, 1)[1]
    if xi==xq
        difqu=0
        difqw=0
    else
    difqu= abs(marongus[abs(convert(Int64,floor(xold)))+1,abs(convert(Int64,floor(zold)))+1]-marongus[abs(convert(Int64,floor(xi)))+1,abs(convert(Int64,floor(zi)))+1])/ 2* abs(xold-xi)
    difqw=abs(marongws[abs(convert(Int64,floor(xold)))+1,abs(convert(Int64,floor(zold)))+1]-marongws[abs(convert(Int64,floor(xi)))+1,abs(convert(Int64,floor(zi)))+1])/ 2*abs(xold-xi)
    end
    #println(difqu," mit ",xold," zold: ",zold)
    
    ui= marongu[abs(convert(Int64,floor(xold)))+1,abs(convert(Int64,floor(zold)))+1] + (1-rl)*tl*difqu
    xi = xi + ui * dt
    wi = marongw[abs(convert(Int64,floor(xold)))+1,abs(convert(Int64,floor(zold)))+1]+rl * wi + sqrt((1 - rl ^ 2)) * sigw * rr +(1-rl)*tl*difqw
    zi = zi + wi * dt

    #print(marongu[abs(convert(Int64,floor(xold)))+1,abs(convert(Int64,floor(zold)))+1], " mit ",xold," zold: ",zold)
    println(marongu[abs(convert(Int64,floor(xold)))+1,abs(convert(Int64,floor(zold)))+1]," alla ",marongw[abs(convert(Int64,floor(xold)))+1,abs(convert(Int64,floor(zold)))+1]," neues Glück ", xi ," ", zi)
    while ((xi<=30 && zi<=60)||(xi>=90 && zi<=60)||(30<=xi<=90 && zi<=0)||zi<1)
        if xi>xold && zold<60 # rechte Wand
            println("alarm rechte Wand")
            ui,wi, br= eckendreck(xi,zi,xold,zold,ui,wi)
            if br==1
                xi=xi-2*(abs(90-xi))
                ui=-ui
            end
        elseif xi<xold && zold<60 && zi>1 #linke Wand
            println("alarm linke Wand")
            ui,wi,br= eckendreck(xi,zi,xold,zold,ui,wi)
            if br==1
            xi=xi+2*(abs(30-xi))
            ui=-ui
            end
        elseif abs(zi-60) < abs(zi-1)  #Decke
            println("alarm Decke")
            wi=-wi
            zi=zi+2*(abs(60-xi))
        else  #Boden
            println("Boden ist aus Lava")
            wi=-wi
            zi=zi+2*(abs(0-zi))
        end
    
end
#println(zi)
return xi, wi, zi,ui

end


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
            br=2
    elseif xi>=90 && typeof(solve(-r*[1,1]+[90,0]-[xold,zold],r))!=Vector{Any} && typeof(solve(-r*[1,1]+[90,0]-[xi,zi],r)) !=Vector{Any}
            ui=-ui
            wi=-wi
            br=2
    elseif xi>=90&& typeof(solve(-r*[1,1]+[90,60]-[xold,zold],r))!=Vector{Any} && typeof(solve(-r*[1,1]+[90,60]-[xi,zi],r)) !=Vector{Any}
        ui=-ui
        wi=-wi
        br=2
else
    ui=ui
    wi=wi
    br=1
end
return ui,wi, br
end 




function gitweis(xi, zi,dt)
    if floor(zi) <0
        zi=0
    end
    xm = abs(convert(Int64,floor((xi+1))))
    zm = abs(convert(Int64,floor((zi+1))))
    gitter[xm, zm] = gitter[xm, zm] + 1
    konk[xm, zm] += 1*((q * dt)/(n * dx * dz))
    return
end
for i in ProgressBar(0:n)
    xi = xq
    zi = zq
    posi = []
    dt=0
    ui=0
    while (ceil(xi + ui * dt) < xgrenz) & (ceil(zi) < zgrenz)
        xold = xi
        zold = zi
        wi = marongw[abs(convert(Int64,floor(xold)))+1,abs(convert(Int64,floor(zold)))+1]
        ui= marongu[abs(convert(Int64,floor(xold)))+1,abs(convert(Int64,floor(zold)))+1]

            tl, dt = prandltl(zi,xi)
            ui = prandl(zi)
            xi, wi, zi,ui = positionen(xi, wi, zi, tl, ui, dt,xold,zold)

            rangecheck(xi, xold, zi, zold,dt)
            #println(zi)
            push!(posi, [xi,zi])
            push!(xlist, xi)
            push!(zlist, zi)
    end
    push!(ges,posi)
end

if  isfile("test.nc") == true
    rm("test.nc",force=true)
end

nccreate("test.nc", "c", "x", collect(0:xgrenz), "z", collect(0:zgrenz))
ncwrite(konk, "test.nc", "c")

plot(xlist,zlist)