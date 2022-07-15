using NetCDF
using Random, Distributions
using ProgressBars
using LinearAlgebra
using Plots
#ncinfo("Übung_5/input_uebung5.nc")
#readdir()
#x = ncread("a.nc", "c")
n = 10^0 # !Anzahl Partikel
ubalken = 5  # !m/s
wbalken =0  # !m/s
xq = 60.5  # !m
zq = 1  # !m
counter = 0
xgrenz = 120  # !m
zgrenz = 120

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

#ncinfo("Übung_5/input_uebung5.nc")
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
        dt = ((k * ustern) / sigw ^ 2) * abs(2)
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
    #print(xi)
    rangex = convert(Int64,floor(xi - xold))
    rangez = convert(Int64,floor(zi - zold))
    #xsi=0
    #zsi=0
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
            #println(zsi)
            push!(toks,(zsi - zold) / (zi - zold))
        end
    end
    tku = sort!(toks)
    #println(tku)
    for i in 2:length(tku)+1
        ti = tku[i]
        told = tku[i - 1]
        t = mean([told, ti])
        posx, posz = gg(xold, zold, xi, zi, t)
        #println(posz)
        gitter[posx, abs(posz)] += ((tku[i] - tku[i-1] )* dt)
        konk[abs(posx), abs(posz)] += ((tku[i] - tku[i-1])* dt*((q * dt)/(n * dx * dz)))
    return
       # end
    end
end

function positionen(xi, wi, zi, tl, ui, dt )
    rl = exp(- dt / tl)
    d = Normal()
    rr = rand(d, 1)[1]
    ui= 5 + (1-rl)*tl*1
    #print(tl)
    xi = xi + ui * dt
    wi = rl * wi + sqrt((1 - rl ^ 2)) * sigw * rr +(1-rl)*tl*1
    zi = zi + wi * dt
    return xi, wi, zi
end

function reflex(xi,zi)
    
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
    print(marongus[61,2])
    xi = xq
    zi = zq
    posi = []
    wi = wbalken
    dt=0
    while (ceil(xi + ubalken * dt) < xgrenz) & (ceil(zi) < zgrenz)
        xold = xi
        zold = zi
        #println(zi)
        if (zi < znull)
            difz= abs(znull-zi)
            zi = zi +2*difz
            wi = -wi
            #println(zi)
            tl, dt = prandltl(zi,xi)
            ui = prandl(zi)
            xi, wi, zi = positionen(xi, wi, zi, tl, ui, dt ) 
            rangecheck(xi, xold, zi, zold,dt)
            push!(posi, [xi,zi])

        else
            tl, dt = prandltl(zi,xi)
            ui = prandl(zi)
            xi, wi, zi = positionen(xi, wi, zi, tl, ui, dt)
            rangecheck(xi, xold, zi, zold,dt)
            push!(posi, [xi,zi])
        end
    end
    push!(ges,posi)
end
"""
gr();
for i in tqdm(1:length(ges))
    x=[]
    z=[]
    for j in 1:length(ges[i])
        push!(x,(ges[i][j][1]))
        push!(z,(ges[i][j][2]))
    end
    plot(x,z)
 
end
savefig("alla.png")
"""
#print(maximum(konk))
if  isfile("test.nc") == true
    print("alla")
    rm("test.nc",force=true)
end

nccreate("test.nc", "c", "x", collect(0:xgrenz), "z", collect(0:zgrenz))
ncwrite(konk, "test.nc", "c")


print(marongu[gurkenlist[1]])
print(marongu[1,1])
