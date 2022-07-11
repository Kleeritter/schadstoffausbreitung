using NetCDF
using Random, Distributions
using ProgressBars
#ncinfo("Übung_5/input_uebung5.nc")
#readdir()
#x = ncread("a.nc", "c")
n = 10 # !Anzahl Partikel
ubalken = 5  # !m/s
wbalken =0  # !m/s
xq = 0  # !m
zq = 0.5  # !m
counter = 0
xgrenz = 110  # !m
zgrenz = 25

dx = 1
dy = 1
dz = 1

ustern = 0.35
k = 0.38
znull = 0.008
sigu = 2.5 * ustern  # m/s
sigw = 1.3 * ustern  # m/s
q = 0.1
gitter=zeros(xgrenz+1,500)
konk=zeros(xgrenz+1,500)


function gg(xold, zold, xi, zi, t)

    xg = xold + t * (xi - xold)
    zg = zold + t * (zi - zold)
    #println(floor(zg))
    #prInt64("der wert",xg, t, xold,xi)
    if xg>xgrenz
        xg=xgrenz
    end
    if floor(zg) <0
        print("olla")
        zg=0
    end
    #println(floor(zg))
    return convert(Int32, floor(xg)+1), convert(Int32, floor(zg)+1)
end

function prandl(zi)
    if zi < znull
        ubalken = 0
    else
        ubalken = (ustern / k) * log(abs(zi) / znull)
    end
    return ubalken
end
function prandltl(zi)
    tl = ((k * ustern) / sigw ^ 2) * abs(zi)
    if (0.1*tl)>((k * ustern) / sigw ^ 2) * abs(2) #falls dt kleiner als tl in 2 m Hoehe
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
        gitweis(xi, zi)
    else
        exaktgitter(xi, xold, zi, zold,dt)
    end
end

function exaktgitter(xi, xold, zi, zold,dt)
    ti = []
    tj = []
    toks = []
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
        println(posz)
        #gitter[posx, abs(posz)] += ((tku[i] - tku[i-1] )* dt)
        konk[posx, abs(posz)] += ((tku[i] - tku[i-1])* dt*((q * dt)/(n * dx * dz)))
    return
       # end
    end
end

function positionen(xi, wi, zi, tl, ui, dt )
    rl = exp(- dt / tl)
    d = Normal()
    rr = rand(d, 1)[1]
    xi = xi + ui * dt
    wi = rl * wi + sqrt((1 - rl ^ 2)) * sigw * rr
    zi = zi + wi * dt
    return xi, wi, zi
end

function gitweis(xi, zi)
    xm = convert(Int64,floor((xi)))
    zm = convert(Int64,floor((zi)))
    gitter[xm, zm] = gitter[xm, zm] + 1
    konk[xm, zm] += 1*((q * dt)/(n * dx * dz))
    return
end
for i in ProgressBar(0:n)
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
            tl, dt = prandltl(zi)
            ui = prandl(zi)
            xi, wi, zi = positionen(xi, wi, zi, tl, ui, dt ) 
            rangecheck(xi, xold, zi, zold,dt)


        else
            tl, dt = prandltl(zi)
            ui = prandl(zi)
            xi, wi, zi = positionen(xi, wi, zi, tl, ui, dt)
            rangecheck(xi, xold, zi, zold,dt)
        end
    end
end
