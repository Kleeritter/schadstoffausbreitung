using NetCDF
using Random, Distributions
#ncinfo("Übung_5/input_uebung5.nc")
readdir()
#x = ncread("a.nc", "c")
n = 10^4 # !Anzahl Partikel
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

function alla()
    d = Normal()
    rr = rand(d, 1)
    #print(5*rr)
    return rr    
end

function gg(xold, zold, xi, zi, t)

    xg = xold + t * (xi - xold)
    zg = zold + t * (zi - zold)
    #print("der wert",xg, t, xold,xi)
    if xg>xgrenz
        xg=xgrenz-1
    end
    return int(xg), int(zg)

function prandl(zi)
    if zi < znull
        ubalken = 0
    else
        ubalken = (ustern / k) * log(abs(zi) / znull)
    end
    return ubalken

function prandltl(zi)
    tl = ((k * ustern) / sigw ^ 2) * abs(zi)
    if (0.1*tl)>((k * ustern) / sigw ^ 2) * abs(2) #falls dt kleiner als tl in 2 m Hoehe
        dt = 0.1*tl
    else
        dt = ((k * ustern) / sigw ^ 2) * abs(2)
    end    
    return tl,  dt
k=alla()
print(k)

#d= InverseGaussian()
