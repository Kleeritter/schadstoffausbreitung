using NetCDF
using Random, Distributions
using ProgressBars
using LinearAlgebra
using PlotlyJS
using GLMakie
using FileIO

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
#znull = 0.008
q = 1

gitter=zeros(xgrenz,zgrenz)
konk=zeros(xgrenz,zgrenz)
function gg(xold, zold, xi, zi, t)

    xg = xold + t * abs(xi - xold)
    zg = zold + t * abs(zi - zold)
    if xg>xgrenz
        xg=xgrenz
    end
    if zg>zgrenz
        zg=zgrenz
    end
    if floor(zg) <1
        zg=1
    end
    return convert(Int64, floor(xg)), convert(Int64, floor(zg))
end

function prandltl(zi,xi,lambda,marongus,marongws,k)
    xii=convert(Int64,floor(xi))
    zii=convert(Int64,floor(zi))
    if floor(zi)==0
zii=1
    end
    tl= 0.05*((k*zii)/(1+k*(zii/lambda)))/(0.23*sqrt(marongus[xii+1,zii+1]+marongws[xii+1,zii+1]))
    if (0.1*tl)>0.05*((k*2)/(1+k*(2/lambda)))/(0.23*sqrt(marongus[xii+1,zii+1]+marongws[xii+1,zii+1])) #falls dt kleiner als tl in 2 m Hoehe
      dt = 0.1*tl
    else
        dt = 0.05*((k*2)/(1+k*(2/lambda)))/(0.23*sqrt(marongus[xii+1,zii+1]+marongws[xii+1,zii+1]))
    end
    return tl,  dt
end  
"""
function rangecheck(xi, xold, zi, zold,dt,n)
    rangex = floor(xi - xold)
    rangez = floor(zi - zold)
    if (rangex + rangez) < 2
        gitweis(xi, zi,dt,n)
    else
        exaktgitter(xi, xold, zi, zold,dt,n)
    end
end
"""
function exaktgitter(xi, xold, zi, zold,dt,n)
    ti = []
    toks = []
    rangex = convert(Int64,floor(abs(xi - xold)))
    rangez = convert(Int64,floor(abs(zi - zold)))

    if rangex >0 && rangez>0
        println("ja")

        for i in 0:rangex
            if xold >xi 
            xsi=ceil(xold) +i
            else
            xsi=ceil(xold) -i
            end
            push!(toks,(xsi - xold) / abs(xi - xold))
        end


        for i in 0:rangez
            if zold >zi 
                zsi=ceil(zold) +i
                else
                zsi=ceil(zold) -i
                end
            push!(toks,(zsi - zold) / abs(zi - zold))

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
else
    konk[convert(Int64,ceil(xi) ),convert(Int64,ceil(zi))] +=  ((q * dt)/(n * dx * dz))
end

end

function positionen(xi, wi, zi, tl, ui, dt,xold,zold,xolder,zolder,marongus,marongws,marongu,marongw )

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
    difqw=abs(marongws[abs(convert(Int64,(ixolder))),abs(convert(Int64,(izolder)))]-marongws[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])/ 2*abs(zolder-zi)

    
    ukack= rl *ui  + sqrt((1 - rl ^ 2)) * sqrt(marongus[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])*rr +(1-rl)*tl*difqu
    ui= marongu[abs(convert(Int64,(ixold))),abs(convert(Int64,(izold)))] +ukack
    xi = xi + ui * dt
    wkack = rl * wi + sqrt((1 - rl ^ 2)) *sqrt(marongws[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])*rr +(1-rl)*tl*difqw
    wi=marongw[abs(convert(Int64,(ixold))),abs(convert(Int64,(izold)))] +wkack
    zi = zi + wi * dt

    while ((xi<=31 && zi<=61)||(xi>=90 && zi <=61)||(30<=xi<=90 && zi<=0))
        if (xi>=90 && zi<=61) # rechte Wand
                xi= xi-2*(abs(90-xi))
                ui=-ui

        elseif (xi<=31 && zi<=61)  #linke Wand

            xi=xi+2*(abs(31-xi))
            ui=-ui

        elseif (30<=xi<=90 && zi<=0)  #Boden
            #println("Boden ist aus Lava")
            wi=-wi
            zi= -zi#zi+2*(abs(0-zi))
        else
            print("Hoppala")
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

function monte(xq,zq,n,lambda,sigma)
    #using SymPy
#n = 10^3# !Anzahl Partikel



marongu =  ncread("Bericht/input_uebung5.nc", "u")
marongw =ncread("Bericht/input_uebung5.nc","w")
marongus = sigma*ncread("Bericht/input_uebung5.nc","u2") 
marongws = sigma*ncread("Bericht/input_uebung5.nc","w2")
gurkenlist=findall(x->x==-9999.0,marongu)
gurkenlistw=findall(x->x==-9999.0,marongw)
for i in 1: length(gurkenlist)
marongu[gurkenlist[i]]=NaN
end 

for i in 1: length(gurkenlistw)
    marongw[gurkenlist[i]]=NaN
end
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
        xold = copy(xi)
        zold = zi

            tl, dt = prandltl(zi,xi,lambda,marongus,marongws,k)
            xi, wi, zi,ui = positionen(xi, wi, zi, tl, ui, dt,xold,zold,xolder,zolder,marongus,marongws,marongu,marongw)

            exaktgitter(xi, xold, zi, zold,dt,n)
            push!(xlist, xi)
            push!(zlist, zi)
    end

end
return  konk 
end



function vergleichmakie(xq,zq,lambda)
na= 10
nb= 100
nc = 100
nd= 100
sigma=1


levels=  [0.001,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2,5]#-1:0.1:1# [0,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2]#[0.001,0.005,0.01,0.05,0.1,0.5,0.75,1,2,5]#-1:0.1:1
 fig = Figure(resolution=(1080,1080))
 xs=LinRange(0, xgrenz,xgrenz)
 ys=LinRange(0, zgrenz, zgrenz)


 
a=monte(xq,zq,na,lambda,sigma)
contourf(fig[1, 1],ys, xs, a,levels=levels)
contour!(fig[1, 1],ys, xs,a ,levels=levels )
b=monte(xq,zq,nb,lambda,sigma)
contourf(fig[1, 2],xs, ys, b,levels=levels,title= "N = " *string(nb))
contour!(fig[1, 2],ys, xs,b ,levels=levels )
c=monte(xq,zq,nc,lambda,sigma)
contourf(fig[2, 1],xs, ys,c,levels=levels,title= "N = " *string(nc))
contour!(fig[2, 1],ys, xs,c ,levels=levels )
d=monte(xq,zq,nd,lambda,sigma)
contourf(fig[2, 2],xs, ys,d,levels=levels,title= "N = " *string(nd))
contour!(fig[2, 2],ys, xs,d ,levels=levels )



Colorbar(fig[1:2,3], limits = (0.1, 5), colormap = :viridis,
    flipaxis = false)

Label(fig[3, :], text = "x(m)")
Label(fig[:, 0], text = "z(m)")
Label(fig[0, :], text = "Konzentrationen im Vergleich fuer x = " *string(xq)*"m z = " *string(zq)*"m", textsize = 30)
save(
"Bericht/Bilder/3_vergleich_x = "*string(xq)*".png", fig)
end

function einzelmakie(xq,zq,lambda)
n=100
levels= [0.00001,0.01,1.0,2,10,100,200,500]#-1:0.1:1
fig = Figure(resolution=(1080,1080))
xs=LinRange(0, xgrenz,xgrenz)
ys=LinRange(0, zgrenz, zgrenz)

a=monte(xq,zq,n,lambda,1)
#al=maximum(a)
#println(al)
contourf(fig[1, 1],ys, xs,a ,levels=levels, xlabel = "x label", ylabel = "y label" )
contour!(fig[1, 1],ys, xs,a ,levels=levels )
Colorbar(fig[1,2], limits = (0, 5), colormap = :viridis,
    flipaxis = false,  label = "g/m3")

Label(fig[0, :], text = "Konzentrationen im Vergleich x = " *string(xq)*"m z = " *string(zq)*"m λ="*string(lambda), textsize = 30)
Label(fig[2, :], text = "x(m)")
Label(fig[:, 0], text = "z(m)")
save(
"Bericht/Bilder/3_single_x = "*string(xq)*"_"*string(lambda)*".png", fig)
end 

function lambdamakie(xq,zq)
    xgrenz = 120  # !m
    zgrenz = 120
    la= 10
    lb= 5
    lc =2
    ld= 1
    n= 1000
    sigma=5
    levels= [0.00001,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2]#-1:0.1:1# [0,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2]#[0.001,0.005,0.01,0.05,0.1,0.5,0.75,1,2,5]#-1:0.1:1
     fig = Figure(resolution=(1080,1080))
     xs=LinRange(0, xgrenz,xgrenz)
     ys=LinRange(0, zgrenz, zgrenz)
    a=monte(xq,zq,n,la,sigma)
    contourf(fig[1, 1],ys, xs, a,levels=levels)
    contour!(fig[1, 1],ys, xs,a ,levels=levels )
    b=monte(xq,zq,n,lb,sigma)
    contourf(fig[1, 2],xs, ys, b,levels=levels)
    contour!(fig[1, 2],ys, xs,b ,levels=levels )
    c=monte(xq,zq,n,lc,sigma)
    contourf(fig[2, 1],xs, ys,c,levels=levels)
    contour!(fig[2, 1],ys, xs,c ,levels=levels )
    d=monte(xq,zq,n,ld,sigma)
    contourf(fig[2, 2],xs, ys,d,levels=levels)
    contour!(fig[2, 2],ys, xs,d ,levels=levels )
    Colorbar(fig[1:2,3], limits = (0.00001, 500), colormap = :viridis,
        flipaxis = false, label = "g/m3")
    Label(fig[3, :], text = "x(m)")
    Label(fig[:, 0], text = "z(m)")
    Label(fig[0, :], text = "Konzentrationen im Vergleich x = " *string(xq)*"m z = " *string(zq)*" m bei verschiedenen λ", textsize = 30)
    save(
    "Bericht/Bilder/3_lambda_x = "*string(xq)*".png", fig)
end 

function sigmamakie(xq,zq)
    xgrenz = 120  # !m
    zgrenz = 120
    l=5
    n= 1000
    siga=1
    sigb=2
    sigc=5
    sigd=10
    levels=[0.001,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2,5]#-1:0.1:1# [0,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2]#[0.001,0.005,0.01,0.05,0.1,0.5,0.75,1,2,5]#-1:0.1:1
    fig = Figure(resolution=(1080,1080))
    xs=LinRange(0, xgrenz,xgrenz)
    ys=LinRange(0, zgrenz, zgrenz)
    a=monte(xq,zq,n,l,siga)
    contourf(fig[1, 1],ys, xs, a,levels=levels)
    contour!(fig[1, 1],ys, xs,a ,levels=levels )
    b=monte(xq,zq,n,l,sigb)
    contourf(fig[1, 2],xs, ys, b,levels=levels)
    contour!(fig[1, 2],ys, xs,b ,levels=levels )
    c=monte(xq,zq,n,l,sigc)
    contourf(fig[2, 1],xs, ys,c,levels=levels)
    contour!(fig[2, 1],ys, xs,c ,levels=levels )
    d=monte(xq,zq,n,l,sigd)
    contourf(fig[2, 2],xs, ys,d,levels=levels)
    contour!(fig[2, 2],ys, xs,d ,levels=levels )
    Colorbar(fig[1:2,3], limits = (0.00001, 5), colormap = :viridis,
        flipaxis = false, label = "g/m3")
    Label(fig[3, :], text = "x(m)")
    Label(fig[:, 0], text = "z(m)")
    Label(fig[0, :], text = "Konzentrationen im Vergleich x = " *string(xq)*"m z = " *string(zq)*" m bei verschiedenen σu & σw", textsize = 30)
    save(
    "Bericht/Bilder/3_sigma_x = "*string(xq)*".png", fig)
end 
    
    function main()

        ### Aufgbabe a
        xq = 60.5  # !m
        zq = 0.5 # !m
        #vergleichmakie(xq,zq,5)
        #lambdamakie(xq,zq)
        sigmamakie(xq,zq)
        ### Aufgbabe b
        xq = 15.5  # !m
        zq =  65.5 # !m
        #vergleichmakie(xq,zq,5)
        #lambdamakie(xq,zq)
        #sigmamakie(zq,xq)

        
    end
    main()
    