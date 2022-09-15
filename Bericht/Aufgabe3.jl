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
#ustern = 0.35
k = 0.38
znull = 0.008
#sigu = 2.5 * ustern  # m/s
#sigw = 1.3 * ustern  # m/s
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
    #println(ukack, " mit ", ixi, " und ", izi, " und ", tl)
    ui= marongu[abs(convert(Int64,(ixold))),abs(convert(Int64,(izold)))] +ukack
    xi = xi + ui * dt
    #print(xi)
    wkack = rl * wi + sqrt((1 - rl ^ 2)) *sqrt(marongws[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])*rr +(1-rl)*tl*difqw
    wi=marongw[abs(convert(Int64,(ixold))),abs(convert(Int64,(izold)))] +wkack
    zi = zi + wi * dt

    #print(marongu[abs(convert(Int64,floor(xold))),abs(convert(Int64,floor(zold)))], " mit ",xold," zold: ",zold)
    #println(dt," alla ",tl," neues Glueck ", xi ," ", zi)
    #println(xi, " ", zi)
    while ((xi<=31 && zi<=61)||(xi>=90 && zi <=61)||(30<=xi<=90 && zi<=1)||zi<1)
        #print("irgendwas")
        if (xi>=90 && zi<=61) # rechte Wand
            #println("alarm rechte Wand")
            #ui,wi, br= eckendreck(xi,zi,xold,zold,ui,wi)
            br=1
            #println(br)
            if br==1
                xi=xi-2*(abs(90-xi))
                ui=-ui
            end
        elseif (xi<=31 && zi<=61) && zi>1 #linke Wand
            #println("alarm linke Wand")
            #ui,wi,br= eckendreck(xi,zi,xold,zold,ui,wi)
            br=1
            if br==1
            xi=xi+2*(abs(31-xi))
            ui=-ui
            end
        #elseif abs(zi-60) < abs(zi-1)  #Decke
         #   #println("alarm Decke")
          #  wi=-wi
           # zi=zi+2*(abs(61-zi))
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
    while (ceil(xi+ui*dt) < xgrenz) #& (ceil(zi) < zgrenz)
        xolder=xold
        zolder=zold
        xold = xi
        zold = zi
        #wi = marongw[abs(convert(Int64,floor(xold))),abs(convert(Int64,floor(zold)))]
        #ui= marongu[abs(convert(Int64,floor(xold))),abs(convert(Int64,floor(zold)))]
        #print(xi)
            tl, dt = prandltl(zi,xi)
            xi, wi, zi,ui = positionen(xi, wi, zi, tl, ui, dt,xold,zold,xolder,zolder)

            rangecheck(xi, xold, zi, zold,dt,n)
            #println(zi)
            push!(xlist, xi)
            push!(zlist, zi)
    end

end
return  konk #xlist, zlist,
end

#if  isfile("test.nc") == true
 #   rm("test.nc",force=true)
#end

#nccreate("test.nc", "c", "x", collect(0:xgrenz), "z", collect(0:zgrenz))
#ncwrite(konk, "test.nc", "c")

#plot(xlist,zlist, xlims=(0, 120), ylims=(0,120))

function einzelgrafen(konk,xq)
    quellenstring=string(xq)
    """
    savefig(plot(
    
    scatter(
    y= zlist,#collect(1:zgrenz),
    x=xlist,
    name="prairie_grass",
    ),
    
    
    Layout(
        title="Montecarlo in der Haeuserschlucht",
        xaxis_title="x(m)",
        yaxis_title="z(m)",
        xaxis_range=[0, 120], 
        yaxis_range=[0, 120]
    )), "Bericht/Bilder/3_xq="*quellenstring*".png")

"""
    savefig(plot(contour(x=collect(0:zgrenz),y=collect(0:xgrenz),z=transpose(konk),
    contours_coloring="lines",
    #showlabels=true,
    #labelfont = attr( # label font properties
     #       size = 12,
      #      color = "white",
       # ),
    #line_width=2,
    #colorscale="Hot",
    contours_start=0.0001,
    contours_end=5,
    contours_size=0.01,
    name="Gauss",),
    

    Layout(
        title="Konzentration (" *units * ") im Vergleich fuer N=50000",
        xaxis_title="x (m)",
        yaxis_title="z (m)",
        legend=attr(
            x=1,
            y=1,
            yanchor="bottom",
            xanchor="right",
            orientation="h"
        )
    )), "Bericht/Bilder/3k_xq="*quellenstring*".png")
     end
    
    
    function vergleichsgrafen(konkc,konkm,konkx,konkl,xq)
        quellenstring=string(xq)
        """
        make_plot(teilchen) = plot(contour(x=collect(0:zgrenz),y=collect(0:xgrenz),z=transpose(teilchen),
        contours_coloring="lines",
        contours_start=0.005,
        contours_end=1,
        contours_size=0.01,),
        )
        savefig([
            make_plot(konkc) make_plot(konkm)
            make_plot(konkx) make_plot(konkl)
        ],   
 "Bericht/Bilder/3_vergleich="*quellenstring*".png")
 """


 p = make_subplots(
    rows=2, cols=2, column_titles=["Without Smoothing"; "With Smoothing"]
)
#add_trace!(p, contour(x=collect(0:zgrenz),y=collect(0:xgrenz),z=transpose(konkc),#contours_coloring="lines",
#contours_start=0.005,
#contours_end=1,
#contours_size=0.01, ), row=1, col=1)
#add_trace!(p, contour(x=collect(0:zgrenz),y=collect(0:xgrenz),z=transpose(konkm),#contours_coloring="lines",
#contours_start=0.005,
#contours_end=1,
#contours_size=0.01, ), row=1, col=2)   

#add_trace!(p, contour(x=collect(0:zgrenz),y=collect(0:xgrenz),z=transpose(konkx), #contours_coloring="lines",
#contours_start=0.005,
#contours_end=1,
#contours_size=0.01,), row=2, col=1)

pa= plot(contour(x=collect(0:zgrenz),y=collect(0:xgrenz),z=transpose(konkc),contours_coloring="lines",
contours_start=0.005,
contours_end=1,
contours_size=0.01, ))

pb= plot(contour(x=collect(0:zgrenz),y=collect(0:xgrenz),z=transpose(konkm),contours_coloring="lines",
contours_start=0.005,
contours_end=1,
contours_size=0.01, ))

pups=[pa pb]
relayout!(pups, title_text="Side by side layout (1 x 2)")
pups
savefig(
pups

,"Bericht/Bilder/tester.png")

#add_trace!(p, contour(x=collect(0:zgrenz),y=collect(0:xgrenz),z=transpose(konkl),contours_coloring="lines",
#contours_start=0.005,
#contours_end=1,
#contours_size=0.01,), row=2, col=2)
#add_trace!(p, scatter(x=[1, 2], y=[1, 2]), row=2, col=2)
#relayout!(p, title_text="Side by side layout (1 x 2)")
#savefig(p, "Bericht/Bilder/3_vergleich="*quellenstring*".png")


"""
p = make_subplots(
    rows=2, cols=2,
    specs=[Spec() Spec(); Spec(colspan=2) missing],
    subplot_titles=["First Subplot" "Second Subplot"; "Third Subplot" missing]
)

add_trace!(p, scatter(x=[1, 2], y=[1, 2]), row=1, col=1)
add_trace!(p, scatter(x=[1, 2], y=[1, 2]), row=1, col=2)
add_trace!(p, scatter(x=[1, 2, 3], y=[2, 1, 2]), row=2, col=1)

relayout!(p, showlegend=false, title_text="Specs with Subplot Title")
savefig(p, "Bericht/Bilder/3_vergleich="*quellenstring*".png")

"""
end

function maggus(konk)
    levels= [0.01,0.05,0.75,1,5]#-1:0.1:1
    xs=LinRange(0, xgrenz,xgrenz)
    ys=LinRange(0, zgrenz, zgrenz)
     a=contour(xs, ys, konk,levels=levels)
    return a
end 
function vergleichmakie(xq,zq)
na= 100
nb= 1000
nc =5000
nd= 10000


 levels= [0.001,0.005,0.01,0.05,0.1,0.5,0.75,1,2,5]#-1:0.1:1
 fig = Figure(resolution=(1080,1080))
 xs=LinRange(0, xgrenz,xgrenz)
 ys=LinRange(0, zgrenz, zgrenz)


contour(fig[1, 1],ys, xs, monte(xq,zq,na),levels=levels)
contour(fig[1, 2],xs, ys, monte(xq,zq,nb),levels=levels,title= "N = " *string(nb))
contour(fig[2, 1],xs, ys, monte(xq,zq,nc),levels=levels,title= "N = " *string(nc))
contour(fig[2, 2],xs, ys, monte(xq,zq,nd),levels=levels,title= "N = " *string(nd))


Label(fig[0, :], text = "Partikelanzahlen im Vergleich für x = " *string(xq)*"m z = " *string(zq)*"m", textsize = 30)
save(
"Bericht/Bilder/3_vergleich x = "*string(xq)*".png", fig)
end


    
    function main()
        xq = 60.5  # !m
        zq = 0.5 # !m
        vergleichmakie(xq,zq)

        #konkl=monte(xq,zq,nl)
        #print(konkc)

        #einzelgrafen(konkm,xq)
        #vergleichsgrafen(konkc,konkm,konkx,konkl,xq)
        xq = 15.5  # !m
        zq =  65.5 # !m
        #konkc=monte(xq,zq,nc)
        #konkm=monte(xq,zq,nm)
        #konkx=monte(xq,zq,nx)
        #konkl=monte(xq,zq,nl)
        #einzelgrafen(konkm,xq)
        #vergleichsgrafen(konkc,konkm,konkx,konkl,xq)
        vergleichmakie(xq,zq)
    end
    main()
    
