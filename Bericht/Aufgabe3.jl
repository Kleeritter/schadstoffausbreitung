using NetCDF
using Random, Distributions
using ProgressBars
using LinearAlgebra
using GLMakie
using SymPy
using FileIO

# Globale Variablen
xgrenz = 120        # Grenze des Modells in x Richtung in m
zgrenz = 120        # Grenze des Modells in z Richtung in m
r = symbols("r")    # Wird benötigt für die Eckenreflexion
dx = 1              # Gitterweite in x-Richtung in m
dz = 1              # Gitterweite in z-Richtung in m
k = 0.38            # Karman-Konstante
q = 1               # Anfangskonzentration in kg/s

## Arrays Initialisieren

#Array wird in der Größe xgrenz*zgrenz aufgespannt und mit 0 gefüllt
konk=zeros(xgrenz,zgrenz)
# Montecarlo-Modell #
##Funktionen zur Berechnung des Montecarlo Modell ##

### Berechnung der Prandtlschicht###
function prandltl(zi,xi,lambda,us,ws,k)
    xii=convert(Int64,floor(xi))
    zii=convert(Int64,floor(zi))
    if floor(zi)==0
    zii=1
    end
    #Berechnung des Lagrangeschen Zeitschrittes tl
    tl= 0.05*((k*zii)/(1+k*(zii/lambda)))/(0.23*sqrt(us[xii+1,zii+1]+ws[xii+1,zii+1]))
    
    #### Berechnung des Zeitschrittes dt
    if (0.1*tl)>0.05*((k*2)/(1+k*(2/lambda)))/(0.23*sqrt(us[xii+1,zii+1]+ws[xii+1,zii+1]))  #falls dt kleiner als tl in 2 m Hoehe
      dt = 0.1*tl                                                                                       # Normalfall
    else
        dt = 0.05*((k*2)/(1+k*(2/lambda)))/(0.23*sqrt(us[xii+1,zii+1]+ws[xii+1,zii+1]))     #falls dt kleiner als tl in 2 m Hoehe wird dt auf tl(2m) gesetzt
    end
    return tl,  dt
end  

### Gitterauswertung ###
function exaktgitter(xi, xold, zi, zold,dt,n)
    #Initialisierung
    ti = [] 
    toks = []
    # Berechnung der Anzahl an Gitterschnitpunkten
    rangex = convert(Int64,floor(abs(xi - xold)))
    rangez = convert(Int64,floor(abs(zi - zold)))

    # Fallunterscheidung: Gibt es Schnittpunkte
    if rangex >0 && rangez>0
    #Falls es sie gibt muss die Konzentration pro Anteil Giterpunkt berechnet werden
        # Berechnung von ti und tj
        for i in 0:rangex
            if xold >xi 
            xsi=ceil(xold) +i
            else
            xsi=ceil(xold) -i
            end
            push!(toks,(xsi - xold) / abs(xi - xold)) # Einhängen in die Liste
        end
        for i in 0:rangez
            if zold >zi 
                zsi=ceil(zold) +i
                else
                zsi=ceil(zold) -i
                end
            push!(toks,(zsi - zold) / abs(zi - zold)) # Einhängen in die Liste
        end
        tku = sort!(toks) #Sortieren der Liste
        for i in 2:length(tku)
            ti = tku[i]
            told = tku[i - 1]
            t = mean([told, ti])
            # Der Mittelwert t der beiden Punkte told und ti kann in einer Geradengleichung verwendet werden
            # So können die Betroffenen Gitterpunkte ermittelt werden
            posx, posz = gg(xold, zold, xi, zi, t) #Ausführen der Geradengleichung
            #Anschließend wird die konzentrationen am mit der GG berechneten Punkt um den entsprechenden Betrag erhöht
            konk[abs(posx), abs(posz)] += ((tku[i] - tku[i-1])* dt*((q * dt)/(n * dx * dz))) 
        end
else
    #Gitterauswertung falls Partikel in dem Giterpunkt verbleibt
    konk[convert(Int64,ceil(xi) ),convert(Int64,ceil(zi))] +=  ((q * dt)/(n * dx * dz)) #Erhöhung um die gesammte Konzentration
end
end

### Berechnung der Geradengleichung ###
function gg(xold, zold, xi, zi, t)
    xg = xold + t * abs(xi - xold)
    zg = zold + t * abs(zi - zold)
    
    #Durch die Geradengleichung können Gitterpunkte jenseits der Grenzen berechnet werden
    #Um einen Abbruch zu verhindern werden in disem Fall die Maximalwerte verwendet
    if xg>xgrenz
        xg=xgrenz
    end
    if zg>zgrenz
        zg=zgrenz
    end
    return convert(Int64, floor(xg)), convert(Int64, floor(zg)) #Zurückgeben der Werte als Int
end

### Funktion zur Berechnung der Positionen ###
function positionen(xi, wi, zi, tl, ui, dt,xold,zold,xolder,zolder,us,ws,u,w )
    #Da Julia mit einer Indexierung ab 1 startet müssen die Indexe der Parameter angepasst werden
    ixolder= floor(xolder )+1 
    izolder=floor(zolder )+1
    ixold= floor(xold )+1
    izold= floor(zold )+1
    izi= floor(zi )+1
    ixi= floor(xi)+1
    #Lagrangesche Autokorellationsfunktion
    rl = exp(- dt / tl)
    Random.seed!()
    d = Normal()
    rr = rand(d, 1)[1]
    # Abfangen von möglichen Grenzfällen für die Indexe
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

    # Zentraler Diffentialquotient
    difqu= abs(us[abs(convert(Int64,(ixolder))),abs(convert(Int64,(izolder)))]-us[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])/ 2* abs(xolder-xi)
    difqw=abs(ws[abs(convert(Int64,(ixolder))),abs(convert(Int64,(izolder)))]-ws[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])/ 2*abs(zolder-zi)

    ## Turbulenter Anteil   
    uturbo = rl * ui + sqrt((1 - rl ^ 2)) * sqrt(us[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])*rr +(1-rl)*tl*difqu
    wturbo = rl * wi + sqrt((1 - rl ^ 2)) * sqrt(ws[abs(convert(Int64,(ixi))),abs(convert(Int64,(izi)))])*rr +(1-rl)*tl*difqw

    #Berechnung der Windgeschwindigkeit (Mittelwert + Turbulenz)
    ui = u[abs(convert(Int64,(ixold))),abs(convert(Int64,(izold)))] +uturbo
    wi = w[abs(convert(Int64,(ixold))),abs(convert(Int64,(izold)))] +wturbo

    #Berechnung der neuen Positionen
    xi = xi + ui * dt
    zi = zi + wi * dt

    #Reflexion überprüfen
    while ((xi<=31 && zi<=61)||(xi>=90 && zi <=61)||(30<=xi<=90 && zi<=0))
        if (xi>=90 && zi<=61) # Reflexion an der rechten Wand
            eckenreflexion(xi,zi,xold,zold,ui,wi)
            if br==1 # Falls dies der Fall ist gab es keine Eckenreflexion
                xi= xi-2*(abs(90-xi))
                ui=-ui
            end
        elseif (xi<=31 && zi<=61)  # Reflexion an der linken Wand
            eckenreflexion(xi,zi,xold,zold,ui,wi)
            if br==1 # Falls dies der Fall ist gab es keine Eckenreflexion
                xi=xi+2*(abs(31-xi))
                ui=-ui
            end
        elseif (30<=xi<=90 && zi<=0)  #Reflexion am Boden
            eckenreflexion(xi,zi,xold,zold,ui,wi)
            if br==1 # Falls dies der Fall ist gab es keine Eckenreflexion
            wi=-wi
            zi= -zi
            end
        else
            print("Fehler in der Reflexion")
        end
    

end
return xi, wi, zi,ui

end


### Funktion zur Überprüfung der Eckenreflexion ###
function  eckenreflexion(xi,zi,xold,zold,ui,wi)
    #Hier wird die Reflexion an den Ecken überprüft.
    #Falls es eine Eckenreflexion gibt muss eine der folgenden Gleichungen lösbar sein, wodurch sich deren Typ verändert
    #Reflexion an der linken oberen Ecke
    if xi<=30&& typeof(solve(r*[1,1]+[30,60]-[xold,zold],r)) !=Vector{Any} && typeof(solve(r*[1,1]+[30,60]-[xi,zi],r)) !=Vector{Any}
        ui=-ui
        wi=-wi
        xi=xold
        zi=zold
        br=2
    #Reflexion an der linken unteren Ecke
    elseif xi<=30&& typeof(solve(r*[1,1]+[30,0]-[xold,zold],r)) !=Vector{Any} && typeof(solve(r*[1,1]+[30,0]-[xi,zi],r)) !=Vector{Any}
            ui=-ui
            wi=-wi
            xi=xold
            zi=zold
            br=2
    #Reflexion an der rechten unteren Ecke
    elseif xi>=90 && typeof(solve(-r*[1,1]+[90,0]-[xold,zold],r))!=Vector{Any} && typeof(solve(-r*[1,1]+[90,0]-[xi,zi],r)) !=Vector{Any}
            ui=-ui
            wi=-wi
            xi=xold
            zi=zold
            br=2
    #Reflexion an der rechten oberen Ecke
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

### Übergreifende Funktion zur Berechnung des Monte-Carlo-Modells ###
function monte(xq,zq,n,lambda,sigma)

# Einlesen der Parameter aus der externen NC-Datei
    u =  ncread("Bericht/input_uebung5.nc", "u")
    w =  ncread("Bericht/input_uebung5.nc","w")
    us = sigma*ncread("Bericht/input_uebung5.nc","u2") 
    ws = sigma*ncread("Bericht/input_uebung5.nc","w2")
    #Ersetzen von den Fehlenden Werten mit NaN
    templist  = findall(x->x==-9999.0,u)
    templistw = findall(x->x==-9999.0,w)
    for i in 1: length(templist)
    u[templist[i]]=NaN
    end 
    for i in 1: length(templistw)
        w[templistw[i]]=NaN
    end

#Loop über alle Partikel
    for i in ProgressBar(1:n)
        #Initialisieren aller Parameter für den nächsten Partikel
        xi = xq
        zi = zq
        dt=0
        ui=0
        wi=0
        xold=xq
        zold=zq
        while (ceil(xi+ui*dt) < xgrenz) 
            #Abspeichern der alten Werte
            xolder=xold
            zolder=zold
            xold = xi
            zold = zi
            #Berechnung der Prandtlschicht
            tl, dt = prandltl(zi,xi,lambda,us,ws,k)
            #Berechnung der neuen Positionen
            xi, wi, zi,ui = positionen(xi, wi, zi, tl, ui, dt,xold,zold,xolder,zolder,us,ws,u,w)
            # Gitterauswertung der Konzentration
            exaktgitter(xi, xold, zi, zold,dt,n)
        end

    end
    return  konk #Zurückgeben der Konzentration als Feld
end


# Visualisierung #

## Funktionen zur Visualisierung ##

### Partikelanzahlplot ###
function partikelanzahlplot(xq,zq)
    #Initialisierung der Parameter
    lambda=5
    sigma=1
    # Partikelanzahlen
    na= 100
    nb= 1000
    nc = 5000
    nd= 10000

    # Grafeneinstellungen
    levels=  [0.00001,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2,5]
    fig = Figure(resolution=(1080,1080))
    xs=LinRange(0, xgrenz,xgrenz)
    ys=LinRange(0, zgrenz, zgrenz)

    # Erster Plot 
    a=monte(xq,zq,na,lambda,sigma)
    contourf(fig[1, 1],ys, xs, a,levels=levels)
    contour!(fig[1, 1],ys, xs,a ,levels=levels)

    #Zweiter Plot
    b=monte(xq,zq,nb,lambda,sigma)
    contourf(fig[1, 2],xs, ys, b,levels=levels,title= "N = " *string(nb))
    contour!(fig[1, 2],ys, xs,b ,levels=levels )

    #Dritter Plot
    c=monte(xq,zq,nc,lambda,sigma)
    contourf(fig[2, 1],xs, ys,c,levels=levels,title= "N = " *string(nc))
    contour!(fig[2, 1],ys, xs,c ,levels=levels )

    #Vierter Plot
    d=monte(xq,zq,nd,lambda,sigma)
    contourf(fig[2, 2],xs, ys,d,levels=levels,title= "N = " *string(nd))
    contour!(fig[2, 2],ys, xs,d ,levels=levels )

    #Farbskala
    Colorbar(fig[1:2,3], limits = (0.1, 5), colormap = :viridis, flipaxis = false, label = "kg/s")

    Label(fig[3, :], text = "x(m)")
    Label(fig[:, 0], text = "z(m)")
    Label(fig[0, :], text = "Konzentrationen im Vergleich fuer x = " *string(xq)*"m z = " *string(zq)*"m", textsize = 30)
    save("Bericht/Bilder/3_vergleich_x = "*string(xq)*".png", fig)

    #Löschen der Grafik aus dem Zwischenspeicher
    empty!(fig)
end

### Einzelplot ###
function einzelplot(xq,zq)
    #Initialisierung der Parameter
    n=5000
    lambda=5
    sigma=1

    # Grafeneinstellungen
    levels= [0.00001,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2,5]
    fig = Figure(resolution=(1080,1080))
    xs=LinRange(0, xgrenz,xgrenz)
    ys=LinRange(0, zgrenz, zgrenz)

    #Plot
    a=monte(xq,zq,n,lambda,1)
    contourf(fig[1, 1],ys, xs,a ,levels=levels, xlabel = "x label", ylabel = "y label" )
    contour!(fig[1, 1],ys, xs,a ,levels=levels )

    #Label
    Colorbar(fig[1,2], limits = (0, 5), colormap = :viridis,flipaxis = false,  label = "kg/s")
    Label(fig[0, :], text = "Konzentration bei x = " *string(xq)*"m z = " *string(zq)*"m λ="*string(lambda), textsize = 30)
    Label(fig[2, :], text = "x(m)")
    Label(fig[:, 0], text = "z(m)")
    save("Bericht/Bilder/3_single_x = "*string(xq)*".png", fig)
    #Löschen der Grafik aus dem Zwischenspeicher
    empty!(fig)
end 

### Funktion zur Erstellung des Lambda-Vergleichsplots
function lambdaplot(xq,zq)
    #Initialisierung der Parameter
    n= 1000
    sigma=5
    #Initialisieren der Lambda Werte
    la= 10
    lb= 5
    lc =2
    ld= 1

    #Grafeneinstellungen
    levels= [0.00001,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2]
    fig = Figure(resolution=(1080,1080))
    xs=LinRange(0, xgrenz,xgrenz)
    ys=LinRange(0, zgrenz, zgrenz)

    # Erster Plot
    a=monte(xq,zq,n,la,sigma)
    contourf(fig[1, 1],ys, xs, a,levels=levels)
    contour!(fig[1, 1],ys, xs,a ,levels=levels )

    #Zweiter Plot
    b=monte(xq,zq,n,lb,sigma)
    contourf(fig[1, 2],xs, ys, b,levels=levels)
    contour!(fig[1, 2],ys, xs,b ,levels=levels )

    # Dritter Plot
    c=monte(xq,zq,n,lc,sigma)
    contourf(fig[2, 1],xs, ys,c,levels=levels)
    contour!(fig[2, 1],ys, xs,c ,levels=levels )

    #Vierter Plot
    d=monte(xq,zq,n,ld,sigma)
    contourf(fig[2, 2],xs, ys,d,levels=levels)
    contour!(fig[2, 2],ys, xs,d ,levels=levels )

    #Label
    Colorbar(fig[1:2,3], limits = (0.00001, 5), colormap = :viridis,flipaxis = false, label = "kg/s")
    Label(fig[3, :], text = "x(m)")
    Label(fig[:, 0], text = "z(m)")
    Label(fig[0, :], text = "Konzentrationen im Vergleich x = " *string(xq)*"m z = " *string(zq)*" m bei verschiedenen λ", textsize = 30)
    save( "Bericht/Bilder/3_lambda_x = "*string(xq)*".png", fig)

    #Löschen der Grafik aus dem Zwischenspeicher
    empty!(fig)
end 

### Funktion zur Erstellung des Sigma-Vergleichsplots
function sigmaplot(xq,zq)
    #Initialisierung der Parameter
    l=5
    n= 1000
    siga=1
    sigb=2
    sigc=5
    sigd=10

    #Grafeneinstellungen
    levels=[0.001,0.01,0.025,0.05,0.5,0.75,1.0,1.25,1.5,2,5]
    fig = Figure(resolution=(1080,1080))
    xs=LinRange(0, xgrenz,xgrenz)
    ys=LinRange(0, zgrenz, zgrenz)

    #Erster Plot
    a=monte(xq,zq,n,l,siga)
    contourf(fig[1, 1],ys, xs, a,levels=levels)
    contour!(fig[1, 1],ys, xs,a ,levels=levels )
    
    #Zweiter Plot
    b=monte(xq,zq,n,l,sigb)
    contourf(fig[1, 2],xs, ys, b,levels=levels)
    contour!(fig[1, 2],ys, xs,b ,levels=levels )

    #Dritter Plot
    c=monte(xq,zq,n,l,sigc)
    contourf(fig[2, 1],xs, ys,c,levels=levels)
    contour!(fig[2, 1],ys, xs,c ,levels=levels )

    #Vierter Plot
    d=monte(xq,zq,n,l,sigd)
    contourf(fig[2, 2],xs, ys,d,levels=levels)
    contour!(fig[2, 2],ys, xs,d ,levels=levels )
    
    #Label
    Colorbar(fig[1:2,3], limits = (0.00001, 5), colormap = :viridis,flipaxis = false, label = "kg/s")
    Label(fig[3, :], text = "x(m)")
    Label(fig[:, 0], text = "z(m)")
    Label(fig[0, :], text = "Konzentrationen im Vergleich x = " *string(xq)*"m z = " *string(zq)*" m bei verschiedenen σu & σw", textsize = 30)
    save("Bericht/Bilder/3_sigma_x = "*string(xq)*".png", fig)
    
    #Löschen der Grafik aus dem Zwischenspeicher
    empty!(fig)
end 
    

# Main Funktion
function main()

    ### Aufgbabe a
    xq = 60.5   # m
    zq = 0.5    # m
    einzelplot(xq,zq)
    partikelanzahlplot(xq,zq)
    lambdaplot(xq,zq)
    sigmaplot(xq,zq)

    ### Aufgbabe b
    xq = 15.5  # m
    zq =  65.5 # m
    einzelplot(xq,zq)
    partikelanzahlplot(xq,zq)
    lambdaplot(xq,zq)
    sigmaplot(xq,zq)

end

#Ausführen der Mainfunktion
main()
    