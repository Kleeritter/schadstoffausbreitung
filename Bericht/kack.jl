using NetCDF
using Random, Distributions
using ProgressBars
using PlotlyJS

# Definition der globalen Variablen


n= 10^3             # Anzahl Partikel
zq = 0              # Quellort z Komponente in m
xq = 0.5            # Quellort x Komponente in m
xgrenz= 110         # Grenze in x Richtung in m
zgrenz=25           # Grenze in z Richtung in m
dx = 1              # Gitterweite x Richtung in m    
dz = 1              # Gitterweite z Richtung in m
ustern = 0.35       # Schubspannungsgeschwindigkeit in m/s
k = 0.38            # Kappa
znull = 0.008       # Rauhigkeitslaenge
sigu = 2.5 * ustern # Standartabweichung u m/sm/s
sigw = 1.3 * ustern # Standartabweichung w m/s m/s
q= 0.7              # Konzentration fuer das Montecarlo Modell in m/s
#rl= exp(- dt/tl)    # Berechnung von rl
units = "g/m^3"     # Einheit fuer Graphen
ubalken=5
## Arrays Initialisieren


#Array wird in der Größe xgrenz*zgrenz aufgespannt und mit 0 gefüllt
konk=zeros(xgrenz,zgrenz)



# Montecarlo-Modell #

##Funktionen zur Berechnung des Montecarlo Modell ##

### Berechnung der Prandtlschicht###
function prandltl(zi,xi)
    #### Berechnung der mittleren Windgeschwindigkeit ####
    if zi < znull 
        ubalken = 0 #kein mittlerer Wind für zi < znull 
    else
        ubalken = (ustern / k) * log(abs(zi) / znull) # Log. Windprofil
    end
    ##Berechnung des Lagrangeschen Zeitschrittes tl
    tl = ((k * ustern) / sigw ^ 2) * abs(zi)

    #### Berechnung des Zeitschrittes dt
    if (0.1*tl)>((k * ustern) / sigw ^ 2) * abs(2) 
        dt = 0.1*tl                             # Normalfall
    else
        dt = ((k * ustern) / sigw ^ 2) * abs(2) #falls dt kleiner als tl in 2 m Hoehe wird dt auf tl(2m) gesetzt
    end
    return tl,  dt, ui
end  


### Berechnung der Positionen###
function positionen(xi, wi, zi,tl, ui, dt)
    #Lagrangesche Autokorellationsfunktion
    rl = exp(- dt / tl)
    Random.seed!()
    d = Normal()
    rr = rand(d, 1)[1]

    xi = xi + ui * dt
    wi=rl*wi + sqrt((1 - rl^2))*sigw* rr
    zi = zi + wi * dt
return xi, wi, zi

end
### Gitterauswertung ###
function exaktgitter(xi, xold, zi, zold,dt)
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

### Geradengleichung zur exakten Gitterauswertung ### 
function gg(xold, zold, xi, zi, t) 
    xg = xold + t * (xi - xold) 
    zg = zold + t * (zi - zold)
    if xg>xgrenz
        xg=xgrenz
    end
    #Die Geradengleichung kann negative Werte liefern.
    #Um einen Abbruch zu verhindern werden zg und xg in diesem Fall auf 1 gesetzt
    if floor(zg) <1 
        zg=1        
    end
    if floor(xg) <1
        xg=1
    end
    return convert(Int64, floor(xg)), convert(Int64, floor(zg)) # Zurückgeben der Werte als Integer
end

### Berechnung des Montecarlo-Modells ###
function monte()
    for i in ProgressBar(1:n+1) #Loop über alle Partikel
        #Initialisieren der Parameter für jeden Partikel
        xi=xq
        zi=zq
        ui=ubalken
        wi=wbalken
        dt=0
        while (ceil(xi+ui*dt) < xgrenz)
            xold=xi #Zwischenspeichern der alten Werte für x 
            zold=zi #Zwischenspeichern der alten Werte für y
            if zi<znull #Berücksichtigen der Totalreflexion
                zi= zi+ 2*abs(zi-znull) #Reflexion
                wi= -wi
                tl, dt,ui = prandltl(zi,xi) #Berechnung der Prandtlschicht
                xi,wi,zi =positionen(xi, wi, zi,tl, ui, dt) #Berechnung der neuen Positionen
                exaktgitter(xi, xold, zi, zold,dt) #Exakte Gitterauswertung
            else
                tl, dt,ui = prandltl(zi,xi) #Berechnung der Prandtlschicht
                xi,wi,zi =positionen(xi, wi, zi,tl, ui, dt) #Berechnung der neuen Positionen
                exaktgitter(xi, xold, zi, zold,dt) #Exakte Gitterauswertung
            end
        end

    end
    return konk #Zurückgeben der Konzentration als Feld
end

# Prairie-Grass-Experiments #
function prairie_grass(konk)
    pg_mod= [] #Initialisieren
    #Ermitteln der benötigten Werte aus dem Montecarlo-Modells
    for i in 100:xgrenz
        for j in 1:zgrenz
            push!(pg_mod, konk[i,j])
        end   
    end
    c0 = 4.63E-02
    gamma = 0.68
    my = 1.3
    zs = 3.4
    z= collect(1:zgrenz)
    pg = zeros(length(z)+1)
for k in 1:zgrenz
   pg[k] = c0 * exp(-gamma * (z[k]/zs)^my)
end
return pg, pg_mod
end

### Visualisierung ###

 function grafen(pg,pg_mod)
    #Darstellung mit PlotlyJS
    title="Vergleich Montecarlo Modell mit dem Prairie-Grass Experiment"
    xaxis_title="Konzentration (" * units * ")"
    plot_prairie=scatter(y=collect(1:zgrenz),x=pg,name="prairie_grass",showlegend=true ,)
    plot_monte= scatter(y=collect(1:zgrenz),x=pg_mod,name=" Montecarlo",showlegend=true ,)
    savefig(plot([plot_prairie,plot_monte],Layout(title= title, xaxis_title=xaxis_title, yaxis_title="z(m)",xaxis_range=[-0.001, 0.05], yaxis_range=[0, 25])), "Bericht/Bilder/2.png")
 end




 ### Main ####
function main()
    #Zunächst wird das Montecarlo-Modell berechnet
    konzentrationen=monte() #Die Konzentration wird in einem Array zwischengespeichert
    #Als nächstes wird das prairie_grass-Experiment mit der berechneten Konzentration durchgeführt
    pg,pg_mod=prairie_grass(konzentrationen) #Die Arrays werden zwischengespeichert
    # Als letztes kommt es zur Visualisierung 
    grafen(pg,pg_mod)
    
end


main() #Ausführen der Main Funktion