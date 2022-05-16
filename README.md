# schadstoffausbreitung
## Compiler Befehle
 gfortran -Ofast -I /usr/include/ -lnetcdf -lnetcdff -c test.f90
 gfortran -o test test.o -Ofast -I /usr/include/ -lnetcdff -lnetcdf
## Ncview unter windows
### Benötigte Pakete
- xorg-server
- xinit
- xorg-docs
- xlaunch
- ncview
### Befehl
xlaunch && ncview test_result.nc


## Alternative Software Panoply
Braucht Java9
