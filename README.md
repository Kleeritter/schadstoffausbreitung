# Schadstoffausbreitung
Errinerungen an das Programierpraktikum Schadstoffausbreitung
## Installation von Netcdf
Wichtig ist es die Datei netcdf.mod zu haben. Die gibt es bei Ubuntu in dem Paket libnetcdff-dev und ist dann in /usr/include. Dort muss mann die Verlinkung hin machen!
## Compiler Befehle
 - gfortran -Ofast -I /usr/include/ -lnetcdf -lnetcdff -c test.f90
 - gfortran -o test test.o -Ofast -I /usr/include/ -lnetcdff -lnetcdf
 - /. test.exe
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
