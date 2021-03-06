# !/usr/bin/env python3
# -*- coding: utf-8 -*-
def grass():
    import matplotlib
    import matplotlib.pyplot as plt
    #matplotlib.use('GTK4Agg')
    import numpy as np
    from numpy import ma
    from matplotlib import colors, ticker, cm
    from netCDF4 import Dataset
    import math

    def nc_read_from_file_2d_all(filename, varname):

        import numpy as np
        import sys

        try:
            f = open(filename)
            f.close()
        #      print("Load: " + filename + ".")
        except FileNotFoundError:
            print("Error: " + filename + ". No such file. Aborting...")
            sys.exit(1)

        nc_file = Dataset(filename, "r", format="NETCDF4")
        tmp_array = np.array(nc_file.variables[varname][:, :], dtype=type(nc_file.variables[varname]))

        return tmp_array

    def nc_read_from_file_1d_all(filename, varname):

        import numpy as np
        import sys

        try:
            f = open(filename)
            f.close()
        #      print("Load: " + filename + ".")
        except FileNotFoundError:
            print("Error: " + filename + ". No such file. Aborting...")
            sys.exit(1)

        nc_file = Dataset(filename, "r", format="NETCDF4")
        tmp_array = np.array(nc_file.variables[varname][:], dtype=type(nc_file.variables[varname]))

        return tmp_array

    filename = "alla.nc"
    fileout = "uebung4.png"
    units = "s/m$^2$"

    conc = nc_read_from_file_2d_all(filename, "c")

    conc = np.where(conc == -9999.0, np.nan, conc)
    x = nc_read_from_file_1d_all(filename, "x")
    z = nc_read_from_file_1d_all(filename, "z")

    for i in range(0, len(x)):
        #print("ja")
        if (x[i] >= 100.0):
            print(x[i])
            print("nein")
            pg_mod = conc[:,i]
            break

    c0 = 4.63E-02
    gamma = 0.68
    my = 1.3
    zs = 3.4
    pg = np.zeros(len(z))
    for k in range(0, len(z)):
        pg[k] = c0 * math.exp(-gamma * (z[k] / zs) ** my)

    print("plotting....")
    fig = plt.figure()
    ax = plt.axes()

    ax.plot(pg, z, label='grass')
    ax.plot(pg_mod, z, label='monte')

    plt.xlabel('Concentration (' + units + ')')
    plt.ylabel('z (m)')
    plt.legend()
    plt.ylim(top=25)
    plt.show()
    return
