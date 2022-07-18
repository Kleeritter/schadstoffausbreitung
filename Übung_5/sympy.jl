using SymPy
r = symbols("r")

print(linsolve(r*[1,1]+[30,60]-[33,63],r))