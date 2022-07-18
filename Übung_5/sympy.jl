using SymPy
r = symbols("r")

print(typeof(solve(r*[1,1]+[30,60]-[33,63],r)))
#print((solve(r*[1,1]+[30,60]-[33,63],r)))