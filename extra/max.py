import numpy as np
import matplotlib.pyplot as plt
# using loadtxt()
lattidude = np.loadtxt("1.csv",
                 delimiter=";", dtype=float, usecols=5 )

temp = np.loadtxt("1.csv",
                 delimiter=";", dtype=float, usecols=0 )


press = np.loadtxt("1.csv",
                 delimiter=";", dtype=float, usecols=1 )                 
#print(max(arr))
lax=[]

print(lattidude)
for i in range(len(lattidude)):
    if lattidude[i] < 50 or temp[i]<0:
        lax.append(i)
        print("ja")
    else:
        lattidude[i]=lattidude[i]/10000
        press[i]=press[i]/100

lattidude=np.delete(lattidude,lax)
temp=np.delete(temp,lax)
press=np.delete(press,lax)
print(lattidude)

plt.plot(lattidude,temp)

plt.show()
plt.close()
plt.plot(lattidude,press)
plt.show()