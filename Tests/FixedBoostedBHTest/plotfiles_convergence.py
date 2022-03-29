
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

N_0 = 128
data_0 = np.loadtxt("output_0_A11broken.txt")

N_1 = 256
data_1 = np.loadtxt("output_0.txt")

#N_2 = 256
#data_2 = np.loadtxt("output.txt")

def get_data(data,N):
    x = np.zeros(N)
    y = np.zeros(N)
    Z1 = np.zeros((N, N))
    Z2 = np.zeros((N, N))
    Z3 = np.zeros((N, N))
    Z4 = np.zeros((N, N))
    for i in range(N) :
        x[i] = data[i, 0]
        y[i] = data[i*N, 1] 
        for j in range(N) :
            Z1[i,j] = data[i*N + j, 2]
            Z2[i,j] = data[i*N + j, 3]
            Z3[i,j] = data[i*N + j, 4]
    return x,np.absolute(Z1[:,N//2])

x0,z0= get_data(data_0,64)
x1,z1= get_data(data_1,128)
#x2,z2= get_data(data_2,256)

plt.plot(x0,z0,label="N=64")
plt.plot(x1,z1,label="N=128")
#plt.plot(x2,z2,label="N=256")

plt.semilogy()
plt.legend()
plt.savefig("convergence_plot.png",dpi=600)
