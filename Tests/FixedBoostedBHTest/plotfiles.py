
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
N = 256
data = np.loadtxt("output_1_A11broken.txt")

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
        #Z4[i,j] = data[i*N + j, 5]

X,Y = np.meshgrid(x,y)
#plt.contourf(X,Y,Z1)

#plt.contourf(X,Y,Z2)
plt.pcolor(x,y,np.abs(Z1), 
            cmap='Reds',norm=colors.LogNorm(vmin=np.abs(Z1).min()+1e-15, vmax=np.abs(Z1).max()))
plt.colorbar()
ax = plt.gca();
#ax.set_xticks(np.arange(0, N, 4));
#ax.set_yticks(np.arange(0, N, 4));
#plt.xlim((0,N))
#plt.ylim((0,N))
ax.grid()
plt.savefig("output_fig1_{}_broken.png".format(N))

plt.clf()

#plt.contourf(X,Y,Z2)
plt.pcolor(x,y,np.abs(Z2), 
            cmap='Reds',norm=colors.LogNorm(vmin=np.abs(Z2).min()+1e-15, vmax=np.abs(Z2).max()))
plt.colorbar()
ax = plt.gca();
#ax.set_xticks(np.arange(0, N, 4));
#ax.set_yticks(np.arange(0, N, 4));
#plt.xlim((0,N))
#plt.ylim((0,N))
ax.grid(which="both")
plt.savefig("output_fig2.png")

plt.figure()
#plt.contourf(X,Y,Z2)
plt.imshow(Z3, extent=[0, N, 0, N], origin='lower',
           cmap='Reds', interpolation='nearest', vmin=-1, vmax=1)
plt.colorbar()
ax = plt.gca();
ax.set_xticks(np.arange(0, N, 4));
ax.set_yticks(np.arange(0, N, 4));
plt.xlim((0,N))
plt.ylim((0,N))
ax.grid()
plt.savefig("output_fig3.png")

#plt.figure()
#plt.contourf(X,Y,Z2)
#plt.imshow(Z4, extent=[0, N, 0, N], origin='lower',
#           cmap='RdGy', interpolation='nearest')#, vmin=-1e-19, vmax=1e-19)
#plt.colorbar()
#ax = plt.gca();
#ax.set_xticks(np.arange(0, N, 4));
#ax.set_yticks(np.arange(0, N, 4));
#plt.xlim((0,N))
#plt.ylim((0,N))
#ax.grid()
#plt.savefig("output_fig4.png")
