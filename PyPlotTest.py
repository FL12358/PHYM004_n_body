import matplotlib.pylab as plt
import numpy as np


def nBody():
    bodyNo = 9
    data = np.genfromtxt("OUTPUT.txt", delimiter="\t", dtype="float")
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    plt.figure(figsize=(10,10))
    a = np.array([0,0.5,1])
    a = np.arange(0,1,1/bodyNo)
    colour = np.tile(a,len(x))
    colour = colour[:len(x)]
    plt.scatter(x,y,s=2,c=colour)
    dim = 5e8
    #plt.ylim(-dim,dim)
    #plt.xlim(-dim,dim)
    plt.axis('equal')
    plt.savefig("testfile.png",dpi=100)
    plt.show
    
nBody()
