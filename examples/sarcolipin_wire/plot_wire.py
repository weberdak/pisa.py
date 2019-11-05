import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d import Axes3D


def xyz2contour(infile):

    # Read X,Y,Z data from file
    x = np.genfromtxt(infile,usecols=(0),dtype=float)
    y = np.genfromtxt(infile,usecols=(1),dtype=float)
    z = np.genfromtxt(infile,usecols=(3),dtype=float)
    
    # Normalize data to minimum score
    minScore = min(z)
    z = [ 1/(x/minScore) for x in z ]
    
    # Convert to grid
    xi = np.unique(x)
    yi = np.unique(y)
    X, Y = np.meshgrid(xi, yi)
    Z = griddata(x,y,z,xi,yi, interp='linear')
    return (X,Y,Z)
    

def readXYZ(infile):
    ''' Read 2D spectrum text file in X,Y,Z format (i.e., from NMRPipe) 
    
    Parameters
    ----------
    infile: text file
        Delimited text file of x, y and intensity.
    
    Returns
    -------
    (x, y, z): list of lists of floats.
    '''
    x = np.genfromtxt(infile,usecols=(0),dtype=float)
    y = np.genfromtxt(infile,usecols=(1),dtype=float)
    z = np.genfromtxt(infile,usecols=(3),dtype=float)
    return (x,y,z)
    

def xyz2contour(x,y,z):
    '''Convert x,y,z points of 2D spectrum to arrays required by matplotlib contour plot.
    
    Parameters
    ----------
    x: list of floats
        List of x values
    y: list of floats
        List of y values
    z: list of floats
        List of intensity values corresponding to x and y.
    '''
    xi = np.unique(x)
    yi = np.unique(y)
    X, Y = np.meshgrid(xi, yi)
    Z = griddata(x,y,z,xi,yi, interp='linear')
    return (X,Y,Z)
    
def limits(data):
    return min(data),max(data)

# Read data, normalize and convert to grid
x0,y0,z0 = readXYZ('pisa_fit.dat')
min_x0,max_x0 = limits(x0)
min_y0,max_y0 = limits(y0)
min_z0,max_z0 = limits(z0)
#z0 = [ 1/(x/min_z0) for x in z0 ]
z0_inv = [ ((max_z0-min_z0)-(x-min_z0))/(max_z0-min_z0) for x in z0 ]
z0 = [ 1-(((max_z0-min_z0)-(x-min_z0))/(max_z0-min_z0)) for x in z0 ]
X0,Y0,Z0 = xyz2contour(x0,y0,z0)
X0,Y0,Z0_inv = xyz2contour(x0,y0,z0_inv)

contour_start = 0.95 # Baseline
contour_num = 30 # number of contour levels
contour_factor = 1.005 # scaling factor between contour levels
levels = contour_start*contour_factor**np.arange(contour_num) 


fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111, projection='3d')

#fig, ax = plt.subplots(figsize=(2.5,2.5))

ax.plot_wireframe(X0,Y0,Z0,linewidth=0.25,color='darkblue',alpha=0.5)
ax.contour(X0,Y0,Z0_inv,linewidths=0.25,colors='darkblue',levels=levels,extent=(min_x0, max_x0, min_y0, max_y0),zdir='z',offset=0.0001)
ax.w_xaxis.pane.fill = False
ax.w_yaxis.pane.fill = False
ax.w_zaxis.pane.fill = False

ax.set_xlim(0,90)
ax.set_ylim(0,360)
ax.set_zlim(0,1.00)
ax.set_xlabel('Tilt angle ($\\tau$,$^\circ$)',size=10)
ax.set_ylabel('Rotation angle ($\\rho$$_{R6}$,$^\circ$)',size=10)
ax.set_zlabel('Normalized fit score',size=10)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
ax.xaxis.set_ticks(np.arange(0, 91, 10))
ax.yaxis.set_ticks(np.arange(0, 361, 45))
for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(6)

plt.tight_layout()
plt.savefig('plot_wire.jpg',dpi=300)
plt.show()
