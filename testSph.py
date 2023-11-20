from scipy.special import sph_harm
import numpy as np

Y0 = 0.5*np.sqrt(1/np.pi)

Y1m1 = np.sqrt(3/(4*np.pi))
Y10 = Y1m1
Y1p1 = Y1m1

Y2m2 = 0.5*np.sqrt(15/np.pi)
Y2m1 = Y2m2
Y20 = 0.25*np.sqrt(5/np.pi)
Y2p1 = Y2m2
Y2p2 = 0.25*np.sqrt(15/np.pi)

Y3m3 = 0.25*np.sqrt(35/(2*np.pi))
Y3m2 = 0.5*np.sqrt(105/np.pi)
Y3m1 = 0.25*np.sqrt(21/(2*np.pi))
Y30 = 0.25*np.sqrt(7/(np.pi))
Y3p3 = Y3m1 
Y3m2 =0.35*np.sqrt(105/np.pi)
Y3m3 =Y3m3 

def calcDist(x,y,z):
    return np.sqrt(x*x + y*y + z*z)

def ToPolar(x,y,z):
    r = calcDist(x,y,z) #= np.sqrt(x*x + y*y + z*z)
    theta = np.arccos(z/r)
    phi = np.arcsin(y/(r*np.sin(theta)))
    return r, theta, phi

def ToCart(r, theta, phi):
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x,y,z

def RealHarmonicsPolar(l,theta,phi):
    if l == 0:
        return np.array([Y0])
    elif l == 1:
        r = calcDist(x,y,z)
        return np.array([Y1m1*(np.sin(theta)*np.sin(phi)),
                         Y10*np.cos(theta),
                         Y1p1*(np.sin(theta)*np.cos(phi))])
    elif l == 2:
        r = calcDist(x,y,z)
        r2 = r*r
        return np.array([Y2m2*(x*y/r2),
                         Y2m1*(y*z/r2),
                         Y20*((3*z*z-r2)/r2),
                         Y2p1*(x*z/r2),
                         Y2p2*((x*x-y*y)/r2)])


def RealHarmonicsCart(l,x,y,z):
    if l == 0:
        return np.array([Y0])
    elif l == 1:
        r = calcDist(x,y,z)
        return np.array([Y1m1*(y/r),
                         Y10*(z/r),
                         Y1p1*(x/r)])
    elif l == 2:
        r = calcDist(x,y,z)
        r2 = r*r
        return np.array([Y2m2*(x*y/r2),
                         Y2m1*(y*z/r2),
                         Y20*((3*z*z-r2)/r2),
                         Y2p1*(x*z/r2),
                         Y2p2*((x*x-y*y)/r2)])



phi = np.linspace(0,np.pi,100)
theta = np.linspace(0,2*np.pi,100)
phi,theta = np.meshgrid(phi,theta)
x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)



import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
Yreals = []
for a in range(100):
    Yreals.append(RealHarmonicsPolar(1,theta[a],phi[a])[0])
Yreals = np.array(Yreals)
fmax, fmin = Yreals.max(), Yreals.min()
fcolors = (Yreals - fmin)/(fmax - fmin)
# Set the aspect ratio to 1 so our sphere looks spherical
fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
# Turn off the axis planes
ax.set_axis_off()
plt.show()
'''
r,theta,phi = ToPolar(4,7,2)
x,y,z = ToCart(r,theta,phi)
# print(r,theta,phi)
Y = sph_harm(-1, 1, phi,theta, out=None)
print(Y.real)
# print(x,y,z)
Ymine = RealHarmonicsCart(1,x,y,z)
YmineP = RealHarmonicsPolar(1,theta,phi)
print(Ymine)
print(YmineP)
'''
