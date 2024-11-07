import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

def plot_heatmap(u_k, x, y):
  plt.clf()
  plt.xlabel('x')
  plt.ylabel('y')
  plt.pcolormesh(x,y,u_k, cmap=plt.cm.jet, vmin=0, vmax = 150)
  plt.colorbar()

##################################################################
#!		                 Constantes			                     #
##################################################################	

ni     = 151
nj     = 51
nt     = 100
re     = 900.
dt     = 1.e-5
dt2    = dt/2
L      = 5.
H      = 1.
sbr    = 0.99
roh2o  = 1.0  #g/cm3
mih2o  = 0.01 #g/cm s
sbr    = 0.99

initial_cond  = 0
cc_bottom     = 0
cc_top        = 0
cc_left       = 1
cc_right      = 0

print(f'Reynolds: {re}')

##################################################################
#!		                 Geração da malha	                     #
##################################################################	

dx     = L/ni
dy     = H/nj

dx2    = dx**2
dy2    = dy**2

xi     = np.linspace(0,L, ni)
yi     = np.linspace(0,H, nj)

##################################################################
#!		      Condições iniciais e de contorno	                 #
##################################################################	

u  = np.zeros((nt, ni, nj))
v  = np.zeros((nt, ni, nj))
p  = np.zeros((nt, ni, nj))
d1 = np.zeros((nt, ni, nj))

# convect_x = np.zeros((nt, ni, nj))
# convect_y = np.zeros((nt, ni, nj))
# difus_x   = np.zeros((nt, ni, nj))
# difus_y   = np.zeros((nt, ni, nj))
# dpdx      = np.zeros((nt, ni, nj))
# dpdy      = np.zeros((nt, ni, nj))

#condição inicial
u[0,:,:] = initial_cond
v[0,:,:] = initial_cond
#condição de contorno x
u[:,0,:]  = cc_top
u[:,ni-1,:] = cc_bottom
u[:,:,0]  = cc_left
u[:,:,nj-1] = cc_right

#condição de contorno y
v[:,0,:]  = cc_top
v[:,ni-1,:] = cc_top
v[:,:,0]  = cc_top
v[:,:,nj-1] = cc_top

##################################################################
#!		          Quantidade de movimento    	                 #
##################################################################	

for it in tqdm(range(nt-1)):

    for i in range(1,ni-1):
        for j in range(1,nj-1):

            convect_x     = (u[it,i+1,j]*u[it,i+1,j] - u[it,i-1,j]*u[it,i-1,j])/(2*dx) + \
                            (u[it,i,j+1]*v[it,i,j+1] - u[it,i,j-1]*v[it,i,j-1])/(2*dy)
            
            convect_y     = (v[it,i,j+1]*v[it,i,j+1] - v[it,i,j-1]*v[it,i,j-1])/(2*dy) + \
                            (v[it,i+1,j]*u[it,i+1,j] - v[it,i-1,j]*u[it,i-1,j])/(2*dx)
            
            difus_x       = (u[it,i+1,j] - 2*u[it,i,j] + u[it,i-1,j])/dx2 + \
                            (u[it,i,j+1] - 2*u[it,i,j] + u[it,i,j-1])/dy2
            
            difus_y       = (v[it,i+1,j] - 2*v[it,i,j] + v[it,i-1,j])/dx2 + \
                            (v[it,i,j+1] - 2*v[it,i,j] + v[it,i,j-1])/dy2   

            dudx          = (u[it,i+1,j] - u[it,i-1,j])/(2*dx)
            dvdy          = (v[it,i,j+1] - v[it,i,j-1])/(2*dy)         
            
            d1[it,i,j]    = dudx + dvdy

            dpdx          = (p[it,i+1,j] - p[it,i-1,j])/(2*dx) 
            dpdy          = (p[it,i,j+1] - p[it,i,j-1])/(2*dy) 

            u[it+1,i,j]   = u[it,i,j] + dt*(-convect_x - dpdx + (1/re)*difus_x)
            v[it+1,i,j]   = v[it,i,j] + dt*(-convect_y - dpdy + (1/re)*difus_y)


    for i in range(1,ni-1):
        for j in range(1,nj-1):
            
            d2d1dx2       = (d1[it,i+1,j] - 2*d1[it,i,j] + d1[it,i-1,j])/dx2
            d2d1dy2       = (d1[it,i,j+1] - 2*d1[it,i,j] + d1[it,i,j-1])/dy2

            dpx           = p[it,i+1,j] + p[it,i-i,j]
            dpy           = p[it,i,j+1] + p[it,i,j-1]

            d2udx2        = (u[it,i+1,j]*u[it,i+1,j] - 2*u[it,i,j]*u[it,i,j] + u[it,i-1,j]*u[it,i-1,j])/dx2
            d2vdy2        = (v[it,i,j+1]*v[it,i,j+1] - 2*v[it,i,j]*v[it,i,j] + v[it,i,j-1]*v[it,i,j-1])/dy2
            d2uvdxy2      = (u[it,i+1,j+1]*v[it,i+1,j+1] - u[it,i+1,j-1]*v[it,i+1,j-1] - u[it,i-1,j+1]*v[it,i-1,j+1] + u[it,i-1,j-1]*v[it,i-1,j-1])/(4*dx*dy)

            conv_esp_pr   = d2udx2 + d2vdy2 + d2uvdxy2

            p[it,i,j]   = ((dx2*dy2)/(2*(dy2 + dx2)))*(dpx/dx2 + dpy/dy2 - (1/re)*(d2d1dx2 + d2d1dy2) + (1/dt)*d1[it,i,j] + conv_esp_pr) 
            #p[it,i,j]   = sbr*p[it,i,j] + (1.- sbr)*p[it+1,i,j]
    
    for j in range(0,ni):
       p[it,j,0]       = p[it,j,2]
       p[it,j,nj-1]    = 1. 

    for i in range(0,nj):
        p[it,0,i]      = .75*p[it,1,i] + .25*p[it,2,i]
        p[it,ni-1,i]   = .75*p[it,ni-2,i] + .25*p[it,ni-3,i]

breakpoint()
plot_heatmap(u[-1,:,:], yi, xi)
plt.show()