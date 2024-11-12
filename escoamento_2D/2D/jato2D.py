import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

def plot_heatmap(u_k, x, y):
  plt.clf()
  plt.xlabel('x')
  plt.ylabel('y')
  plt.pcolormesh(x,y,u_k, cmap=plt.cm.jet)
  plt.colorbar()

##################################################################
#!		                 Constantes			                     #
##################################################################	

ni     = 151
nj     = 51
nt     = 1000
re     = 900.
dt     = 1.e-3
L      = 5.
H      = 1.
sbr    = 0.99
roh2o  = 1.0  #g/cm3
mih2o  = 0.01 #g/cm s

initial_cond  = 1
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
#!		       Condições iniciais e de contorno	                 #
##################################################################	

u  = np.zeros((nt, nj, ni))
v  = np.zeros((nt, nj, ni))
p  = np.zeros((nt, nj, ni))
d1 = np.zeros((nt, nj, ni))

#condição inicial
u[0,:,:] = initial_cond
v[0,:,:] = 0
p[0,:,:] = initial_cond
#condição de contorno x
u[:,0,:]  = cc_top
u[:,nj-1,:] = cc_bottom
u[:,:,0]  = cc_left
u[:,:,ni-1] = cc_right

#condição de contorno y
v[:,0,:]  = cc_top
v[:,nj-1,:] = cc_top
v[:,:,0]  = cc_top
v[:,:,ni-1] = cc_top

##################################################################
#!		          Quantidade de movimento    	                 #
##################################################################	

for it in range(nt-1):

    for i in range(1,ni-1):
        for j in range(1,nj-1):

            convect_x     = (u[it,j,i+1]*u[it,j,i+1] - u[it,j,i-1]*u[it,j,i-1])/(2*dx) + (u[it,j+1,i]*v[it,j+1,i] - u[it,j-1,i]*v[it,j-1,i])/(2*dy) 
            convect_y     = (v[it,j+1,i]*v[it,j+1,i] - v[it,j-1,i]*v[it,j-1,i])/(2*dy) + (v[it,j,i+1]*u[it,j,i+1] - v[it,j,i-1]*u[it,j,i-1])/(2*dx)
            
            difus_x       = (u[it,j,i+1] - 2*u[it,j,i] + u[it,j,i-1])/dx2 + (u[it,j+1,i] - 2*u[it,j,i] + u[it,j-1,i])/dy2            
            difus_y       = (v[it,j,i+1] - 2*v[it,j,i] + v[it,j,i-1])/dx2 + (v[it,j+1,i] - 2*v[it,j,i] + v[it,j-1,i])/dy2   

            dudx          = (u[it,j,i+1] - u[it,j,i-1])/(2*dx)
            dvdy          = (v[it,j+1,i] - v[it,j-1,i])/(2*dy)         
            
            d1[it,j,i]    = dudx + dvdy

            dpdx          = (p[it,j,i+1] - p[it,j,i-1])/(2*dx) 
            dpdy          = (p[it,j+1,i] - p[it,j-1,i])/(2*dy) 

            u[it+1,j,i]   = u[it,j,i] + dt*(-convect_x - dpdx + (1/re)*difus_x)
            v[it+1,j,i]   = v[it,j,i] + dt*(-convect_y - dpdy + (1/re)*difus_y)

    
    for i in range(1,ni-1):
        for j in range(1,nj-1):
            
            d2d1dx2       = (d1[it,j,i+1] - 2*d1[it,j,i] + d1[it,j,i-1])/dx2
            d2d1dy2       = (d1[it,j+1,i] - 2*d1[it,j,i] + d1[it,j-1,i])/dy2

            dpx           = (p[it,j,i+1] + p[it,j,i-1])/dx2
            dpy           = (p[it,j+1,i] + p[it,j-1,i])/dy2 
            dp            = dpx + dpy

            d2udx2        = (u[it,j,i+1]*u[it,j,i+1] - 2*u[it,j,i]*u[it,j,i] + u[it,j,i-1]*u[it,j,i-1])/dx2
            d2vdy2        = (v[it,j+1,i]*v[it,j+1,i] - 2*v[it,j,i]*v[it,j,i] + v[it,j-1,i]*v[it,j-1,i])/dy2
            d2uvdxy2      = (u[it,j+1,i+1]*v[it,j+1,i+1] - u[it,j-1,i+1]*v[it,j-1,i+1] - u[it,j+1,i-1]*v[it,j+1,i-1] + u[it,j-1,i-1]*v[it,j-1,i-1])/(4*dx*dy)

            conv_esp_pr   = d2udx2 + d2vdy2 + 2*d2uvdxy2

            p[it,j,i]   = ((dx2*dy2)/(2*(dy2 + dx2)))*(dp + conv_esp_pr - (1/re)*(d2d1dx2 + d2d1dy2) - (1/dt)*d1[it,j,i] ) 
            #p[it,j,i]   = sbr*p[it,j,i] + (1. - sbr)*p[it+1,j,i]

    for j in range(0,nj-1):
       p[it,j,0]       = p[it,j,1]
       p[it,j,ni-1]    = 1. 

    for i in range(0,ni-1):
        p[it,0,i]      = .75*p[it,1,i] + .25*p[it,2,i]
        p[it,nj-1,i]   = .75*p[it,nj-2,i] + .25*p[it,nj-3,i]

    if it%10 == 0:
        dmax=0
        jmax=0
        imax=0
        for i in range(ni):
            for j in range(nj):
                if d1[it,j,i] > dmax:
                    dmax = d1[it,j,i]
                    imax = i
                    jmax = j
        print(f'it: {it} -- i: {imax} -- j: {jmax} -- Dilatação: {dmax}')


breakpoint()

plot_heatmap(u[-1,:,:], xi, yi)
plt.show()

