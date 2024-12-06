import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
matplotlib.rcParams['animation.embed_limit'] = 2**128

def plot_heatmap(u_k, x, y):
  plt.clf()
  plt.xlabel('x')
  plt.ylabel('y')
  plt.pcolormesh(x,y,u_k, cmap=plt.cm.jet)
  plt.colorbar()

def plot_mesh(x,y):
    
    xi, yj = np.meshgrid(x, y)

    plt.plot(xi, yj, 'k.')
    plt.xlabel('x')
    plt.ylabel('y')

def animate(k):
  plot_heatmap(u[k, :,:], xi, yi)

##################################################################
#!		                 Constantes			                     #
##################################################################	

ni     = 151
nj     = 51
nt     = 10000
re     = 900.
dt     = 1.e-3
L      = 5.
H      = 1.
sbr    = 0.95
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

u  = np.zeros((nj, ni))
v  = np.zeros((nj, ni))
p  = np.zeros((nj, ni))
d1 = np.zeros((nj, ni))

#condição inicial
u[:,:] = initial_cond
v[:,:] = 0
p[:,:] = initial_cond

#condição de contorno x
u[0,:]  = cc_top
u[nj-1,:] = cc_bottom
u[:,0]  = cc_left
u[:,ni-1] = cc_right

#condição de contorno y
v[0,:]  = cc_top
v[nj-1,:] = v[nj-2,:]
v[:,0]  = cc_top
v[:,ni-1] = v[:,ni-2]

##################################################################
#!		          Quantidade de movimento    	                 #
##################################################################	

for it in range(nt-1):

    un = u.copy()
    vn = v.copy()
    pn = p.copy()
    
    for i in range(1,ni-1):
        for j in range(1,nj-1):
        
            convect_x     = (un[j,i+1]*un[j,i+1] - un[j,i-1]*un[j,i-1])/(2*dx) + (un[j+1,i]*vn[j+1,i] - un[j-1,i]*vn[j-1,i])/(2*dy) 
            convect_y     = (vn[j+1,i]*vn[j+1,i] - vn[j-1,i]*vn[j-1,i])/(2*dy) + (vn[j,i+1]*un[j,i+1] - vn[j,i-1]*un[j,i-1])/(2*dx)
            
            difus_x       = (un[j,i+1] - 2*un[j,i] + un[j,i-1])/dx2 + (un[j+1,i] - 2*un[j,i] + un[j-1,i])/dy2            
            difus_y       = (vn[j,i+1] - 2*vn[j,i] + vn[j,i-1])/dx2 + (vn[j+1,i] - 2*vn[j,i] + vn[j-1,i])/dy2   

            dudx          = (un[j,i+1] - un[j,i-1])/(2*dx)
            dvdy          = (vn[j+1,i] - vn[j-1,i])/(2*dy)         
            
            d1[j,i]       = dudx + dvdy

            dpdx          = (pn[j,i+1] - pn[j,i-1])/(2*dx) 
            dpdy          = (pn[j+1,i] - pn[j-1,i])/(2*dy) 

            u[j,i]   = un[j,i] + dt*(-convect_x - dpdx + (1/re)*difus_x)
            v[j,i]   = vn[j,i] + dt*(-convect_y - dpdy + (1/re)*difus_y)

    
    for i in range(1,ni-1):
        for j in range(1,nj-1):
            
            d2d1dx2       = (d1[j,i+1] - 2*d1[j,i] + d1[j,i-1])/dx2
            d2d1dy2       = (d1[j+1,i] - 2*d1[j,i] + d1[j-1,i])/dy2

            dpx           = (pn[j,i+1] + pn[j,i-1])/dx2
            dpy           = (pn[j+1,i] + pn[j-1,i])/dy2 
            dp            = dpx + dpy

            d2udx2        = (un[j,i+1]*un[j,i+1] - 2*un[j,i]*un[j,i] + un[j,i-1]*un[j,i-1])/dx2
            d2vdy2        = (vn[j+1,i]*vn[j+1,i] - 2*vn[j,i]*vn[j,i] + vn[j-1,i]*vn[j-1,i])/dy2
            d2uvdxy2      = (un[j+1,i+1]*vn[j+1,i+1] - un[j-1,i+1]*vn[j-1,i+1] - un[j+1,i-1]*vn[j+1,i-1] + un[j-1,i-1]*vn[j-1,i-1])/(4*dx*dy)

            conv_esp_pr   = d2udx2 + d2vdy2 + 2*d2uvdxy2

            coef          = (dx2*dy2)/(2*(dy2 + dx2))

            p[j,i]        = coef*(dp + conv_esp_pr - (1/re)*(d2d1dx2 + d2d1dy2) - (1/dt)*d1[j,i]) 
            p[j,i]        = sbr*pn[j,i] + (1. - sbr)*p[j,i]

    for j in range(0,nj-1):
        p[j,0]            = 1#p[j,1]
        p[j,-1]           = p[j,-2] 
    
    for i in range(0,ni-1):
        p[0,i]            = p[1,i]#.75*p[1,i] + .25*p[2,i]
        p[-1,i]           = p[-2,i]#.75*p[nj-2,i] + .25*p[nj-3,i]
    
    
    if it%10 == 0:
        dmax=0
        jmax=0
        imax=0
        for i in range(ni):
            for j in range(nj):
                if d1[j,i] > dmax:
                    dmax = d1[j,i]
                    imax = i
                    jmax = j
        print(f'it: {it} -- i: {imax} -- j: {jmax} -- Dilatação: {dmax}')

#anim = animation.FuncAnimation(plt.figure(), animate, interval = 1, frames = nt, repeat=False)
#anim.save(filename='./escoamento_2D/2D/flow.html', writer="html")
plot_heatmap(u[:,:], xi, yi)
plt.savefig('./escoamento_2D/2D/final_it_result.png', bbox_inches='tight')

