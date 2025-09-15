import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
import pickle
matplotlib.rcParams['animation.embed_limit'] = 2**128

def plot_heatmap(u_k, x, y):
  plt.clf()
  plt.xlabel('x')
  plt.ylabel('y')
  plt.pcolormesh(x,y,u_k, cmap=plt.cm.jet)
  plt.colorbar()

def plot_heatmapB(v_k, x, y):
  plt.clf()
  plt.xlabel('x')
  plt.ylabel('y')
  plt.pcolormesh(x,y,v_k, cmap=plt.cm.jet)
  plt.colorbar()

def plot_mesh(x,y):

    xi, yj = np.meshgrid(x, y)

    plt.plot(xi, yj, 'k.')
    plt.xlabel('x')
    plt.ylabel('y')
def animate(k):
  plot_heatmap(u[k, :,:], xi, yi)


def format_path(path):
    """""Formats the path string in order to avoid conflicts."""

    if path[-1]!='/':
        path = path + '/'

    if not os.path.exists(path):
        os.makedirs(path)

    return path

###Constantes###
ni     = 101
nj     = 101
nt     = 30000 # laço velocidade
re     = 100
dt     = 1.e-3
L      = 1.
H      = 1.
N = (ni)*(nj)
npr =  80 # laço Poisson
tolp = 1e-8 #critério de parada poisson
tolv = 1e-6/N # critério de parada velocidade
#sbr    = 0.99
roh2o  = 1.0  #g/cm3
mih2o  = 0.01 #g/cm s
w = 1
initial_cond  = 1
cc_bottom     = 0
cc_top        = 1
cc_left       = 0
cc_right      = 0
###Geração da malha###

dx     = L/ni
dy     = H/nj

dx2    = dx**2
dy2    = dy**2
xi     = np.linspace(0,L, ni)
yi     = np.linspace(0,H, nj)

### Condições Iniciais e de Contorno###

u  = np.zeros(( nt, nj, ni))
v  = np.zeros(( nt, nj, ni))
p  = np.zeros(( nt, nj, ni))
F  = np.zeros(( nt, nj, ni))
G  = np.zeros(( nt, nj, ni))

#Geração do grafico dos vetores de velocidade (V=ui+vj)
uvt = np.zeros((nj,ni))
vvt = np.zeros((nj,ni))

### Condições Iniciais ###

u[0,:,:]=0
u[0,nj-2,2:ni-1]=initial_cond
v[0,:,:]=0
p[0,:,:]=0
F[0,:,:]=0
G[0,:,:]=0
###Condições de contorno###
#Técnicas de reflexão aplicadas as paredes:#
u[:,nj-1,2:(ni-1)] = 2*cc_top - u[:,nj-2,2:(ni-1)]  # superior#
u[:,0,2:(ni-1)] = - u[:,1,2:(ni-1)] # inferior#
v[:,2:(nj-1),0] = - v[:,2:(nj-1),1] #esquerda#
v[:,2:(nj-1),ni-1] = - v[:,2:(nj-1),ni-2] #direita#

F[:,1:nj-1,1] = F[:,1:nj-1,ni-1] = 0
G[:,1,1:ni-1] = G[:,nj-1,1:ni-1] = 0

# dp/dx = 0 nas paredes:#
p[:,1:nj-1,0] = p[:,1:nj-1,1] #esquerda#
p[:,1:nj-1,ni-1] = p[:,1:nj-1,ni-2] #direita#
#dp/dy = 0 nas paredes:#
p[:,0,1:ni-1] = p[:,1,1:ni-1] #inferior#
p[:,nj-1,1:ni-1] = p[:,nj-2,1:ni-1]#superior#

###Quantidade de movimento###

it = 0
Rv = 1
while it<=(nt-2) and Rv>tolv:
    itp=0
    Rp=1
    while itp<=npr and Rp>tolp :
        #Cálculo de F nos pontos Ds da malha#
        for i in np.arange(1,ni-2):
            for j in np.arange(1,nj-1):
                #Cálculo de du*2/dx#
                u1m = (u[it,j,i+1]+u[it,j,i+2])*0.5
                u2m = (u[it,j,i]+u[it,j,i+1])*0.5
                if u1m>=0. :
                    u1 = u[it,j,i+1]
                else :
                    u1 = u[it,j,i+2]
                if u2m>=0 :
                    u2 = u[it,j,i]
                else :
                    u2 = u[it,j,i+1]

                du2dxD = (u1m*u1 - u2m*u2)/(dx)
                #Cálculo de d(uv)/dy#
                v1m = (v[it,j+1,i+1]+v[it,j+1,i])*0.5
                v2m = (v[it,j,i+1]+v[it,j,i])*0.5
                if v1m>=0 :
                    u1 = u[it,j,i+1]
                else :
                    u1 = u[it,j+1,i+1]
                if v2m>=0 :
                    u2 = u[it,j-1,i+1]
                else :
                    u2 = u[it,j,i+1]

                duvdyD = (v1m*u1 - v2m*u2)/(dy)

                convD = du2dxD + duvdyD

                #Cálculo de d*2(u)/dx*2+d*2(u)/dy*2#
                d2udx2D = (u[it,j,i] - 2*u[it,j,i+1] + u[it,j,i+2])/dx2
                d2udy2D = (u[it,j-1,i+1] - 2*u[it,j,i+1] + u[it,j+1,i+1])/dy2

                viscD = d2udx2D + d2udy2D

                F[it,j,i+1] = u[it,j,i+1] + (dt*(-convD+(1/re)*viscD))

        #Cálculo de G nos pontos Es da malha#
        for i in np.arange(1,ni-1):
            for j in np.arange(1,nj-2):
                #Cálculo de dv*2/dy#
                v1m = (v[it,j+1,i] + v[it,j+2,i])*0.5
                v2m = (v[it,j,i] + v[it,j+1,i])*0.5

                if v1m>=0 :
                    v1 = v[it,j+1,i]
                else :
                    v1 = v[it,j+2,i]
                if v2m>=0 :
                    v2 = v[it,j,i]
                else :
                    v2 = v[it,j+1,i]

                dv2dyE = v1m*v1 - v2m*v2

                #Cálculo de d(uv)/dx#
                u1m = (u[it,j,i+1]+u[it,j+1,i+1])*0.5
                u2m = (u[it,j,i]+u[it,j+1,i])*0.5

                if u1m>=0 :
                    v1 = v[it,j+1,i]
                else :
                    v1 = v[it,j+1,i+1]

                if u2m>=0 :
                    v2 = v[it,j+1,i-1]
                else :
                    v2 = v[it,j+1,i]

                duvdxE = u1m*v1 - u2m*v2

                convE = dv2dyE + duvdxE

                #Cálculo de  d*2(v)/dx*2+d*2(v)/dy*2#
                d2vdy2E = (v[it,j+2,i] - 2*v[it,j+1,i] + v[it,j,i])/dy2
                d2vdx2E =  (v[it,j+1,i-1] - 2*v[it,j+1,i] + v[it,j+1,i+1])/dx2

                viscE = d2vdx2E + d2vdy2E

                G[it,j+1,i] =  v[it,j+1,i] + dt*(-convE+(1/re)*viscE)

        #Cálculo de pressão nos pontos internos da malha com relaxação#
        for i in np.arange(1,ni-1):
            for j in np.arange(1,nj-1):
                z = p[it,j,i]
                p[it,j,i]   = (dy2*(p[it,j,i+1]+p[it,j,i-1]) + dx2*(p[it,j+1,i]+p[it,j-1,i]) -
                                (dx2*dy2/dt)*((F[it,j,i+1]-F[it,j,i])/dx+(G[it,j+1,i]-G[it,j,i])/dy))/(2*(dy2+dx2))
                p[it,j,i] = w*p[it,j,i] + (1-w)*z

        #Cálculo das pressões nas paredes esquerda, direita, inferior e superior, respectivamente nos laços, considerando
        # taxas de variação  nulas nessas paredes#
        for j in range(1,nj-1):
            p[it,j,0]       = p[it,j,1]
            p[it,j,ni-1]    = p[it,j,ni-2]


        for i in range(1,ni-1):
            p[it,0,i]      = p[it,1,i]
            p[it,nj-1,i]   = p[it,nj-2,i]

        R=0
        for i in range(1,ni-1):
            for j in range (1,nj-1):
                R = R + (((F[it,j,i+1]-F[it,j,i])/dx+(G[it,j+1,i]-G[it,j,i])/dy)/dt -((p[it,j,i+1]-
                            2*p[it,j,i]+p[it,j,i-1])/dx2+(p[it,j+1,i]-2*p[it,j,i]+p[it,j-1,i])/dy2))**2
        Rp = R**0.5
        if it%100==0 and itp%80==0:
            print(f'Iteração total: {it} -- Iteração Poisson: {itp} -- Resíduo: {Rp}')
        itp = itp+1

    #Cálculo de p*(n+1) na malha interna e fronteiras#
    for i in np.arange(1,ni-1):
        for j in np.arange(1,nj-1):
            p[it+1,j,i]=p[it,j,i]

    for j in range(1,nj-1):
        p[it+1,j,0] = p[it+1,j,1]
        p[it+1,j,ni-1] = p[it+1,j,ni-2]

    for i in range(1,ni-1):
        p[it+1,0,i] = p[it+1,1,i]
        p[it+1,nj-1,i] = p[it+1,nj-2,i]

    #Cálculo de u*(n+1) na malha interna#
    for i in np.arange(1,ni-2):
            for j in np.arange(1,nj-1):
                u[it+1,j,i+1]   = (F[it,j,i+1]) - dt*(p[it+1,j,i+1]-p[it+1,j,i])/dx

    #Cálculo de v*(n+1) na malha interna#
    for i in np.arange(1,ni-1):
            for j in np.arange(1,nj-2):
                v[it+1,j+1,i]   = (G[it,j+1,i]) - dt*(p[it+1,j+1,i]-p[it+1,j,i])/dy

    #Atualização de u e v, nas fronteiras com a técnica de reflexão:#
    for i in np.arange(2,ni-1):
        u[it+1,nj-1,i] = 2*cc_top - u[it+1,nj-2,i]
        #u[it+1,nj-1,i] = 2*cc_top - u[it+1,nj-2,i]
        u[it+1,0,i] = -u[it+1,1,i]

    for j in np.arange(2,nj-1):
        v[it+1,j,0] = - v[it+1,j,1]
        v[it+1,j,ni-1] = - v[it+1,j,ni-2]

    R = 0
    for i in np.arange(1,ni-1):
        for j in np.arange(1,nj-1):
            R = R + (np.abs(u[it+1,j,i]-u[it,j,i])+np.abs(v[it+1,j,i]-v[it,j,i]))
    Rv = R/dt
    cont = 0
    if it%100 == 0:
      print(f'Iteração Total : {it} -- Resíduo velocidade: {Rv}')
    it = it+1

for i in np.arange(0,ni):
    for j in np.arange(0,nj):
            uvt[j,i] = u[it,j,i]
            vvt[j,i] = v[it,j,i]

print(f'Iteração Total: {it} -- Resíduo: {Rv}')
save_path = format_path(f'./results/{ni}x{nj}/')

save_name = [f'{save_path}u_df_pred.pkl', f'{save_path}/v_df_pred.pkl']
u_pred    = u[-1,:,:]
v_pred    = v[-1,:,:]
vels      = [u_pred, v_pred]
for i, name in enumerate(save_name):
    with open(name, 'wb') as fp:
        pickle.dump(vels[i], fp) 

# plot_heatmap(u[it,:,:], xi, yi)
# plt.show()
# plot_heatmap(v[it,:,:], xi, yi)

# X, Y = np.meshgrid(xi,yi)
# plt.quiver(X[::1, ::1], Y[::1, ::1], uvt[::1, ::1], vvt[::1, ::1])
# plt.xlabel('x')
# plt.ylabel('y')
