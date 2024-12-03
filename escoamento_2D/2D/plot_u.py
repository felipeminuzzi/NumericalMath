import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_heatmap(u_k, x, y):
  plt.clf()
  plt.xlabel('x')
  plt.ylabel('y')
  plt.pcolormesh(x,y,u_k, cmap=plt.cm.jet)
  plt.colorbar()

file_path = './escoamento_2D/2D/u.dat'
df = pd.read_csv(file_path, sep='   ', decimal='.')
df = df.rename(columns={'160':'x', 'Unnamed: 1':'y', 'Unnamed: 2':'u'})

ni = 161
nj = 51

xi = df.iloc[0:161,0].values.astype(float)
yj = df.iloc[::161,1].values.astype(float)
u_values = df['u'].values.reshape(nj, ni).astype(float)

plot_heatmap(u_values, xi, yj)
plt.show()
