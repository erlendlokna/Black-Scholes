import numpy as np
import matplotlib.pyplot as plt
newparams = {'figure.figsize': (8.0, 4.0), 'axes.grid': True,
             'lines.markersize': 8, 'lines.linewidth': 2,
             'font.size': 14}
from matplotlib import cm
plt.rcParams.update(newparams)

from mpl_toolkits.mplot3d.axes3d import Axes3D



def plot_solution(x, t, U, txt='Solution / Option price'):
   #3d plot probided by lecturer.
    fig = plt.figure(figsize=(12, 8))
    T, X = np.meshgrid(t,x)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_wireframe(T, X, U)
    ax.plot_surface(T, X, U, cmap=cm.coolwarm)
    ax.view_init(azim=45)              # Rotate the figure
    plt.xlabel('time')
    plt.ylabel('stock price $')
    plt.title(txt)

