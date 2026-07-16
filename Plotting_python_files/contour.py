import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

import matplotlib as mpl
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 16,
    "axes.labelsize": 20,
    "axes.titlesize": 18,
    "legend.fontsize": 15,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
})

x = np.loadtxt("xcoords_3.csv", delimiter=",")
y = np.loadtxt("ycoords.csv", delimiter=",")
psi_1 = np.loadtxt("resultGS.csv", delimiter=",")
psi_2 = np.loadtxt("resultSOR.csv", delimiter=",")
psi_3 = np.loadtxt("resultBiCGSTAB.csv", delimiter=",")
psi_4 = np.loadtxt("resultBiCGSTAB_ILU(0).csv", delimiter=",")
psi_5 = np.loadtxt("resultBiCGSTAB_SIP.csv", delimiter=",")
X, Y = np.meshgrid(x, y, indexing="ij")


plt.rcParams.update({
    "font.size": 25,
    "axes.labelsize": 40,
    "axes.titlesize": 35,
    "xtick.labelsize": 35,
    "ytick.labelsize": 35,
    "legend.fontsize": 15,
})
plt.figure(figsize=(15,12))
ax = plt.gca()
cs = ax.contour(
        X, Y, psi_5,
        levels=100,
        colors="#0080FF",      
        linewidths=1
    )

# Obstacles
ax.add_patch(plt.Rectangle((5,0),1,1,
                               facecolor='white',
                               edgecolor='black',
                               linewidth=1.0,
                               zorder=10))
ax.add_patch(plt.Rectangle((8,0),1,1,
                               facecolor='white',
                               edgecolor='black',
                               linewidth=1.0,
                               zorder=10))
# Remove top and right borders
ax.tick_params(
        direction='in',
        length=4,
        width=0.8,
        top=False,
        right=False
    )

ax.set_xlabel(r"$x$",labelpad = 10)
ax.set_ylabel(r"$y$", rotation=0, labelpad=30)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(0.8)
ax.spines['left'].set_linewidth(0.8)
ax.set_aspect('equal')
ax.tick_params(axis='x', pad=12)   
ax.tick_params(axis='y', pad=12)   
ax.set_ylim(0,2)
ax.set_xlim(4,10)
plt.tight_layout()
plt.show()


data = np.loadtxt("residualGS.csv", delimiter=",")
data2 = np.loadtxt("residualSOR.csv",delimiter = ",")
data3 =  np.loadtxt("residualBiCGSTAB.csv",delimiter = ",")
data4 =  np.loadtxt("residualBiCGSTAB_SIP.csv",delimiter = ",")
data5 =  np.loadtxt("residualBiCGSTAB_ILU(0).csv",delimiter = ",")
x = data[:,0]   # first column
y = data[:,1]   # second column
X = data2[:,0]
Y = data2[:,1]
X_1 =data3[:,0]
Y_1 = data3[:,1] 
X_2 = data4[:,0]
Y_2 =data4[:,1]
X_3 = data5[:,0]
Y_3 =data5[:,1]

plt.rcParams.update({
    "font.size": 18,
    "axes.labelsize": 22,
    "axes.titlesize": 22,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    "axes.linewidth": 1.0,
})

fig, ax = plt.subplots(figsize=(10, 6))

# Publication colour palette
colors = [
    "#1f77b4",  # blue
    "#d62728",  # red
    "#2ca02c",  # green
    "#9467bd",  # purple
    "#ff7f0e"   # orange
]

ax.plot(x,   y,   lw=1.5, color=colors[0], label="Gauss--Seidel")
ax.plot(X,   Y,   lw=1.5, color=colors[1], label="SOR")
ax.plot(X_1, Y_1, lw=1.5, color=colors[2], label="BiCGSTAB")
ax.plot(X_3, Y_3, lw=1.5, color=colors[3], label="BiCGSTAB + ILU(0)")
ax.plot(X_2, Y_2, lw=1.5, color=colors[4], label="BiCGSTAB + SIP")

ax.set_xlabel(r"$\log(i)$")
ax.set_ylabel(r"$\log(r)$", rotation=0, labelpad=30)

# Publication-style ticks
ax.tick_params(
    direction="in",
    length=5,
    width=1.0
)


# Remove unnecessary spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# No whitespace
ax.margins(x=0)

# Compact legend
ax.legend(
    loc="lower left",
    frameon=False,
    handlelength=2.5
)

fig.tight_layout()

fig.savefig(
    "plot.pdf",
    bbox_inches="tight"
)

fig.savefig(
    "plot.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()