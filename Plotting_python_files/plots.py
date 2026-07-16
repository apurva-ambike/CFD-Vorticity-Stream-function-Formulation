import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 16,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "legend.fontsize": 15,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "axes.linewidth": 0.8,
})


xc =pd.read_csv("xcoords_3.csv").to_numpy()
x = np.arange(1,xc.size+1,1)

plt.scatter(x,xc)
plt.show()

for xi in xc:
    plt.axvline(x=xi, color='blue', linewidth=0.5)
plt.show()

dx = pd.read_csv("deltax_3.csv").to_numpy()
plt.scatter(xc,dx)
plt.show()

yc =pd.read_csv("ycoords.csv").to_numpy()
y = np.arange(1,yc.size+1,1)
plt.scatter(y,yc)
plt.show()

for yi in yc:
    plt.axhline(y=yi, color='blue', linewidth=0.5)
plt.show()

 # Larger figure
fig, ax = plt.subplots(figsize=(10, 5))
for xi in xc:
    ax.axvline(x=xi,color ='#4A90E2',linewidth=0.35)

for yi in yc:
    ax.axhline(y=yi,color ='#4A90E2', linewidth=0.35)
ax.margins(x = 0,y= 0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel(r"$x$", fontsize=40,fontweight = 'bold')
ax.set_ylabel(r"$y$", fontsize=40,rotation =0,labelpad = 20,fontweight = 'bold')
ax.tick_params(axis='both', labelsize=30)
#ax.set_xlim(left=4)
ax.set_ylim(0)
ax.set_aspect('equal')        # Important for computational grids
plt.tight_layout()
from matplotlib.patches import Rectangle

ax.add_patch(Rectangle((5, 0), 1, 1,
                       facecolor='white',
                       edgecolor='black',
                       linewidth=1,
                       zorder=10))

ax.add_patch(Rectangle((8, 0), 1, 1,
                       facecolor='white',
                       edgecolor='black',
                       linewidth=1,
                       zorder=10))
#plt.savefig("gridzoom.png", bbox_inches="tight",dpi =1000)
plt.xticks(fontweight='bold')
plt.yticks(fontweight='bold')
ax.tick_params(
    direction="in",
    length=5,
    width=0.8
)

# Thin frame
for spine in ax.spines.values():
    spine.set_linewidth(0.8)
plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle


fig, ax = plt.subplots(figsize=(10, 5))

# Grid lines
for xi in xc:
    ax.axvline(
        x=xi,
        color="#4A90E2",
        linewidth=0.4,
        alpha=0.8
    )

for yi in yc:
    ax.axhline(
        y=yi,
        color="#4A90E2",
        linewidth=0.4,
        alpha=0.8
    )

# Obstacles
ax.add_patch(Rectangle(
    (5, 0), 1, 1,
    facecolor="white",
    edgecolor="black",
    linewidth=0.8,
    zorder=10
))

ax.add_patch(Rectangle(
    (8, 0), 1, 1,
    facecolor="white",
    edgecolor="black",
    linewidth=0.8,
    zorder=10
))

# Labels
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$", rotation=0, labelpad=18)

# Limits
# ax.set_xlim(3, 11)
ax.set_ylim(bottom=0)

# Preserve mesh geometry
ax.set_aspect("equal")

# Publication ticks
ax.tick_params(
    direction="in",
    length=5,
    width=0.8
)

# Thin frame
for spine in ax.spines.values():
    spine.set_linewidth(0.8)

ax.margins(x=0, y=0)

plt.tight_layout()

plt.savefig(
    "grid_zoom.pdf",
    bbox_inches="tight"
)

plt.savefig(
    "grid_zoom.png",
    dpi=600,
    bbox_inches="tight"
)

plt.show()