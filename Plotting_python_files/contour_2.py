import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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

psi1 = np.loadtxt("resultGS.csv", delimiter=",")
psi2 = np.loadtxt("resultSOR.csv", delimiter=",")
psi3 = np.loadtxt("resultBiCGSTAB.csv", delimiter=",")
psi4 = np.loadtxt("resultBiCGSTAB_ILU(0).csv", delimiter=",")
psi5 = np.loadtxt("resultBiCGSTAB_SIP.csv", delimiter=",")

X, Y = np.meshgrid(x, y, indexing="ij")

fig = plt.figure(figsize=(10, 12))
gs = GridSpec(3, 4, figure=fig)

axes = [
    fig.add_subplot(gs[0, 0:2]),
    fig.add_subplot(gs[0, 2:4]),
    fig.add_subplot(gs[1, 0:2]),
    fig.add_subplot(gs[1, 2:4]),
    fig.add_subplot(gs[2, 1:3])
]

psis = [psi1, psi2, psi3, psi4, psi5]
titles = [
    "(a) Gauss--Seidel",
    "(b) SOR",
    "(c) BiCGSTAB",
    "(d) BiCGSTAB + ILU(0)",
    "(e) BiCGSTAB + SIP"

]

for ax, psi, title in zip(axes, psis, titles):
   
    cs = ax.contour(
        X, Y, psi,
        levels=40,
        colors="#0080FF",      
        linewidths=0.7
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

    ax.set_title(title, pad=8)

    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$", rotation=0, labelpad=16)

    ax.set_xlim(3,11)
    ax.set_ylim(0,4)

    ax.set_aspect("equal")

    ax.tick_params(
        direction='in',
        length=4,
        width=0.8,
        top=False,
        right=False
    )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.spines['left'].set_linewidth(0.8)

plt.tight_layout()
plt.savefig("five_contours.png", dpi=600, bbox_inches="tight")
plt.savefig("five_contours.pdf", dpi=600, bbox_inches="tight")
plt.show()

