import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import matplotlib as mpl


mpl.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsmath}",
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 16,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "legend.fontsize": 15,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
})



x = np.loadtxt("xcoords_3.csv", delimiter=",")
y = np.loadtxt("ycoords.csv", delimiter=",")
X, Y = np.meshgrid(x, y, indexing="ij")
psi = np.loadtxt("psi_final_1.csv", delimiter=",")
omega = np.loadtxt("omega_final_1.csv", delimiter=",")
u = np.loadtxt('u_final_1.csv', delimiter=',')
v = np.loadtxt('v_final_1.csv', delimiter=',')
info = np.loadtxt('VTE_data_1.csv',delimiter = ',')
X_4 = info[:,1]
Y_4 = info[:,7]
Y_5 = info[:,2]
Y_6 = info[:,4]
X_5 = info[:,0]

plt.rcParams.update({
    "font.size": 25,
    "axes.labelsize": 25,
    "axes.titlesize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 15,
})
min_psi = np.min(psi)

#neg_levels  = np.linspace(-0.057, -0.001, 25)
neg_levels = -np.geomspace(0.8*abs(min_psi), 1e-5, 20)
neg_levels = np.sort(neg_levels)
near_zero   = np.linspace(-1e-5, 0.0, 3)
pos_levels  = np.linspace(0.005, 6.0,40)
levels = np.unique(np.concatenate([neg_levels, near_zero, pos_levels]))

fig, ax = plt.subplots(figsize=(15,12))
ax.contour(X, Y, psi, levels=levels[levels < 0], colors="#0080FF",
           linewidths=0.8, linestyles='solid')   # negative — recirculation
ax.contour(X, Y, psi, levels=levels[levels > 0], colors="#0080FF",
           linewidths=0.8)                         # positive — main flow

ax.add_patch(Rectangle((5, 0), 1, 1, color='black', fill=False, linewidth=1.2))
ax.add_patch(Rectangle((8, 0), 1, 1, color='black', fill=False, linewidth=1.2))

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$',rotation =0, labelpad =20)
ax.set_aspect('1.35')
ax.margins(x=0)
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
ax.tick_params(axis='x', pad=12)  
ax.tick_params(axis='y', pad=10)  
plt.tight_layout()
plt.show()



plt.rcParams.update({
    "font.size": 40,
    "axes.labelsize": 45,
    "axes.titlesize": 45,
    "xtick.labelsize": 38,
    "ytick.labelsize": 38,
    "legend.fontsize": 15,
})
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
fig_width = 25       # was 16 -> bigger canvas
fig_height = fig_width * (5/14) * 1.15   # keeps proportions, bigger overall

fig, ax = plt.subplots(figsize=(fig_width, fig_height))


# ---------------- Colour limits ----------------
vmin = np.nanpercentile(omega, 2)
vmax = np.nanpercentile(omega, 98)

norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
levels = np.linspace(vmin, vmax, 100)

# ---------------- Filled contours ----------------
cf = ax.contourf(
    X, Y, omega,
    levels=levels,
    cmap="RdBu_r",
    norm=norm,
    extend="both"
)

cf.set_rasterized(True)

# ---------------- Contour lines ----------------
levels_neg = levels[levels < 0]
levels_pos = levels[levels > 0]

line_levels = np.concatenate([levels_neg[::10], levels_pos[::5]])  # denser sampling on the smaller positive side
line_levels = line_levels[np.abs(line_levels) > 1e-12]

cs = ax.contour(
    X, Y, omega,
    levels=line_levels,
    colors="k",
    linewidths=0.4,
    alpha=0.35,
    linestyles='solid'
)
# ---------------- Axes ----------------
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$",rotation =0,labelpad = 30)
ax.set_ylim(0, 5)
ax.set_xlim(0,14)
ax.set_aspect('equal')
ax.tick_params(axis='x', pad=12)   # Move x-axis labels downward
ax.tick_params(axis='y', pad=12)   # Move y-axis labels leftward

ax.add_patch(Rectangle((5, 0), 1, 1, color='black', fill=False, linewidth=0.8))
ax.add_patch(Rectangle((8, 0), 1, 1, color='black', fill=False, linewidth=0.8))


# ---------------- Colorbar (matched to actual rendered ax height) ----------------
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0)
cbar = fig.colorbar(cf, cax=cax)
ticks = np.linspace(vmin, vmax, 7)   # or however many ticks you want
cbar.set_ticks(ticks)
cbar.set_ticklabels([f"{t:.2f}" for t in ticks])
cbar.ax.tick_params(pad=20)

cbar.set_label(r"$\omega$",rotation = 0,labelpad = 30)
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig("vorticity_contours.pdf", dpi=400, bbox_inches="tight", pad_inches=0.05)
plt.show()

from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(15, 5))
gs = GridSpec(1, 2, figure=fig, wspace=0.35)

# ---------------- First plot ----------------
ax1 = fig.add_subplot(gs[0, 0])

ax1.plot(X_4, Y_4,linewidth = 1.0)

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.margins(x=0)

# ---------------- Second plot ----------------
ax2 = fig.add_subplot(gs[0, 1])

ax2.plot(X_5, Y_5,label = " ",
         linewidth=1.0)

ax2.plot(X_5, Y_6,label = " ",
         linewidth=1.0)

ax2.legend(frameon=False)

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.margins(x=0)
plt.show()


import matplotlib.gridspec as gridspec
min_psi = np.min(psi)

neg_levels = -np.geomspace(0.8*abs(min_psi), 1e-5, 20)
neg_levels = np.sort(neg_levels)
near_zero  = np.linspace(-1e-5, 0.0, 3)
pos_levels = np.linspace(0.005, 6.0, 75)
levels = np.unique(np.concatenate([neg_levels, near_zero, pos_levels]))

fig = plt.figure(figsize=(18, 8))
gs = gridspec.GridSpec(2, 1, figure=fig, hspace=0.25)

# ---- Common plotting function to avoid repeating styling code ----
def plot_zone(ax, xlim, ylim):
    ax.contour(X, Y, psi, levels=levels[levels < 0], colors="#0080FF",
               linewidths=0.8, linestyles='solid')   # negative — recirculation
    ax.contour(X, Y, psi, levels=levels[levels > 0], colors="#0080FF",
               linewidths=0.8)                         # positive — main flow

    ax.add_patch(Rectangle((5, 0), 1, 1, color='black', fill=False, linewidth=1.2))
    ax.add_patch(Rectangle((8, 0), 1, 1, color='black', fill=False, linewidth=1.2))

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$', rotation=0, labelpad=20)
    ax.set_aspect('equal')
    ax.margins(x=0)
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
    ax.tick_params(axis='x', pad= 2)
    ax.tick_params(axis='y', pad=10)

# ---- Zone 1: Upstream eddy + inter-block cavity vortex ----
ax1 = fig.add_subplot(gs[0, 0])
plot_zone(ax1, xlim=(2, 13), ylim=(0, 3))

# ---- Zone 2: Downstream wake bubble ----
ax2 = fig.add_subplot(gs[1, 0])
plot_zone(ax2, xlim=(8, 19), ylim=(0, 3))

plt.tight_layout()
plt.show()
print(vmin)
print(vmax)