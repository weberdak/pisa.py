import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

# Read in points for line
line_file_1 = 'sln_explore_wave.dat'
line_data_1 = np.genfromtxt(line_file_1,usecols=(0,2,3),dtype=float)
line_x_1, line_y_1, line_z_1 = zip(*line_data_1)

# Read points of PISA simulations
sim_file_1 = 'sln_explore_log.dat'
sim_x_1 = np.genfromtxt(sim_file_1,usecols=(1),dtype=int)
sim_y_1, sim_z_1 = zip(*np.genfromtxt(sim_file_1,usecols=(5,6),dtype=float))
rho_start = 6

# Plot data
fig = plt.figure(figsize=(2.75,2.75))
spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0])
ax2 = fig.add_subplot(spec[1])

ax1.set_xlim(0,32)
ax2.set_xlim(0,32)
ax1.set_ylim(75,140)
ax2.set_ylim(0,5)

ax1.plot(line_x_1, line_y_1,c='darkred',linewidth=0.5,alpha=1)
ax2.plot(line_x_1, line_z_1,c='darkred',linewidth=0.5,alpha=1)

ax1.scatter(sim_x_1+rho_start,sim_y_1,c='darkblue',s=10)
ax2.scatter(sim_x_1+rho_start,sim_z_1,c='darkblue',s=10)

#ax1.set_xlabel('PISA Index',size=10)
ax2.set_xlabel('PISA Index',size=10)
ax1.set_ylabel('$^{15}$N CS (ppm)',size=10)
ax2.set_ylabel('$^{15}$N-$^{1}$H DC (kHz)',size=10)

for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(8)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
for tick in ax2.yaxis.get_major_ticks():
    tick.label.set_fontsize(8)

plt.tight_layout()
plt.savefig('sln_wave.jpg',dip=300)
plt.show()
