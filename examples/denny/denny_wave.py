import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

# Read points to plot line from wave file
line_file = 'denny_wave.dat'
line_data = np.genfromtxt(line_file,usecols=(0,2,3),dtype=float)
line_x, line_y, line_z = zip(*line_data)

# Read simulated crosspeaks
sim_file = 'denny_log.dat'
sim_x, sim_y, sim_z = zip(*np.genfromtxt(sim_file,usecols=(0,3,4),dtype=float))

# Plot data
fig = plt.figure(figsize=(3.25,3.0))
spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig)
ax1 = fig.add_subplot(spec[0])
ax2 = fig.add_subplot(spec[1])

# Axis limits
ax1.set_xlim(-1,18)
ax2.set_xlim(-1,18)
ax1.set_ylim(100,200)
ax2.set_ylim(-2,11)

# Plot chemical shifts
ax1.plot(line_x, line_y,c='black',linewidth=0.5,alpha=1)
ax1.scatter(sim_x,sim_y,c='black',s=10,marker='s')

# Plot dipolar couplings
ax2.plot(line_x, line_z,c='black',linewidth=0.5,alpha=1)
ax2.scatter(sim_x,sim_z,c='black',s=10,marker='s')

# Axis labels
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
plt.savefig('denny_wave.jpg',dip=300)
plt.show()
