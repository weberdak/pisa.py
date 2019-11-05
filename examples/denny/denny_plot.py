import matplotlib.pyplot as plt
import numpy as np

# Read points to plot line from wave file
line_file = 'denny_wave.dat'
line_data = np.genfromtxt(line_file,usecols=(2,3),dtype=float)
line_x, line_y = zip(*line_data)

# Read simulated crosspeaks
sim_file = 'denny_log.dat'
sim_names = np.genfromtxt(sim_file,usecols=(0),dtype=str)
sim_x, sim_y = zip(*np.genfromtxt(sim_file,usecols=(3,4),dtype=float))

# Plot data
fig,ax = plt.subplots(figsize=(3.5,3.5))

# Axis limits
ax.set_xlim(200,100)
ax.set_ylim(-2,11)

# Plot crosspeaks
ax.plot(line_x, line_y,c='black',linewidth=0.5,alpha=1)
ax.scatter(sim_x,sim_y,c='black',s=10,marker='s')

# Axis labels
ax.set_xlabel('$^{15}$N CS (ppm)',size=10)
ax.set_ylabel('$^{15}$N-$^{1}$H DC (kHz)',size=10)

# Tick sizes
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(8) 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(8)

# Text labels
texts = []
for x,y,s in zip(sim_x,sim_y,sim_names):
    texts.append(ax.text(x,y,s,size=8,color='black'))

# Output
plt.tight_layout()
plt.savefig('denny_plot.jpg',dip=300)
plt.show()
