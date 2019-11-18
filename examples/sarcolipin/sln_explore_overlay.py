import matplotlib.pyplot as plt
import nmrglue as ng
import numpy as np

# Specify results files
spectrum = 'hcsesampi4.ft'
log_file = 'sln_explore_log.dat'
wave_file = 'sln_explore_wave.dat'

# Read spectrum with nmrglue
dic,data = ng.pipe.read(spectrum)

# Read files for lineplot
lx, ly = zip(*np.genfromtxt(wave_file,usecols=(2,3),dtype=float))

# Set contour levels
contour_start = 144214*5 # Baseline
contour_num = 30 # number of contour levels
contour_factor = 1.1 # scaling factor between contour levels
levels = contour_start*contour_factor**np.arange(contour_num) 

# Set ppm scales
uc_1h = ng.pipe.make_uc(dic,data,dim=0)
ppm_1h = uc_1h.ppm_scale()
ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
uc_15n = ng.pipe.make_uc(dic,data,dim=1)
ppm_15n = uc_15n.ppm_scale()
ppm_15n_0, ppm_15n_1 = uc_15n.ppm_limits()

# Plot spectrum
fig,ax = plt.subplots(figsize=(3.5,3.0))

s = -1
ax.contour(data,levels=levels,linewidths=0.25,colors='darkblue',
           extent=(ppm_15n_0, ppm_15n_1, ppm_1h_0*s, ppm_1h_1*s))
ax.set_xlim(120,70)
ax.set_ylim(0,5)
ax.set_xlabel('$^{15}$N CS (ppm)',size=10)
ax.set_ylabel('$^{15}$N-$^{1}$H DC (kHz)',size=10)

# Plot line
ax.plot(lx, ly, c='darkred',alpha=1,linewidth=0.5)

# Plot points
xc,yc = zip(*np.genfromtxt(log_file,usecols=(3,4),dtype=float))
names = np.genfromtxt(log_file,usecols=(0),dtype=str)
ax.scatter(xc,yc,c='darkred',alpha=1)

# Text labels
texts = []
for x,y,s in zip(xc,yc,names):
    texts.append(ax.text(x,y,s,size=6,color='darkred'))

plt.tight_layout()
plt.savefig('sln_explore_overlay.jpg',dpi=300)
plt.show()

