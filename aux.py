#===============
# Disk selection
#===============

import matplotlib.pylab as plt
import numpy as np

def select_disks_interactively(width, height, disk_geoms, figsize=(5,5)):
	fig, ax = plt.subplots(figsize=figsize)
	ax.set(xlim=(0, width), ylim = (0, height))

	def _add_disk(event):
	    x, y = event.xdata, event.ydata
	    disk_geom = plt.Circle((x, y), radius=1, alpha=0.3, edgecolor='b')
	    ax.add_artist(disk_geom)
	    disk_geoms.append(disk_geom)
	    fig.canvas.draw()

	fig.canvas.mpl_connect('button_press_event', _add_disk)
	fig.show()


def select_disks_randomly(width, height, n_points, disk_geoms, figsize=(5,5)):
	xs = np.random.uniform(low=0.0, high=width, size=n_points)
	ys = np.random.uniform(low=0.0, high=height, size=n_points)
	centroids = [(x,y) for (x,y) in zip(xs,ys)]

	fig, ax = plt.subplots(figsize=figsize)
	ax.set(xlim=(0, width), ylim = (0, height))

	for c in centroids:
		geom = plt.Circle(c, radius=1, alpha=0.3, edgecolor='b')
		disk_geoms.append(geom)
		ax.add_artist(geom)
	
	fig.show()    
