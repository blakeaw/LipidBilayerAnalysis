import matplotlib.pyplot as plt
import numpy as np

def plot_step_vectors(vectors, colors=None):
	x = vectors[:,0]
	y=vectors[:,1]
	vx = vectors[:,2]
	vy = vectors[:,3]
	step_vec_plot = plt.figure()
	if colors is not None:
		plt.quiver(x,y,vx,vy,color=colors)
	else:
		plt.quiver(x,y,vx,vy)
	plt.title('Lateral Displacement Vectors')
	return plt.show(step_vec_plot)
