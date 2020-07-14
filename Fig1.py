import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

reader = csv.reader(open("/home/buri/PhD work/senescence project/senescenceModel/3DplotData.csv", "r"), delimiter=",")
x = list(reader)
result = np.array(x).astype("float")
x = np.arange(20,82,2)
y = np.arange(101)
X,Y = np.meshgrid(x,y)
Z = result

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, Z, color='black', lw=.5)

# make the panes transparent
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# make the grid lines transparent
ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

ax.set_zlabel("T cells (%)")
ax.set_ylabel("Number of Cell Divisions, X")
ax.set_xlabel("Age, t (years)")
plt.show()
