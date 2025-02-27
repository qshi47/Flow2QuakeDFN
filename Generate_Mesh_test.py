from Generate_Mesh_Buijze import is_point_in_triangle
import numpy as np
import matplotlib.pyplot as plt



p1 = [2036.4011040061653, 2800]
p2 = [1963.5988959938347, 3000]
p3 = [2036.4011040061653, 3000]

x_min, x_max = min(p1[0],p2[0],p3[0]), max(p1[0],p2[0],p3[0])
z_min, z_max = min(p1[1],p2[1],p3[1]), max(p1[1],p2[1],p3[1])
cube_size = 4.0

# Generate grid points within the bounding box
x_coords = np.arange(x_min, x_max + cube_size, cube_size)
z_coords = np.arange(z_min, z_max + cube_size, cube_size)

cubes_inside = []

for x in x_coords:
    for z in z_coords:
        # Check if the cube center is inside the trapezoidal cross-section
        if is_point_in_triangle( (x,z), p1, p2, p3):
            cubes_inside.append([x,z])

cubes_inside = np.array(cubes_inside)

plt.plot(cubes_inside[:,0], cubes_inside[:,1], '.')
plt.plot(np.vstack([p1,p2,p3,p1])[:,0], np.vstack([p1,p2,p3,p1])[:,1], 'r-')
plt.show()
