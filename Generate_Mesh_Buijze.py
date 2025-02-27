import numpy as np
import matplotlib.pyplot as plt
from Generate_Mesh_Obj_Cubes import save_cubes_to_obj, save_cubes_to_vtk, save_plane_to_vtk
import h5py


def generate_triangle_cylinder_vertices(cube_size, p1, p2, p3, thickness):
    """
    Generate cube elements that approximate a trapezoidal cylinder and return their vertices.
    
    Parameters:
    - bottom_base: Width of the bottom base of the trapezoid.
    - top_base: Width of the top base of the trapezoid.
    - depth: Depth of the trapezoid (in x-direction).
    - height: Height of the cylinder.
    - cube_size: Size of each cube element.
    
    Returns:
    - List of vertices for each cube element.
    """
    # Define the bounding box of the trapezoidal cylinder
    x_min, x_max = min(p1[0],p2[0],p3[0]), max(p1[0],p2[0],p3[0])
    z_min, z_max = min(p1[1],p2[1],p3[1]), max(p1[1],p2[1],p3[1])

    # Generate grid points within the bounding box
    x_coords = np.arange(x_min, x_max+cube_size, cube_size)
    z_coords = np.arange(z_min, z_max+cube_size, cube_size)

    cubes_vertices = []

    for x in x_coords:
        for z in z_coords:
            # Check if the cube center is inside the trapezoidal cross-section
            if is_point_in_triangle( (x,z), p1, p2, p3):
                # Get the 8 vertices of the cube
                vertices = get_cube_vertices(x, z, cube_size/2, cube_size/2, thickness)
                cubes_vertices.append(vertices)
    
    cubes_vertices = np.array(cubes_vertices)
    return cubes_vertices


def is_point_in_triangle(pt, v1, v2, v3):
    """
    Checks if a point is inside a triangle using the barycentric coordinate method.
    Written by ChatGPT4.o, created by Qian 2025/1/7

    :param pt: The point to check (x, y).
    :param v1: First vertex of the triangle (x, y).
    :param v2: Second vertex of the triangle (x, y).
    :param v3: Third vertex of the triangle (x, y).
    :return: True if the point is inside the triangle, False otherwise.
    """
    EPSILON = 1e-9  # Tolerance for floating-point comparisons

    x, y = pt
    x1, y1 = v1
    x2, y2 = v2
    x3, y3 = v3

    # Compute the area of the entire triangle
    denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))

    if abs(denominator) < EPSILON:
        # Degenerate triangle (area is zero)
        return False

    # Compute barycentric coordinates
    a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denominator
    b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denominator
    c = 1 - a - b

    # Check if point is inside the triangle using a tolerance
    return -EPSILON <= a <= 1 + EPSILON and -EPSILON <= b <= 1 + EPSILON and -EPSILON <= c <= 1 + EPSILON


def get_cube_vertices(x, z, cube_len_x, cube_len_z, thickness):
    """
    Get the 8 vertices of a cube centered at (x, y, z) with the given size.
    """
    vertices = [
        [x - cube_len_x/2, 0.0,       z - cube_len_z/2],
        [x + cube_len_x/2, 0.0,       z - cube_len_z/2],
        [x - cube_len_x/2, thickness, z - cube_len_z/2],
        [x + cube_len_x/2, thickness, z - cube_len_z/2],
        [x - cube_len_x/2, 0.0,       z + cube_len_z/2],
        [x + cube_len_x/2, 0.0,       z + cube_len_z/2],
        [x - cube_len_x/2, thickness, z + cube_len_z/2],
        [x + cube_len_x/2, thickness, z + cube_len_z/2]
    ]
    return np.array(vertices)


if __name__ == "__main__":
    pi = 3.1415
    deg2rad = pi/180
    Phi = 70
    Phi_rad = Phi*deg2rad # 70.0
    l1 = 2000 - 100 * 1/np.tan(Phi_rad)
    l2 = 2000 + 100 * 1/np.tan(Phi_rad)
    Thickness = 2000
    Cube_Size_2D = 4
    offset = 0. # 50.0 
    # f - front, b-back
    
    p4f = [4000 - l1, 2800]
    p5f = [4000 - l2, 3000]
    p6f = [4000 - l1, 3000]

    p1f = [l1 - offset*1/np.tan(Phi_rad), 2800.0 + offset]
    p2f = [l1- offset*1/np.tan(Phi_rad), 3000.0 + offset]
    p3f = [l2- offset*1/np.tan(Phi_rad), 2800.0 + offset]
    
    # cube_big_left = [get_cube_vertices(l1/2- offset*1/np.tan(Phi_rad), 2900+ offset, cube_len_x=l1,cube_len_z=200, thickness=Thickness)]
    # cube_big_right = [get_cube_vertices(4000-l1/2, 2900, cube_len_x=l1,cube_len_z=200, thickness=Thickness)]
    # cubes_left = generate_triangle_cylinder_vertices(Cube_Size_2D, p1f, p2f, p3f, Thickness)
    # cubes_right = generate_triangle_cylinder_vertices(Cube_Size_2D, p4f, p5f, p6f, Thickness)
    # cubes_all = np.vstack( (cube_big_left[0][np.newaxis,:,:], cubes_left, cubes_right, cube_big_right[0][np.newaxis,:,:] ) )

    # Dip = 70 deg, Offset = 0
    cube_big_left =  [get_cube_vertices(1000, 2900, cube_len_x=2000, cube_len_z=200, thickness=Thickness)]
    cube_big_right = [get_cube_vertices(3000, 2900, cube_len_x=2000, cube_len_z=200, thickness=Thickness)]
    cubes_all = np.vstack( (cube_big_left[0][np.newaxis,:,:], cube_big_right[0][np.newaxis,:,:] ) )

    # Print the vertices of the first few cubes
    # from Plot_with_Cubes import Plot_BLOCKS
    # 
    # for i, cube in enumerate(cubes[:1]):
    #     print(f"Cube {i+1} vertices:")
    #     for vertex in cube:
    #         print(vertex)
    #     print("\n l1 is", l1)

    # fig = plt.figure(figsize=(12, 8))
    # ax = fig.add_subplot(111, projection="3d")

    # fig, ax = Plot_BLOCKS(fig, ax, cubes_left)
    # fig, ax = Plot_BLOCKS(fig, ax, cube_big_left)
    # fig, ax = Plot_BLOCKS(fig, ax, cubes_right, "red")
    # fig, ax = Plot_BLOCKS(fig, ax, cube_big_right, "red")
    # ax.set_xlim((1900, 2100))
    # plt.show()
    
    
    

    # >>>>>> save_cubes_to_obj(cubes_all,"Buijze19_50.obj")
    
    save_cubes_to_vtk("Buijze19_50_cubes.vtk", cubes_all )


    FaultCenterX,FaultCenterY,FaultCenterZ = (p2f[0]+p4f[0])/2, Thickness/2, (p2f[1]+p4f[1])/2
    FaultNormal = [np.sin(Phi_rad), 0, np.cos(Phi_rad)]
    
    LengthDip = 800.0
    LengthStrike = Thickness
    save_plane_to_vtk( "Buijze19_50_plane.vtk", 
                      [FaultCenterX,FaultCenterY,FaultCenterZ], FaultNormal, 
                       plane_i_size = LengthDip, plane_j_size = LengthStrike ) # plane_i_size = length along the dip direction!
    
    
    
    

    # Save as HDF5
    with h5py.File("Input_Buijze19_Cubes.h5", 'w') as f:
        f.create_dataset('cubes',  data = cubes_all)


    with h5py.File("Input_Buijze19_Geometry.jld2", "w") as f:
        f.create_dataset('Phi',  data = Phi)
        f.create_dataset('Thickness', data = Thickness)
        f.create_dataset('Offset',  data = offset)
        f.create_dataset('Cube_Size_2D',  data = Cube_Size_2D)
        f.create_dataset('FaultCenterX',  data = FaultCenterX)
        f.create_dataset('FaultCenterY',  data = FaultCenterY)
        f.create_dataset('FaultCenterZ',  data = FaultCenterZ)
        f.create_dataset('LengthDip',  data = LengthDip)
        f.create_dataset('LengthStrike',  data = LengthStrike)

    

    print(">>>> Files All Saved >>>>")
    