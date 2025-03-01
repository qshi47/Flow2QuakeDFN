import numpy as np
import matplotlib.pyplot as plt
import h5py
from tqdm import tqdm

from Functions_Generate_Cubes import is_point_in_triangle, get_cube_vertices, save_cubes_to_vtk, save_plane_to_vtk

def Generate_LeftTriangle_Cylinder_Vertices(cube_size, p1, p2, p3, thickness):
    # Generate grid points within the bounding box
    x_coords = np.arange(p1[0]+cube_size/2, p3[0]+cube_size/2, cube_size)
    z_coords = np.arange(p1[1]+cube_size/2, p2[1]+cube_size/2, cube_size)

    cubes_vertices = []
    with tqdm(total = len(x_coords)*len(z_coords), desc="Total cubes") as pbar:
        for x in x_coords:
            for z in z_coords:
                # Check if point_to_be_examined is inside the trapezoidal cross-section
                point_to_be_examined = (x+cube_size/2, z+cube_size/2)
                if is_point_in_triangle( point_to_be_examined, p1, p2, p3):
                    # Get the 8 vertices of the cube
                    vertices = get_cube_vertices(x, z, cube_size/2, cube_size/2, thickness)
                    cubes_vertices.append(vertices)
                pbar.update(1)
    
    return np.array(cubes_vertices)

def Generate_RightTriangle_Cylinder_Vertices(cube_size, p4, p5, p6, thickness):
    # Generate grid points within the bounding box
    x_coords = np.arange(p6[0]-cube_size/2, p5[0]-cube_size/2, -cube_size)
    z_coords = np.arange(p6[1]-cube_size/2, p4[1]-cube_size/2, -cube_size)

    cubes_vertices = []
    with tqdm(total = len(x_coords)*len(z_coords), desc="Total cubes") as pbar: 
        for x in x_coords:
            for z in z_coords:
                # Check if point_to_be_examined is inside the trapezoidal cross-section
                point_to_be_examined = (x-cube_size/2, z-cube_size/2)
                if is_point_in_triangle( point_to_be_examined, p4, p5, p6):
                    # Get the 8 vertices of the cube
                    vertices = get_cube_vertices(x, z, cube_size/2, cube_size/2, thickness)
                    cubes_vertices.append(vertices)
                pbar.update(1)
    
    cubes_vertices = np.array(cubes_vertices)
    return cubes_vertices



if __name__ == "__main__":
    OutputCubesVerticesFilename = "Input_Buijze19_Cubes.h5"
    OutputBasicGeometryFilename = "Input_Buijze19_Geometry.jld2"
    OutputCubesVTKFilename = "Buijze19_00_cubes.vtk"
    OutputPlaneVTKFilename = "Buijze19_00_plane.vtk"
    
    # Reservoir parameters
    Phi = 70 # degree
    Thickness = 2000 # meter
    Height = 200 # meter
    offset = 0.0 # 0.0 | 50.0 # meter
    Cube_Size_2D = 4 # meter
    
    Phi_rad = Phi*np.pi/180
    BigCubeLength = 2000-0.5*((Height+offset)/np.tan(Phi_rad))
    LeftBigCubeCenterX = BigCubeLength/2
    RightBigCubeCenterX = 4000 - BigCubeLength/2
    LeftBigCubeCenterZ = 2800 + Height/2 + offset
    RightBigCubeCenterZ = 2800 + Height/2
    
    # Six vertices of two triangles
    LeftTriangleP1 = [BigCubeLength, 2800+offset]   
    LeftTriangleP2 = [BigCubeLength, 2800+offset+Height]
    LeftTriangleP3 = [BigCubeLength + Height/np.tan(Phi_rad), 2800+offset]
    
    RightTriangleP4 = [4000-BigCubeLength, 2800]
    RightTriangleP5 = [4000-BigCubeLength - Height/np.tan(Phi_rad), 3000]
    RightTriangleP6 = [4000-BigCubeLength, 3000]
    
    
    # Vertices of all cube 
    cube_big_left =  [get_cube_vertices(LeftBigCubeCenterX, LeftBigCubeCenterZ, BigCubeLength, Height, Thickness)]
    cube_big_right = [get_cube_vertices(RightBigCubeCenterX, RightBigCubeCenterZ, BigCubeLength, Height, Thickness)]
    cubes_left =  Generate_LeftTriangle_Cylinder_Vertices(Cube_Size_2D, LeftTriangleP1, LeftTriangleP2, LeftTriangleP3, Thickness)
    cubes_right = Generate_RightTriangle_Cylinder_Vertices(Cube_Size_2D, RightTriangleP4, RightTriangleP5, RightTriangleP6, Thickness)
    cubes_all = np.vstack( (cube_big_left[0][np.newaxis,:,:], cubes_left, cubes_right, cube_big_right[0][np.newaxis,:,:] ) )
    

    """
    # Two blocks
    # Dip = 70 deg, Offset = 0
    cube_big_left =  [get_cube_vertices(1000, 2900, cube_len_x=2000, cube_len_z=200, thickness=Thickness)]
    cube_big_right = [get_cube_vertices(3000, 2900, cube_len_x=2000, cube_len_z=200, thickness=Thickness)]
    cubes_all = np.vstack( (cube_big_left[0][np.newaxis,:,:], cube_big_right[0][np.newaxis,:,:] ) )
    """

    """
    # One big block
    # Dip = 70 deg, Offset = 0
    cube_big_one =  [get_cube_vertices(2000, 2900, cube_len_x=4000, cube_len_z=200, thickness=Thickness)]
    cubes_all = np.array(cube_big_one)
    """

    # Fault parameters 
    FaultCenterX,FaultCenterY,FaultCenterZ = (LeftTriangleP2[0]+RightTriangleP4[0])/2, Thickness/2, (LeftTriangleP2[1]+RightTriangleP4[1])/2
    FaultNormal = [np.sin(Phi_rad), 0, np.cos(Phi_rad)]
    LengthDip = 800.0
    LengthStrike = Thickness



    # Save as VTK
    print("------ Saving VTK files -------")
    save_cubes_to_vtk(OutputCubesVTKFilename, cubes_all)

    save_plane_to_vtk( OutputPlaneVTKFilename, 
                      [FaultCenterX,FaultCenterY,FaultCenterZ], FaultNormal, 
                       plane_i_size = LengthDip, plane_j_size = LengthStrike ) # plane_i_size = length along the dip direction!
    # save_cubes_to_obj(cubes_all,"Buijze19_50.obj")


    # Save as HDF5
    print("------ Save HDF5 files -------")
    with h5py.File(OutputCubesVerticesFilename, 'w') as f:
        f.create_dataset('cubes',  data = cubes_all)


    with h5py.File(OutputBasicGeometryFilename, "w") as f:
        f.create_dataset('Phi',           data = Phi)
        f.create_dataset('Thickness',     data = Thickness)
        f.create_dataset('Offset',        data = offset)
        f.create_dataset('Cube_Size_2D',  data = Cube_Size_2D)
        f.create_dataset('FaultCenterX',  data = FaultCenterX)
        f.create_dataset('FaultCenterY',  data = FaultCenterY)
        f.create_dataset('FaultCenterZ',  data = FaultCenterZ)
        f.create_dataset('LengthDip',     data = LengthDip)
        f.create_dataset('LengthStrike',  data = LengthStrike)

    
    print(f"Offset {offset} m")
    print(f"Cube numer : {cubes_all.shape[0]}")
    print(">>>> Files All Saved >>>>")


# Notations: f - front, b-back
# l1 = 2000 - 100 * 1/np.tan(Phi_rad)
# l2 = 2000 + 100 * 1/np.tan(Phi_rad)

# p4f = [4000 - l1, 2800]
# p5f = [4000 - l2, 3000]
# p6f = [4000 - l1, 3000]
# p1f = [l1 - offset*1/np.tan(Phi_rad), 2800.0 + offset]
# p2f = [l1- offset*1/np.tan(Phi_rad), 3000.0 + offset]
# p3f = [l2- offset*1/np.tan(Phi_rad), 2800.0 + offset]
    
# def generate_triangle_cylinder_vertices(cube_size, p1, p2, p3, thickness):
# # Define the bounding box of the trapezoidal cylinder
# x_min, x_max = min(p1[0],p2[0],p3[0]), max(p1[0],p2[0],p3[0])
# z_min, z_max = min(p1[1],p2[1],p3[1]), max(p1[1],p2[1],p3[1])

# # Generate grid points within the bounding box
# x_coords = np.arange(x_min, x_max+cube_size, cube_size)
# z_coords = np.arange(z_min, z_max+cube_size, cube_size)

# cubes_vertices = []

# for x in x_coords:
#     for z in z_coords:
#         # Check if the cube center is inside the trapezoidal cross-section
#         if is_point_in_triangle( (x,z), p1, p2, p3):
            
#                 # Get the 8 vertices of the cube
#             vertices = get_cube_vertices(x, z, cube_size/2, cube_size/2, thickness)
#             cubes_vertices.append(vertices)

# cubes_vertices = np.array(cubes_vertices)
# return cubes_vertices