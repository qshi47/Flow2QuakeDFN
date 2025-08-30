import numpy as np 
import pyvista as pv
from tqdm import tqdm

def save_cubes_to_vtk(filename, cubes):
    """
    Save cubes and a plane to a VTK file for Paraview.

    Parameters:
    - filename: The name of the VTK file.
    - cubes: List of cubes, each with 8 vertices.
    """
    # Create a PyVista PolyData object for the cubes
    cube_mesh = pv.PolyData()
    for cube in tqdm(cubes, unit="Cubes"):
        center   = ( cube[0]+cube[-1] )/2

        cube_mesh += pv.Cube(center   = [ center[0], center[1], -center[2] ], 
                             x_length = np.abs( (cube[1]-cube[0])[0] ),
                             y_length = np.abs( (cube[2]-cube[0])[1] ), 
                             z_length = np.abs( (cube[4]-cube[0])[2] ) )
    # Save the mesh as a VTK file
    cube_mesh.save(filename)


def save_cubes_to_txt(filename, cubes):
    with open(filename, 'w') as f:
        f.write("Index, CenterX, CenterY, CenterZ, LengthX, LengthY, LengthZ\n")
    
    for idx, cube in enumerate(cubes):
        center   = ( cube[0]+cube[-1] )/2
        x_length = np.abs( (cube[1]-cube[0])[0] )
        y_length = np.abs( (cube[2]-cube[0])[1] ) 
        z_length = np.abs( (cube[4]-cube[0])[2] )
        with open(filename, 'a') as f:
            f.write(f"{idx+1}, {center[0]}, {center[1]}, {center[2]}, {x_length}, {y_length}, {z_length}\n")


def save_planes_to_vtk(filename, planes_center, planes_direction, planes_i_size, planes_j_size):
    """
    Save cubes and a plane to a VTK file for Paraview.

    Parameters:
    - filename: The name of the to-save VTK file.
    - planes_i_size : length along the dip direction!
    - planes_j_size : length along the strike direction!
    """
    planes_mesh = pv.PolyData()
    
    for idx_plane in tqdm( range(len(planes_center)), unit="planes"):
        planes_mesh += pv.Plane(center = planes_center[idx_plane,:], # [1, 1, 1]
                             direction = planes_direction[idx_plane], # [1,1,1]
                                i_size = planes_i_size[idx_plane], # 4
                                j_size = planes_j_size[idx_plane], # 2
                                i_resolution= 8, j_resolution = 8)

    # Save the mesh as a VTK file
    planes_mesh.save(filename)


def save_single_plane_to_vtk(filename, plane_center, plane_direction, plane_i_size, plane_j_size):
    """
    Save cubes and a plane to a VTK file for Paraview.

    Parameters:
    - filename: The name of the VTK file.
    - planes_i_size : length along the dip direction!
    - planes_i_size : length along the strike direction!
    """
    
    # Create a PyVista PolyData object for the plane
    plane_mesh = pv.Plane(center = plane_center, # [1, 1, 1]
                          direction = plane_direction, # [1,1,1]
                          i_size = plane_i_size, # 4
                          j_size = plane_j_size, # 2
                          i_resolution= 32, j_resolution = 32)

    # Save the mesh as a VTK file
    plane_mesh.save(filename)


def get_cube_vertices_general(x, y, z, cube_len_x,cube_len_y, cube_len_z):
    """
    Get the 8 vertices of a cube centered at (x, y, z) with the given size.
    """
    vertices = [
        [x - cube_len_x/2, y - cube_len_y/2,       z - cube_len_z/2],
        [x + cube_len_x/2, y - cube_len_y/2,       z - cube_len_z/2],
        [x - cube_len_x/2, y + cube_len_y/2,       z - cube_len_z/2],
        [x + cube_len_x/2, y + cube_len_y/2,       z - cube_len_z/2],
        [x - cube_len_x/2, y - cube_len_y/2,       z + cube_len_z/2],
        [x + cube_len_x/2, y - cube_len_y/2,       z + cube_len_z/2],
        [x - cube_len_x/2, y + cube_len_y/2,       z + cube_len_z/2],
        [x + cube_len_x/2, y + cube_len_y/2,       z + cube_len_z/2]
    ]
    return np.array(vertices)
