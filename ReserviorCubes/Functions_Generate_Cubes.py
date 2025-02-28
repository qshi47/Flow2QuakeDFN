import numpy as np 
import pyvista as pv

# >>>>>>>>>>>>>>>> Legacy Code >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# def save_cubes_to_obj(cubes, filename):
#     """
#     Save a list of cubes (each with 8 vertices) to an OBJ file for visualization.

#     Parameters:
#     - cubes: List of cubes, where each cube is a list of 8 vertices (each vertex is a list of [x, y, z]).
#     - filename: Output OBJ file name.
#     """
#     with open(filename, "w") as obj_file:
#         obj_file.write("o Cubes\n")
#         vertex_count = 0

#         for cube in cubes:
#             # Write vertices
#             for v in cube:
#                 obj_file.write(f"v {v[0]} {v[1]} {v[2]}\n")

#             # Write faces (using 1-based indexing for OBJ format)
#             faces = [
#                 [1, 2, 4, 3],  # Bottom face
#                 [5, 6, 8, 7],  # Top face
#                 [1, 2, 6, 5],  # Side face 1
#                 [2, 4, 8, 6],  # Side face 2
#                 [4, 3, 7, 8],  # Side face 3
#                 [3, 1, 5, 7],  # Side face 4
#             ]
#             for face in faces:
#                 obj_file.write(f"f {' '.join(str(vertex_count + i) for i in face)}\n")

#             vertex_count += 8


# def add_plane_to_obj(plane_vertices, filename):
#     """
#     Add a double-sided square plane to an existing OBJ file.

#     Parameters:
#     - filename: The name of the existing OBJ file.
#     - plane_vertices: A list of 4 vertices defining the square plane.
#     """
#     with open(filename, "a") as obj_file:
#         obj_file.write("o Plane\n")
#         # Get the current number of vertices in the OBJ file
#         with open(filename, "r") as f:
#             vertex_count = sum(1 for line in f if line.startswith("v "))

#         # Write the plane vertices
#         for v in plane_vertices:
#             obj_file.write(f"v {v[0]} {v[1]} {v[2]}\n")

#         # Write the plane face (using 1-based indexing for OBJ format)
#         # Front face
#         obj_file.write(f"f {vertex_count + 1} {vertex_count + 2} {vertex_count + 3} {vertex_count + 4}\n")
#         # Back face (inverted face order for double-sided visibility)
#         obj_file.write(f"f {vertex_count + 4} {vertex_count + 3} {vertex_count + 2} {vertex_count + 1}\n")


# def save_colored_cubes_and_plane_to_obj(filename_obj, filename_mtl, cubes, plane_vertices):
#     """
#     Save cubes and a plane to an OBJ file with different colors using an MTL file.

#     Parameters:
#     - filename_obj: The name of the OBJ file.
#     - filename_mtl: The name of the MTL file.
#     - cubes: List of cubes, each with 8 vertices.
#     - plane_vertices: List of 4 vertices defining the plane.
#     """
#     # Write the MTL file to define materials
#     with open(filename_mtl, "w") as mtl_file:
#         mtl_file.write("newmtl cube_material\n")
#         mtl_file.write("Kd 1.0 0.0 0.0\n")  # Red color for cubes
#         mtl_file.write("newmtl plane_material\n")
#         mtl_file.write("Kd 0.0 0.0 1.0\n")  # Blue color for the plane

#     # Write the OBJ file
#     with open(filename_obj, "w") as obj_file:
#         obj_file.write(f"mtllib {filename_mtl}\n")

#         # Write cubes with cube_material
#         obj_file.write("usemtl cube_material\n")
#         vertex_count = 0
#         for cube in cubes:
#             for v in cube:
#                 obj_file.write(f"v {v[0]} {v[1]} {v[2]}\n")
#             faces = [
#                 [1, 2, 4, 3],
#                 [5, 6, 8, 7],
#                 [1, 2, 6, 5],
#                 [2, 4, 8, 6],
#                 [4, 3, 7, 8],
#                 [3, 1, 5, 7],
#             ]
#             for face in faces:
#                 obj_file.write(f"f {' '.join(str(vertex_count + i) for i in face)}\n")
#             vertex_count += 8

#         # Write the plane with plane_material
#         obj_file.write("usemtl plane_material\n")
#         for v in plane_vertices:
#             obj_file.write(f"v {v[0]} {v[1]} {v[2]}\n")
#         obj_file.write(f"f {vertex_count + 1} {vertex_count + 2} {vertex_count + 3} {vertex_count + 4}\n")


# def save_cubes_and_plane_to_vtk(filename, cubes, plane_vertices):
#     """
#     Save cubes and a plane to a VTK file for Paraview.

#     Parameters:
#     - filename: The name of the VTK file.
#     - cubes: List of cubes, each with 8 vertices.
#     - plane_vertices: List of 4 vertices defining the plane.
#     """
#     # Create a PyVista PolyData object for the cubes
#     cube_mesh = pv.PolyData()
#     for cube in cubes:
#         cube_mesh += pv.Cube(center=[0, 0, 0], x_length=1, y_length=1, z_length=1)

#     # Create a PyVista PolyData object for the plane
#     plane_mesh = pv.Plane(center=[0.5, 0.5, 0], i_size=1, j_size=1)

#     # Combine the meshes
#     combined_mesh = cube_mesh + plane_mesh

#     # Save the mesh as a VTK file
#     combined_mesh.save(filename)


# def save_cubes_and_plane_to_vtk2(filename, cubes, plane_vertices):
#     """
#     Save cubes and a plane to a VTK file with separate datasets for Paraview.

#     Parameters:
#     - filename: The name of the VTK file.
#     - cubes: List of cubes, each with 8 vertices.
#     - plane_vertices: List of 4 vertices defining the plane.
#     """
#     # Create a PyVista MultiBlock dataset to handle separate parts
#     multiblock = pv.MultiBlock()

#     # Add cubes to the MultiBlock dataset
#     cube_blocks = []
#     for cube in cubes:
#         cube_mesh = pv.PolyData()
#         for v in cube:
#             cube_mesh += pv.Cube(center=v, x_length=1, y_length=1, z_length=1)
#         cube_blocks.append(cube_mesh)

#     # Merge all cubes into a single block
#     cubes_mesh = pv.merge(cube_blocks)
#     multiblock["Cubes"] = cubes_mesh

#     # Create a plane mesh and add it to the MultiBlock dataset
#     plane_mesh = pv.Plane(center=[0.5, 0.5, 0], i_size=1, j_size=1)
#     multiblock["Plane"] = plane_mesh

#     # Save the MultiBlock dataset as a VTK file
#     multiblock.save(filename)


# def save_cubes_and_plane_to_ply(filename, cubes, plane_vertices):
#     """
#     Save cubes and a plane to a PLY file with different colors for Paraview.

#     Parameters:
#     - filename: The name of the PLY file.
#     - cubes: List of cubes, each with 8 vertices.
#     - plane_vertices: List of 4 vertices defining the plane.
#     """
#     vertices = []
#     faces = []
#     colors = []

#     # Add cubes to the PLY file
#     for cube in cubes:
#         start_idx = len(vertices)
#         vertices.extend(cube)
#         faces.extend([
#             [start_idx, start_idx + 1, start_idx + 3, start_idx + 2],  # Bottom face
#             [start_idx + 4, start_idx + 5, start_idx + 7, start_idx + 6],  # Top face
#             [start_idx, start_idx + 1, start_idx + 5, start_idx + 4],  # Side face 1
#             [start_idx + 1, start_idx + 3, start_idx + 7, start_idx + 5],  # Side face 2
#             [start_idx + 3, start_idx + 2, start_idx + 6, start_idx + 7],  # Side face 3
#             [start_idx + 2, start_idx, start_idx + 4, start_idx + 6],  # Side face 4
#         ])
#         colors.extend([[255, 0, 0]] * 6)  # Red color for the cubes

#     # Add plane to the PLY file
#     start_idx = len(vertices)
#     vertices.extend(plane_vertices)
#     faces.append([start_idx, start_idx + 1, start_idx + 2, start_idx + 3])
#     colors.append([0, 0, 255])  # Blue color for the plane

#     # Write the PLY file
#     with open(filename, "w") as ply_file:
#         ply_file.write("ply\n")
#         ply_file.write("format ascii 1.0\n")
#         ply_file.write(f"element vertex {len(vertices)}\n")
#         ply_file.write("property float x\n")
#         ply_file.write("property float y\n")
#         ply_file.write("property float z\n")
#         ply_file.write(f"element face {len(faces)}\n")
#         ply_file.write("property list uchar int vertex_index\n")
#         ply_file.write("property uchar red\n")
#         ply_file.write("property uchar green\n")
#         ply_file.write("property uchar blue\n")
#         ply_file.write("end_header\n")

#         # Write vertices
#         for v in vertices:
#             ply_file.write(f"{v[0]} {v[1]} {v[2]}\n")

#         # Write faces with colors
#         for i, face in enumerate(faces):
#             ply_file.write(f"{len(face)} {' '.join(map(str, face))} {colors[i][0]} {colors[i][1]} {colors[i][2]}\n")

        
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

def save_cubes_to_vtk(filename, cubes):
    """
    Save cubes and a plane to a VTK file for Paraview.

    Parameters:
    - filename: The name of the VTK file.
    - cubes: List of cubes, each with 8 vertices.
    - plane_vertices: List of 4 vertices defining the plane.
    """
    # Create a PyVista PolyData object for the cubes
    cube_mesh = pv.PolyData()
    for cube in cubes:
        cube_mesh += pv.Cube(center   = ( cube[0]+cube[-1] )/2, 
                             x_length = (cube[1]-cube[0])[0],
                             y_length = (cube[2]-cube[0])[1], 
                             z_length = (cube[4]-cube[0])[2] )

    # Create a PyVista PolyData object for the plane
    # plane_mesh = pv.Plane(center=[0.5, 0.5, 0], i_size=1, j_size=1)

    # Combine the meshes
    # combined_mesh = cube_mesh + plane_mesh

    # Save the mesh as a VTK file
    cube_mesh.save(filename)



def save_plane_to_vtk(filename, plane_center, plane_direction, plane_i_size, plane_j_size):
    """
    Save cubes and a plane to a VTK file for Paraview.

    Parameters:
    - filename: The name of the VTK file.
    - plane_i_size: Length along the X axis, passing by the center 
    """
    
    # Create a PyVista PolyData object for the plane
    plane_mesh = pv.Plane(center = plane_center, # [1, 1, 1]
                          direction = plane_direction, # [1,1,1]
                          i_size = plane_i_size, # 4
                          j_size = plane_j_size, # 2
                          i_resolution= 32, j_resolution = 32)

    # Save the mesh as a VTK file
    plane_mesh.save(filename)



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


