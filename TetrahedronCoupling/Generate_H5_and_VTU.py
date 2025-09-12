import h5py 
import meshio
import numpy as np
from tqdm import tqdm


def read_vertices_from_txt(filename):
    """
    Read tetrahedron vertices from a text file.
    Returns an (N, 4, 3) array where N is the number of tetrahedrons.
    """
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    N = data.shape[0] // 4
    tetra_coords = data[:, 1:].reshape(N, 4, 3)
    return tetra_coords

def transfer_vertices_to_faces(tetra_coords):
    # Define the 4 faces of a tetrahedron (Gmsh ordering)
    face_indices = [
    ([0, 1, 2], 3),
    ([0, 1, 3], 2),
    ([0, 2, 3], 1),
    ([1, 2, 3], 0)
    ]

    # Signed volume function
    def signed_tet_volume(a, b, c, d):
        sign = -1 # +1 for right-handed, -1 for left-handed
        return np.dot(sign * np.cross(b - a, c - a), d - a) / 6.0 

    # Allocate result
    N = tetra_coords.shape[0]
    tetra_faces = np.zeros((N, 4, 3, 3))

    # Loop over tetrahedra and fix face orientations
    for i, tet in enumerate(tetra_coords):
        for j, (face_idx, opp_idx) in enumerate(face_indices):
            v0, v1, v2 = tet[face_idx]
            v_opposite = tet[opp_idx]

            # Check signed volume to determine orientation
            vol = signed_tet_volume(v0, v1, v2, v_opposite)
            if vol > 0:
                # Flip face to correct outward orientation
                face = [v0, v2, v1]
            else:
                face = [v0, v1, v2]

            tetra_faces[i, j] = face
    return tetra_faces

def write_msh_to_vtu(filename_msh, filename_vtu):
    """Convert a .msh Gmsh file to .vtu format for ParaView, removing problematic cell sets."""
    # Read the MSH file
    mesh = meshio.read(filename_msh)

    # Explicitly remove cell sets from cell_data (to avoid the out-of-bounds error)
    mesh.cell_sets = {}  # Clear cell sets
    mesh.cell_data = {}  # Remove all cell data if necessary

    # Save as VTU
    meshio.write(filename_vtu, mesh)



def main(InputTetrahedronVerticesFilename, OutputTetrahedronFacesFilename):

    tetra_vertices_coords = read_vertices_from_txt(f"{InputTetrahedronVerticesFilename}.txt")
    tetra_faces_coords = transfer_vertices_to_faces(tetra_vertices_coords)

    # Save as .vtu
    write_msh_to_vtu(f"{InputTetrahedronVerticesFilename}.msh", f"{OutputTetrahedronFacesFilename}.vtu")
    print(f"---- VTU files Saved: {OutputTetrahedronFacesFilename}.vtu ----")
        
        
    # Save as .h5 
    with h5py.File(f"{OutputTetrahedronFacesFilename}.h5", 'w') as f:
        f.create_dataset('vertices_coord',  data = tetra_faces_coords)

    print(f"---- HDF5 files Saved {OutputTetrahedronFacesFilename}.h5 ----")

    print("Tetrahedrons Faces shape is ", tetra_faces_coords.shape)


if __name__ == "__main__":
    Input_Tetrahedron_Vertices_Filename = "./TetrahedronCoupling/Input_Tetrahedron"
    Output_Tetrahedron_Faces_Filename = "./TetrahedronCoupling/Input_Tetrahedron"

    main(Input_Tetrahedron_Vertices_Filename, Output_Tetrahedron_Faces_Filename)
