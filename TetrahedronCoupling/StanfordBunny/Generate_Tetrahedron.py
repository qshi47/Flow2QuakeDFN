import gmsh, math, os, tempfile
import numpy as np


def get_vertices(gmsh):
    # Get all node coordinates
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = np.array(node_coords).reshape(-1, 3)

    # Get tetrahedral elements (type 4)
    tetra_tags, tetra_node_tags = gmsh.model.mesh.getElementsByType(4)
    tetra_node_tags = np.array(tetra_node_tags).reshape(-1, 4)

    # Create a mapping from node tags to coordinates
    node_dict = {tag: coord for tag, coord in zip(node_tags, node_coords)}

    # Extract tetrahedral coordinates
    tetra_coords = np.array([[node_dict[tag] for tag in tetra] for tetra in tetra_node_tags]) # (N,4,3)

    return tetra_coords
    

def save_vertices_to_txt(tetra_coords, filename):
    """
    tetra_coords: (N, 4, 3) array: coords of the 4 vertices of N tetrahedrons
    """
    N = tetra_coords.shape[0]
    with open(filename, "w") as f:
        f.write(f"Index, X, Y, Z \n")
        for i in range(N):
            for j in range(4):
                f.write(f"{i+1}, {tetra_coords[i, j, 0]}, {tetra_coords[i, j, 1]}, {tetra_coords[i, j, 2]}\n")


stl = "./TetrahedronCoupling/StanfordBunny/Stanford_Bunny.stl" 
lc  = 100.0  # target element size

gmsh.initialize()
gmsh.model.add("bunny")
gmsh.merge(stl)

# Clean obvious import artifacts
gmsh.model.mesh.removeDuplicateNodes([])
gmsh.model.mesh.removeDuplicateElements([])

# Build topology from the discrete surface (no reparametrization)
gmsh.model.mesh.createTopology()

# Create one closed volume bounded by all surface patches
surf_tags = [s[1] for s in gmsh.model.getEntities(2)]
sl = gmsh.model.geo.addSurfaceLoop(surf_tags)
gmsh.model.geo.addVolume([sl])
gmsh.model.geo.synchronize()

# 3D tetra meshing (pick one algorithm; 4=Frontal/Netgen, 1=Delaunay)
gmsh.option.setNumber("Mesh.Algorithm3D", 4)
# Control size (tune to taste)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 25000)

gmsh.model.mesh.generate(3)

gmsh.write("./TetrahedronCoupling/Input_Tetrahedron.msh")   # Gmsh v4 format
# gmsh.fltk.run()             # uncomment to view


tetra_vertices_coords = get_vertices(gmsh)
save_vertices_to_txt(tetra_vertices_coords, "./TetrahedronCoupling/Input_Tetrahedron.txt")

gmsh.finalize()