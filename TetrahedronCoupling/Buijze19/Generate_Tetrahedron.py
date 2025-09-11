import gmsh
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


def main( OutputTetrahedronVerticesFilename ):
    # Reservoir parameters
    Phi = 70 # degree
    Thickness = 2000 # meter
    Height = 200 # meter
    offset = 0.0 # 0.0 | 50.0 # meter
    
    
    # Processing the geometry   
    Phi_rad = Phi*np.pi/180

    # Get coordinates of the 8 vertices of the trapezoids
    # Right-handed
    LeftTrapsoidX1 = 0
    LeftTrapsoidX2 = 0
    LeftTrapsoidX4 = 2000-0.5*((Height+offset)/np.tan(Phi_rad))
    LeftTrapsoidX3 = LeftTrapsoidX4 + Height/np.tan(Phi_rad)
    
    BothTrapsoidYmin = 0
    BothTrapsoidYmax = Thickness

    RightTrapsoidX2 = 4000 - LeftTrapsoidX4
    RightTrapsoidX1 = RightTrapsoidX2 - Height/np.tan(Phi_rad)
    RightTrapsoidX3 = 4000
    RightTrapsoidX4 = 4000

    LeftTrapsoidZmin = - ( 3000 + offset )
    LeftTrapsoidZmax = - ( 2800 + offset )
    RightTrapsoidZmin = -3000
    RightTrapsoidZmax = -2800

    gmsh.initialize()
    gmsh.model.add("two_trapezoid_prisms")
    
    p1Left = gmsh.model.occ.addPoint(LeftTrapsoidX1, 0, LeftTrapsoidZmin)
    p2Left = gmsh.model.occ.addPoint(LeftTrapsoidX1, 0, LeftTrapsoidZmax)
    p3Left = gmsh.model.occ.addPoint(LeftTrapsoidX3, 0, LeftTrapsoidZmax)
    p4Left = gmsh.model.occ.addPoint(LeftTrapsoidX4, 0, LeftTrapsoidZmin)

    p1Right = gmsh.model.occ.addPoint(RightTrapsoidX1, 0, RightTrapsoidZmin)
    p2Right = gmsh.model.occ.addPoint(RightTrapsoidX2, 0, RightTrapsoidZmax)
    p3Right = gmsh.model.occ.addPoint(RightTrapsoidX3, 0, RightTrapsoidZmax)
    p4Right = gmsh.model.occ.addPoint(RightTrapsoidX4, 0, RightTrapsoidZmin)

    def add_trapzoid_prism(p1, p2, p3, p4, length): 
        l1 = gmsh.model.occ.addLine(p1, p2)
        l2 = gmsh.model.occ.addLine(p2, p3)
        l3 = gmsh.model.occ.addLine(p3, p4)
        l4 = gmsh.model.occ.addLine(p4, p1)

        cl = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
        surface = gmsh.model.occ.addPlaneSurface([cl])
        prism  = gmsh.model.occ.extrude([(2, surface)], 0, length, 0, numElements=[], recombine=False)
        
        return prism

    add_trapzoid_prism(p1Left, p2Left, p3Left, p4Left, Thickness)
    add_trapzoid_prism(p1Right, p2Right, p3Right, p4Right, Thickness)

    # ------------------------
    # Step 3: Mesh settings for adaptive tetrahedral mesh
    gmsh.option.setNumber("Mesh.Algorithm3D", 4)  # Frontal-Delaunay
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 100)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1000)

    # (100, 1000) --> 489
    # (200, 2000) --> 
    # (500, 1000) no

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)

    gmsh.write(f"{OutputTetrahedronVerticesFilename}.msh")
    print(f"---- MSH File Saved: {OutputTetrahedronVerticesFilename}.msh ----")
    

    tetra_vertices_coords = get_vertices(gmsh)
    save_vertices_to_txt(tetra_vertices_coords, f"{OutputTetrahedronVerticesFilename}.txt")
    print(f"---- TXT File Saved: {OutputTetrahedronVerticesFilename}.txt ----")
    
    # Finalize Gmsh
    gmsh.finalize()



if __name__ == "__main__":

    Output_Tetrahedron_Filename = "./TetrahedronCoupling/Input_Tetrahedron"
    
    main(Output_Tetrahedron_Filename)