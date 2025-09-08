import numpy as np
import matplotlib.pyplot as plt
import h5py
from tqdm import tqdm

from Functions_Generate_Cubes import get_cube_vertices_general, save_cubes_to_vtk, save_cubes_to_txt

def get_cubes_vertices_example(inject_point, cube_side_len, cube_height, cube_num_side ):
    cube_centers_xy = []
    cubes_all = []
    cube_centers_x_diff = np.arange(-(cube_num_side/2-0.5)*cube_side_len, (cube_num_side/2+0.5)*cube_side_len, cube_side_len)
    cube_centers_y_diff = cube_centers_x_diff.copy()        
    cube_center_z = inject_point[2]

    for cube_center_x in inject_point[0] + cube_centers_x_diff:
        for cube_center_y in inject_point[1] + cube_centers_y_diff:        
            cubes_all.append( get_cube_vertices_general(cube_center_x, cube_center_y, cube_center_z, 
                                                        cube_side_len, cube_side_len, cube_height )  )
    
    return np.array(cubes_all)


if __name__ == "__main__":
    """
    # Convention: +X  = East, +Y = North, +Z = Down (Left-handed)
    """
    # -- cubes vertices --
    OutputCubesTXTFilename = "Input_Example_Cubes.txt"
    # -- for Julia -- 
    OutputCubesVerticesFilename = "Input_Example_Cubes.h5"
    # -- for Paraview --
    OutputCubesVTKFilename = "Input_Example_Cubes.vtk"

    If_save_to_VTK = True

    # Reservoir parameters
    InjectionCenter = [0, 0, 4000]
    CubeHeight = 100 # meter
    CubeSideLength = 200 # meter
    CubeSideNumber = 40
    
    cubes_all = get_cubes_vertices_example(InjectionCenter, CubeSideLength, CubeHeight, CubeSideNumber)
    print(f"Cube number : {cubes_all.shape[0]}")

    # Save to .txt file
    print(f"------ Saving TXT files {OutputCubesTXTFilename}-------")
    save_cubes_to_txt(OutputCubesTXTFilename, cubes_all)

    # Save to .vtk File
    if If_save_to_VTK == True:
        print(f"------ Saving VTK files {OutputCubesVTKFilename}-------")
        save_cubes_to_vtk(OutputCubesVTKFilename, cubes_all)

    # Save to .h5 File
    print(f"------ Save HDF5 files {OutputCubesVerticesFilename} -------")
    with h5py.File(OutputCubesVerticesFilename, 'w') as f:
        f.create_dataset('cubes',  data = cubes_all)

    print(">>>> Files All Saved >>>>")