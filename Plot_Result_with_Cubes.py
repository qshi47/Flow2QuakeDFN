# %%
import h5py

import os
import numpy as np 

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.animation import FuncAnimation, PillowWriter

matplotlib.rcParams.update({
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

LINEWIDTH3D = 0.5

FolderName = "./Results/demo" 
#  ExternalStressCalculationExamples/Plates
MovieName = "1_Movie"

def delete_folder_contents(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                os.rmdir(file_path)
        except Exception as e:
            print(f'Failed to delete {file_path}. Reason: {e}')

delete_folder_contents(FolderName)


with h5py.File("./Results/Result.jld2") as file:
    # print(file.keys())
    History_Time = np.array(file["History_Time"]).T
    History_Time = History_Time.reshape( (len(History_Time), ) )
    History_V = np.array( file["History_V"] ).T
    History_NormalStress = np.array(file["History_NormalStress"]).T
    History_Disp = np.array(file["History_Disp"]).T
    History_Theta = np.array(file["History_Theta"]).T

# %%
with h5py.File("./Input_ExternalStressChange.jld2") as file:
    # print(file.keys())
    ExternalStress_TimeArray = np.array(file["ExternalStress_TimeArray"])
    Cm_Delta_P = np.array(file["Cm_Delta_P"]).T
    Delta_P = np.array( file["Delta_P"] ).T
    ExternalStress_Normal = np.array(file["ExternalStress_Normal"]).T
    ExternalStress_Shear = np.array(file["ExternalStress_Shear"]).T
    # Disp_1 = np.array( file["Disp_1"] ).T
    # Disp_2 = np.array( file["Disp_2"] ).T
    # Disp_3 = np.array( file["Disp_3"] ).T
    # Sigma_11 = np.array( file["SIGMA_11"] ).T
    # Sigma_22 = np.array( file["SIGMA_22"] ).T
    # Sigma_33 = np.array( file["SIGMA_33"] ).T
    # Sigma_12 = np.array( file["SIGMA_12"] ).T
    # Sigma_13 = np.array( file["SIGMA_13"] ).T
    # Sigma_23 = np.array( file["SIGMA_23"] ).T

# %%
with h5py.File("./Input_Discretized.jld2") as file:
    # print(file.keys())
    FaultCenter = np.array( file["FaultCenter"] ).T
    FaultLengthStrike = np.array(file["FaultLengthStrike"]) 
    FaultLengthDip = np.array(file["FaultLengthDip"])
    FaultStrikeAngle = np.array(file["FaultStrikeAngle"]) 
    FaultDipAngle = np.array(file["FaultDipAngle"]) 
    FaultLLRR = np.array(file["FaultLLRR"])


# %%
with h5py.File("./Input_Reservoir_Vertices.h5", 'r') as f:
    # Access the matrices group
    matrices_group = f['matrices']
    
    # Read all matrices
    Reservoir_BLOCKS_Vertices_Pos = np.array( [matrices_group[key][:] for key in matrices_group.keys()] )

    # Transport (N, 3, 8) to (N, 8, 3)
    Reservoir_BLOCKS_Vertices_Pos = Reservoir_BLOCKS_Vertices_Pos.swapaxes(1,2)

# %%
with h5py.File("./Input_Reservoir_Center.h5", 'r') as f:
    # Access the matrices group
    vectors_group = f['vectors']
    
    # Read all matrices
    Reservoir_BLOCKS_Origin_Pos = np.array( [vectors_group[key][:] for key in vectors_group.keys()] )

# %%
# with h5py.File("../Input_Reservoir_Shape.jld2") as file:
#     print(file.keys())
#     Reservoir_BLOCKS_Origin_Pos = np.array(file["Reservoir_BLOCKS_Origin_Pos"])
#     Reservoir_BLOCKS_Vertices_Pos = np.array(file["Reservoir_BLOCKS_Vertices_Pos"])

def CalculatebPlotBoundaries():
    # Calculate buffer and plot boundaries
    buffer = np.max(FaultLengthStrike) / 3
    xmin = np.min(FaultCenter[:, 0]) - buffer
    xmax = np.max(FaultCenter[:, 0]) + buffer
    ymin = np.min(FaultCenter[:, 1]) - buffer
    ymax = np.max(FaultCenter[:, 1]) + buffer
    zmax = 0
    zmin = -np.max(FaultCenter[:, 2]) - buffer
    
    xCenter = (xmin + xmax) / 2
    yCenter = (ymin + ymax) / 2
    zCenter = (zmin + zmax) / 2
    
    maxlength = np.max([xmax - xmin, ymax - ymin, zmin - zmax])
    xmin = xCenter - maxlength / 2
    xmax = xCenter + maxlength / 2
    ymin = yCenter - maxlength / 2
    ymax = yCenter + maxlength / 2
    zmin = zmax - maxlength - buffer

    return xmin, xmax, ymin, ymax, zmin, zmax

    
# %%
def FaultPlot_3D_v2(fig, ax, FaultStress4Plot, sm, FaultCenter, 
                    FaultLengthStrike, FaultLengthDip, FaultStrikeAngle, FaultDipAngle, FaultLLRR):
    """
    Translated from Julia to Python by Qian Shi 2024/10/28
    Orinigally in Julia by KJ
    """
    # Use the colormap from the scalar mappable
    cmap = sm.cmap

    # Get color scaling values
    MaxValue = sm.norm.vmax
    MinValue = sm.norm.vmin

    # Iterate through fault indices
    for FaultIdx in range(len(FaultLengthStrike)):
        # Rotation matrices (converting from Julia's cosd/sind to numpy)
        RotMatStrike = np.array([
            [np.cos(np.deg2rad(FaultStrikeAngle[FaultIdx])), -np.sin(np.deg2rad(FaultStrikeAngle[FaultIdx])), 0],
            [np.sin(np.deg2rad(FaultStrikeAngle[FaultIdx])), np.cos(np.deg2rad(FaultStrikeAngle[FaultIdx])), 0],
            [0, 0, 1]
        ])
        
        RotMatDip = np.array([
            [1, 0, 0],
            [0, np.cos(np.deg2rad(FaultDipAngle[FaultIdx])), -np.sin(np.deg2rad(FaultDipAngle[FaultIdx]))],
            [0, np.sin(np.deg2rad(FaultDipAngle[FaultIdx])), np.cos(np.deg2rad(FaultDipAngle[FaultIdx]))]
        ])
        
        # Calculate vertex points (Julia's matrix multiplication equivalent in Python)
        p1 = RotMatStrike @ RotMatDip @ np.array([FaultLengthStrike[FaultIdx]/2, -FaultLengthDip[FaultIdx]/2, 0]) + \
            np.array([FaultCenter[FaultIdx, 0], FaultCenter[FaultIdx, 1], -FaultCenter[FaultIdx, 2]])
        
        p2 = RotMatStrike @ RotMatDip @ np.array([-FaultLengthStrike[FaultIdx]/2, -FaultLengthDip[FaultIdx]/2, 0]) + \
            np.array([FaultCenter[FaultIdx, 0], FaultCenter[FaultIdx, 1], -FaultCenter[FaultIdx, 2]])
        
        p3 = RotMatStrike @ RotMatDip @ np.array([-FaultLengthStrike[FaultIdx]/2, FaultLengthDip[FaultIdx]/2, 0]) + \
            np.array([FaultCenter[FaultIdx, 0], FaultCenter[FaultIdx, 1], -FaultCenter[FaultIdx, 2]])
        
        p4 = RotMatStrike @ RotMatDip @ np.array([FaultLengthStrike[FaultIdx]/2, FaultLengthDip[FaultIdx]/2, 0]) + \
            np.array([FaultCenter[FaultIdx, 0], FaultCenter[FaultIdx, 1], -FaultCenter[FaultIdx, 2]])
        
        # Create polygon vertices
        verts = [tuple(map(tuple, [p1, p2, p3, p4]))]
        
        # Create 3D polygon collection
        p3c = Poly3DCollection(verts, linewidths=1, alpha=0.8)
        
        # Add to the plot
        collection = ax.add_collection3d(p3c)
        collection.set_linewidth(LINEWIDTH3D)  # Thinner lines
        collection.set_antialiased(True) 

        
        # Calculate color based on normalized stress value
        """
        PlotValue = (FaultStress4Plot[FaultIdx] - MinValue) / (MaxValue - MinValue)
        """
        PlotValue = sm.norm( FaultStress4Plot[FaultIdx] )
        
        # Get face color from colormap
        face_color = cmap(PlotValue)
        
        # Edge color (light gray)
        edge_color = [128/256, 128/256, 128/256] # grey
        
        # Set polygon colors
        p3c.set_facecolor(face_color)
        p3c.set_edgecolor(edge_color)

    return fig, ax

# %%
def Plot_Stress_on_Fault(figX, axX, StressForPlot):
    sm = plt.cm.ScalarMappable( cmap = "Spectral", 
                                norm=plt.Normalize(
                                vmin =  -3e7, # np.quantile(StressForPlot, 0.02), 
                                vmax =  +1e3 )  # np.quantile(StressForPlot, 0.98) )
                                )

    FaultPlot_3D_v2(figX, axX, StressForPlot, sm, FaultCenter, FaultLengthStrike, FaultLengthDip,
    FaultStrikeAngle, FaultDipAngle, FaultLLRR)

    cbar = plt.colorbar(sm, ax = axX, location="right", shrink=0.8, pad=0.25, aspect = 20)

    return figX, axX, cbar


def Plot_BLOCKS(figX, axX, BLOCKS, linecolor = "black"):
    for block in BLOCKS:
        edges = [
            (0, 1), (1, 3), (3, 2), (2, 0),  # bottom square 
            (4, 5), (5, 7), (7, 6), (6, 4),  # top square 
            (0, 4), (1, 5), (2, 6), (3, 7)   # vertical edges
        ]
        x = block[:, 0]
        y = block[:, 1]
        z = -block[:, 2]

        for (i, j) in edges:
            axX.plot([x[i], x[j]], [y[i], y[j]], [z[i], z[j]], color = linecolor, alpha = 0.5)
        
    
    return figX, axX


def Plot_BLOCKS_with_Pressure(figX, axX, delta_pressure, blocks ):
    ColorMap_cube = matplotlib.colormaps.get_cmap("coolwarm")  # cm.get_cmap("coolwarm")
    ColorNorm_cube = colors.Normalize(vmin = -3,  # np.min(delta_pressure), 
                                      vmax = 8 ) # np.max(delta_pressure))
    sm_cube = plt.cm.ScalarMappable(cmap = "coolwarm", norm = ColorNorm_cube)
    
    ColorBar_cube = plt.colorbar(sm_cube, ax = axX, location="left", shrink=0.8, pad=0.15, aspect = 20)
    
    for (idxBlock,vertices) in enumerate(blocks):
        vertices_np = np.array(vertices)
        vertices_np[:,2] *= -1

        # Define faces (adjust indexing)
        faces = [
            vertices_np[[0,1,3,2], :],  # bottom face
            vertices_np[[4,5,7,6], :],  # top face
            vertices_np[[0,1,5,4], :],  # front face
            vertices_np[[1,3,7,5], :],  # right face
            vertices_np[[3,2,6,7], :],  # back face
            vertices_np[[2,0,4,6], :]   # left face
        ]
        

        cube = Poly3DCollection(faces, 
                                facecolors = ColorMap_cube( ColorNorm_cube(delta_pressure[idxBlock]) ), 
                                edgecolors = "black", 
                                alpha = 0.1)
        
        collection = axX.add_collection3d(cube)
        collection.set_linewidth(LINEWIDTH3D)  # Thinner lines
        collection.set_antialiased(True)
    
    return figX, axX, ColorBar_cube



# %%
def create_movie_from_images(image_files):
    """
    Convert image files to a movie
    
    Parameters:
    -----------
    image_files : list
        List of image filenames to include in the movie
    """
    # Sort filenames to ensure correct order
    image_files.sort()
    
    # Read first image to get dimensions
    first_image = plt.imread(image_files[0])
    
    # Create figure and axis for animation
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None) # remove the blank

    # Initialize plot
    im = ax.imshow(first_image)
    ax.axis('off')
    
    # Animation update function
    def update(frame):
        im.set_array(plt.imread(image_files[frame]))
        return [im]
    
    # Create animation
    anim = FuncAnimation(
        fig, 
        update, 
        frames=len(image_files), 
        interval=200,  # 200ms between frames
        blit=True
    )
    
    # Save as movie
    writer = PillowWriter(fps=5)
    anim.save(f'{FolderName}/{MovieName}.gif',  writer=writer)
    
    """
    # Clean up temporary image files
    for filename in image_files:
         os.remove(filename)
    """

    plt.close(fig)


def generate_single_figure(time_idx):
    """
    Generate a single figure
    
    Parameters:
    -----------
    time_idx : int
        Index of the figure to generate
    
    Returns:
    --------
    filename : str
        Path to the saved figure
    """
    # Create a unique filename
    filename = f"{FolderName}/fig_{time_idx:04d}.png"

    fig3 = plt.figure( figsize=(12, 8) ) # figsize=(20, 15), dpi=300
    ax3 = fig3.add_subplot(111, projection="3d")

    
    fig3, ax3, cbar3 = Plot_Stress_on_Fault(fig3, ax3, ExternalStress_Normal[time_idx,:]) #  np.log10(History_V[time_idx,:]) ExternalStress_Normal[time_idx,:]
    # fig3, ax3, cbar_cubes = Plot_BLOCKS_with_Pressure(fig3, ax3, np.log10( 1e-10 + Delta_P[time_idx,:] ), Reservoir_BLOCKS_Vertices_Pos) 
    fig3, ax3 = Plot_BLOCKS(fig3, ax3, Reservoir_BLOCKS_Vertices_Pos)
    
    
    xmin, xmax, ymin, ymax, zmin, zmax = CalculatebPlotBoundaries()
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(ymin, ymax)
    ax3.set_zlim(zmin, zmax)

    ax3.set_title( f"Time = {ExternalStress_TimeArray[time_idx]/(24*3600):.2e} Days" ) # History_Time
    ax3.set_xlabel("X (m)")
    ax3.set_ylabel("Y (m)")
    ax3.set_zlabel("Z (m)")

    # cbar_cubes.set_label("log10 (Pressure Change) (Pa)") # cbar.ax.set_ylabel("Pressure Change(Pa)")
    cbar3.set_label(" External Normal Strress (Pa) ") # External Normal Strress (Pa)
    

    ax3.view_init(elev=30, azim=-45, roll=0)
    # plt.show()

    fig3.savefig(filename,
                 bbox_inches='tight')  # Remove extra white space
    plt.close()
    
    return filename


def create_movie_parallel(TimeIdx4Plot):
    """
    Create a movie by generating figures in parallel
    """
    num_cores = 16 #cpu_count()
    
    # Create a process pool
    # 'filenames' is a list consisting of file string
    """
    # Method1: MultiProcessing
    with Pool( processes = num_cores ) as pool:
        # Generate figures in parallel
        filenames = pool.map(generate_single_figure, TimeIdx4Plot)
    """
    # Method2: tqdm
    filenames = process_map(generate_single_figure, TimeIdx4Plot,  max_workers = num_cores)
    
    # figures are already generated, then create the movies 
    print("----Plot Figures Done, Now Making Movies----")
    create_movie_from_images(filenames)


def main():
    create_movie_parallel( range( 0,  len(History_Time) ) ) # len(History_Time) len(ExternalStress_TimeArray)
    print(f"Movie created successfully: {MovieName}.gif")


if __name__ == '__main__':
    main()