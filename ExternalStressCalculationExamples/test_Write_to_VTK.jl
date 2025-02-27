using WriteVTK
using StaticArrays

# >>>> Example
"""
# Define structured grid dimensions
nx, ny = 10, 10
x = collect(range(0, 1, length=nx))  # Convert range to vector
y = collect(range(0, 1, length=ny))  # Convert range to vector
# z = zeros(nx)  # 2D plane on z = 0

# Create the structured grid and write to VTK
vtk = vtk_grid("example_structured_output.vtr", x, y )

# Add magnitudes (e.g., scalar values for each grid point)
magnitudes = x .* ones(nx, ny) # rand(nx, ny)

vtk["magnitude"] = magnitudes

# Write the file
vtk_save(vtk)

println("DONE!")
"""
# >>>> Unstructured Grid 
"""
# Define unstructured grid points as SVectors
coordinates = [
    SVector(0.0, 0.0, 0.0),
    SVector(1.0, 0.0, 0.0),
    SVector(0.0, 1.0, 0.0),
    SVector(1.0, 1.0, 0.0)
]

# Create the unstructured grid
vtk = vtk_grid("example_unstructured_output", coordinates)

# Add scalar data for each point
vtk["magnitude"] = [1.0, 2.0, 3.0, 4.0]

# Save the file
vtk_save(vtk)
"""


LoadingInputFileName="Input_Discretized.jld2"
LoadingInputExternalStressFile = "Input_ExternalStressChange.jld2"
OutputFileName = "Buijze19_50_FaultStress"

FaultCenter= load(LoadingInputFileName, "FaultCenter")

ExternalStress_Shear = load(LoadingInputExternalStressFile, "ExternalStress_Shear" )
ExternalStress_Normal = load(LoadingInputExternalStressFile, "ExternalStress_Normal" )

coordinates = [SVector(row...) for row in eachrow(FaultCenter)]

vtk = vtk_grid(OutputFileName, coordinates )
# Add scalar data for each point
vtk["ShearStress"] = ExternalStress_Shear
vtk["NormalStress"] = ExternalStress_Normal


# Save the file
vtk_save(vtk)