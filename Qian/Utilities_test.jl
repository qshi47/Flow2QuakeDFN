using DelimitedFiles
using JLD2
include("Utilities_GRG.jl")

# println("[GO!!!]")

# BLOCKS_centers_pos = readdlm("./Qian/GRG_BLOCK_center_pos.txt") # (6351, 3)
# PRESSURE_pos = readdlm("./Qian/GRG_PRESSURE_pos.txt")
# BLOCKS = load("./Qian/GRG_BLOCK_vertices.jld2", "BLOCKS")

# Correspondense_Matrix = Calculate_Correspondense_Matrix(BLOCKS_centers_pos, PRESSURE_pos)

BLOCKS_center_pos_XYZ = BuildReservoirGeometry_44([0,0,1950], 400, 100)