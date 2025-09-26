# Convention: +X = East, +Y = North, +Z = Down (Left-handed)
using Printf
using CSV
using DataFrames

include("../Functions_Kuvshinov_Cuboids.jl")

function save_cuboids_to_txt(filename::AbstractString, cuboids::Vector{<:AbstractMatrix})
    open(filename, "w") do io
        println(io, "Index, CenterX, CenterY, CenterZ, LengthX, LengthY, LengthZ")
    end

    for (idx, cube) in enumerate(cuboids)
        centerx = (cube[1, 1] + cube[end, 1]) / 2
        centery = (cube[1, 2] + cube[end, 2]) / 2
        centerz = (cube[1, 3] + cube[end, 3]) / 2

        x_length = abs(cube[2, 1] - cube[1, 1])
        y_length = abs(cube[3, 2] - cube[1, 2])
        z_length = abs(cube[5, 3] - cube[1, 3])

        open(filename, "a") do io
            @printf(io, "%d, %g, %g, %g, %g, %g, %g\n",
                    idx, centerx, centery, centerz, x_length, y_length, z_length)
        end
    end
end

function get_all_cubes_vertices(inject_point, cube_side_len, cube_height, cube_num_side)
    start = -(cube_num_side/2 - 0.5) * cube_side_len
    diffs = start:cube_side_len:(start + (cube_num_side - 1) * cube_side_len)

    cz = inject_point[3]

    cubes = [Get_Cuboid_Vertices(inject_point[1] + dx,
                                inject_point[2] + dy,
                                cz,
                                cube_side_len, cube_side_len, cube_height)
             for dx in diffs, dy in diffs]

    return vec(cubes)
end


function main(OutputCubesTXTFilename)
    println("---- Generating Single Cuboid Mesh ----")

    # Reservoir parameters
    Center    = (-1000, 0, 3000.0)
    Cuboid_side_len   = 100.0
    Cuboid_side_number = 40
    Thickness = 400.0

    cuboids_all = get_all_cubes_vertices(Center, Cuboid_side_len, Thickness, Cuboid_side_number)

    save_cuboids_to_txt(OutputCubesTXTFilename, cuboids_all)

    println("---- TXT file Saved: $(OutputCubesTXTFilename)  ----")
end

OutputCubesTXTFilename = "CuboidCoupling/Input_Cuboids.txt"

main(OutputCubesTXTFilename)
