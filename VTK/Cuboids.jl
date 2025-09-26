using Printf
using CSV
using DataFrames

include("../CuboidCoupling/Functions_Kuvshinov_Cuboids.jl")

function save_cuboids_to_vtk(filename::AbstractString, cubes::Vector{<:AbstractMatrix})
    n_cubes   = length(cubes)
    n_points  = 8 * n_cubes
    n_quads   = 6 * n_cubes
    total_ints = n_quads * (4 + 1)  # (count + 4 vertex indices) per quad

    open(filename, "w") do io
        # Header
        println(io, "# vtk DataFile Version 3.0")
        println(io, "Cubes")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")

        # Points (flip Z like the PyVista call did with -center[2])
        println(io, "POINTS $n_points float")
        for cube in cubes
            for i in 1:8
                x, y, z = cube[i, 1], cube[i, 2], cube[i, 3]
                @printf(io, "%g %g %g\n", x, y, -z)
            end
        end

        # Faces as quads (6 per cube). Indices are 0-based in VTK legacy format.
        println(io, "POLYGONS $n_quads $total_ints")
        # Vertex ordering compatible with get_cube_vertices_general:
        # 0: (-x,-y,-z), 1: (+x,-y,-z), 2: (-x,+y,-z), 3: (+x,+y,-z),
        # 4: (-x,-y,+z), 5: (+x,-y,+z), 6: (-x,+y,+z), 7: (+x,+y,+z)
        quad_faces = (
            (0, 1, 3, 2),  # bottom (z-)
            (4, 5, 7, 6),  # top (z+)
            (0, 1, 5, 4),  # y-
            (2, 3, 7, 6),  # y+
            (0, 2, 6, 4),  # x-
            (1, 3, 7, 5),  # x+
        )
        for ci in 0:(n_cubes - 1)
            base = ci * 8
            for (i1, i2, i3, i4) in quad_faces
                println(io, "4 $(base+i1) $(base+i2) $(base+i3) $(base+i4)")
            end
        end
    end
end

function main(InputCuboidsFile, OutputVTKFile)
    Cuboids_Count, Cuboids_Center, Cuboids_Length = Load_Reservoir_Cuboids(InputCuboidsFile)

    Cuboids_Vertices = Vector{Matrix{Float64}}()
    for idx in 1:Cuboids_Count
        vertices = Get_Cuboid_Vertices(Cuboids_Center[idx,1],Cuboids_Center[idx,2],Cuboids_Center[idx,3], Cuboids_Length[idx,1],Cuboids_Length[idx,2],Cuboids_Length[idx,3])
        push!(Cuboids_Vertices, Matrix(vertices))
    end

    save_cuboids_to_vtk(OutputVTKFile, Cuboids_Vertices)
    # println("---- VTK file Saved: $(OutputVTKFile)  ----")
    @info "VTK file Saved: $(OutputVTKFile)"

end


InputCuboidsFile = "CuboidCoupling/Input_Cuboids.txt"
OutputVTKFile = "VTK/Input_Cuboids.vtk"

main(InputCuboidsFile, OutputVTKFile)
