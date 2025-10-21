using DelimitedFiles
using JLD2
using Interpolations
using ProgressBars
using WriteVTK
using Printf

include("../CuboidCoupling/Functions_Kuvshinov_Cuboids.jl")


InputCuboidsFile = "CuboidCoupling/Input_Cuboids.txt"
TemperatureFile = "CuboidCoupling/Input_Temperature.txt"
ResultFile = "Results/Result.jld2"


# load cuboids
Cuboids_Count, Cuboids_Center, Cuboids_Length = Load_Reservoir_Cuboids(InputCuboidsFile)
Cuboids_Center[:,3] .= -Cuboids_Center[:,3]

# Cuboids_Vertices = Vector{Matrix{Float64}}()
# for idx in 1:Cuboids_Count
#     vertices = Get_Cuboid_Vertices(Cuboids_Center[idx,1],Cuboids_Center[idx,2],Cuboids_Center[idx,3], Cuboids_Length[idx,1],Cuboids_Length[idx,2],Cuboids_Length[idx,3])
#     push!(Cuboids_Vertices, Matrix(vertices))
# end

# load rupture time array
ResultData = load(ResultFile)
HistoryTime = Vector( ResultData["History_Time"][:,1] )

# load external time array
ExternalStress_TimeArray = Float64.(readdlm(TemperatureFile, ',' )[1,2:end])

# load temperature data
TemperatureChange = Float64.(readdlm(TemperatureFile, ',')[2:end,2:end]')

# interpolate temperature to rupture time array
Temperature_Interpolated = zeros(length(HistoryTime), Cuboids_Count)
@info "Interpolating temperature for $(length(HistoryTime)) time steps"

for idx_Cuboid in tqdm( 1:Cuboids_Count, unit = "cuboid" )
    xs = ExternalStress_TimeArray
    ys = TemperatureChange[:, idx_Cuboid]
    interp_linear = linear_interpolation(xs, ys)
    Temperature_Interpolated[:, idx_Cuboid] = interp_linear(HistoryTime)
end






function write_cuboid_series_vtu(
    centers::AbstractMatrix{<:Real},
    halfsizes::AbstractMatrix{<:Real},
    temperature::AbstractMatrix{<:Real},
    time::AbstractVector{<:Real};
    outdir::AbstractString = "cuboids_vtu",
    compress::Bool = true,
    write_pvd::Bool = true
)
    # ---- checks ----
    nc = size(centers, 1)
    @assert size(centers, 2) == 3           "centers must be (nc, 3)"
    @assert size(halfsizes) == (nc, 3)      "halfsizes must be (nc, 3)"
    nt, nc_p = size(temperature)
    @assert nc_p == nc                       "temperature must be (nt, nc)"
    @assert length(time) == nt               "time length must equal number of time steps"
    
    isdir(outdir) || mkpath(outdir)
    if isdir(outdir)
        for f in readdir(outdir)
            rm(joinpath(outdir, f); force=true)
        end
    end

    # ---- storage ----
    npts   = 8 * nc
    points = Matrix{Float64}(undef, 3, npts)
    cells  = Vector{MeshCell{VTKCellTypes.VTKCellType, Vector{Int}}}(undef, nc)

    # write a column (no broadcasting/tuples)
    @inline function setcol!(A::AbstractMatrix{<:Real}, j::Int, x::Float64, y::Float64, z::Float64)
        @inbounds begin
            A[1, j] = x
            A[2, j] = y
            A[3, j] = z
        end
    end

    # VTK_HEXAHEDRON node order for the 8 corners:
    # 0:(-x,-y,-z) 1:(+x,-y,-z) 2:(+x,+y,-z) 3:(-x,+y,-z)
    # 4:(-x,-y,+z) 5:(+x,-y,+z) 6:(+x,+y,+z) 7:(-x,+y,+z)

    # ---- build geometry ----
    @inbounds for i in 1:nc
        cx = Float64(centers[i,1]); cy = Float64(centers[i,2]); cz = Float64(centers[i,3])
        hx = Float64(halfsizes[i,1]); hy = Float64(halfsizes[i,2]); hz = Float64(halfsizes[i,3])
        @assert hx > 0 && hy > 0 && hz > 0 "halfsizes must be positive; got ($hx,$hy,$hz) for cell $i"

        base = 8*(i-1)

        # z- (bottom)
        setcol!(points, base+1, cx - hx, cy - hy, cz - hz)
        setcol!(points, base+2, cx + hx, cy - hy, cz - hz
        )
        setcol!(points, base+3, cx + hx, cy + hy, cz - hz)
        setcol!(points, base+4, cx - hx, cy + hy, cz - hz)

        # z+ (top)
        setcol!(points, base+5, cx - hx, cy - hy, cz + hz)
        setcol!(points, base+6, cx + hx, cy - hy, cz + hz)
        setcol!(points, base+7, cx + hx, cy + hy, cz + hz)
        setcol!(points, base+8, cx - hx, cy + hy, cz + hz)

        # connectivity for this hex: 8 consecutive point indices
        conn = collect(base .+ (1:8))  # Vector{Int}
        cells[i] = MeshCell(VTKCellTypes.VTK_HEXAHEDRON, conn)
    end

    # ---- write one VTU per time step ----
    vtu_files = String[]
    @info "Writing VTU files"
    @inbounds for it in tqdm( 1:nt, unit = "time step" )
        fname = joinpath(outdir, @sprintf("cuboids_%04d.vtu", it))
        push!(vtu_files, fname)

        vtk_grid(fname, points, cells; compress=compress) do vtk
            vtk["temperature", VTKCellData()] = Float64.(view(temperature, it, :))
            vtk["TimeValue", VTKFieldData()] = [Float64(time[it])]
        end
    end
    @info "Wrote $(nt) VTU files in '$outdir'."

    # ---- optional .pvd collection ----
    if write_pvd
        pvdfile = joinpath(outdir, "series.pvd")
        open(pvdfile, "w") do io
            println(io, """<?xml version="1.0"?>""")
            println(io, """<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">""")
            println(io, """  <Collection>""")
            for (it, f) in enumerate(vtu_files)
                @printf(io, "    <DataSet timestep=\"%.16g\" group=\"\" part=\"0\" file=\"%s\"/>\n",
                        Float64(time[it]), basename(f))
            end
            println(io, """  </Collection>""")
            println(io, """</VTKFile>""")
        end
        @info "Wrote collection '$pvdfile'."
    end

    return nothing
end

Year_in_Second = 365*24*3600

# interpolated
write_cuboid_series_vtu(Cuboids_Center, Cuboids_Length./2, Temperature_Interpolated, HistoryTime/Year_in_Second; outdir="VTK/CuboidsTemperature_vtu")

# non-interpolated (original)
# write_cuboid_series_vtu(Cuboids_Center, Cuboids_Length./2, TemperatureChange, ExternalStress_TimeArray, outdir="VTK/CuboidsTemperature_vtu")