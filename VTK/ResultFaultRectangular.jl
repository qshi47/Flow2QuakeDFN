using StaticArrays
using WriteVTK
using Printf
using JLD2
using ProgressBars

deg2radd(θ) = (π/180) * θ

# ENU (z up). Strike from +x→+y (CCW in EN-plane).
# dip_variant:
#   :perp_to_strike => dip in the vertical plane ⟂ to strike, dipping downward (−z)
#   :x_to_z         => dip in the x–z plane from +x toward −z (independent of strike)
function strike_dip_frame_ENU(strike_deg::Real, dip_deg::Real; dip_variant::Symbol=:perp_to_strike)
    ϕ = deg2radd(strike_deg)
    δ = deg2radd(dip_deg)

    ŝ = SVector(cos(ϕ), sin(ϕ), 0.0)  # along strike (horizontal)

    if dip_variant === :perp_to_strike
        ĥ = SVector(sin(ϕ), -cos(ϕ), 0.0)          # horizontal, right of strike
        d̂ = SVector(ĥ[1]*cos(δ), ĥ[2]*cos(δ), -sin(δ))  # tilt downward (−z in ENU)
    elseif dip_variant === :x_to_z
        d̂ = SVector(cos(δ), 0.0, -sin(δ))          # +x toward −z
    else
        error("Unknown dip_variant = $dip_variant")
    end
    return ŝ, d̂
end

# Build quad corners from center (E,N,Z), lengths along strike (Ls) and down-dip (Ld)
function quad_corners(center::NTuple{3,Float64}, Ls::Float64, Ld::Float64, ŝ, d̂)
    cx, cy, cz = center
    a = 0.5 * Ls; b = 0.5 * Ld
    p1 = (cx - a*ŝ[1] - b*d̂[1], cy - a*ŝ[2] - b*d̂[2], cz - a*ŝ[3] - b*d̂[3])
    p2 = (cx + a*ŝ[1] - b*d̂[1], cy + a*ŝ[2] - b*d̂[2], cz + a*ŝ[3] - b*d̂[3])
    p3 = (cx + a*ŝ[1] + b*d̂[1], cy + a*ŝ[2] + b*d̂[2], cz + a*ŝ[3] + b*d̂[3])
    p4 = (cx - a*ŝ[1] + b*d̂[1], cy - a*ŝ[2] + b*d̂[2], cz - a*ŝ[3] + b*d̂[3])
    return (p1, p2, p3, p4)
end

function export_fault_patches_vtu_ENU(FaultCenter::AbstractMatrix{<:Real},
                                      Angle::AbstractMatrix{<:Real},
                                      Length::AbstractMatrix{<:Real},
                                      data::AbstractMatrix{<:Real};
                                      outprefix::AbstractString="fault_ENU",
                                      times::AbstractVector{<:Real}=collect(1:size(data,1)),
                                      dip_variant::Symbol=:perp_to_strike)

    N = size(FaultCenter,1)
    T = size(data,1)
    @assert size(FaultCenter,2) == 3
    @assert size(Angle,2) == 2
    @assert size(Length,2) == 2
    @assert size(Angle,1) == N == size(Length,1) == size(data,2)
    @assert length(times) == T

    # Geometry (time-invariant)
    quads = Vector{NTuple{4,NTuple{3,Float64}}}(undef, N)
    for i in 1:N
        cx, cy, cz_ENZ = FaultCenter[i,1], FaultCenter[i,2], FaultCenter[i,3]
        cz = -Float64(cz_ENZ)  # flip sign: ENZ → ENU
        strike, dip = Angle[i,1], Angle[i,2]
        Ls, Ld      = Length[i,1], Length[i,2]
        ŝ, d̂ = strike_dip_frame_ENU(strike, dip; dip_variant=dip_variant)
        quads[i] = quad_corners((Float64(cx), Float64(cy), cz), Float64(Ls), Float64(Ld), ŝ, d̂)
    end

    # Points as 3×(4N) matrix, 1-based connectivity
    pts = Array{Float64}(undef, 3, 4*N)
    cells = Vector{MeshCell}(undef, N)
    for i in 1:N
        p1,p2,p3,p4 = quads[i]
        b = 4*(i-1)
        pts[:, b+1] .= (p1[1], p1[2], p1[3])
        pts[:, b+2] .= (p2[1], p2[2], p2[3])
        pts[:, b+3] .= (p3[1], p3[2], p3[3])
        pts[:, b+4] .= (p4[1], p4[2], p4[3])
        cells[i] = MeshCell(VTKCellTypes.VTK_QUAD, (b+1, b+2, b+3, b+4))
    end

    # Time steps
    @info "Writing $(T) time steps to VTU files with prefix '$outprefix'."
    
    for (ti,t) in tqdm( enumerate(times), unit="Time steps")
        fn = @sprintf("%s_%04d", outprefix, ti)  # .vtu added by WriteVTK
        vtk = vtk_grid(fn, pts, cells; compress=false)
        vtk["SlipRate", VTKCellData()] = vec(data[ti, :])
        vtk["TIME", VTKFieldData()] = [t]
        vtk_save(vtk)
    end

    # PVD collection
    open(outprefix * ".pvd", "w") do io
        println(io, """<?xml version="1.0"?>""")
        println(io, """<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">""")
        println(io, "  <Collection>")
        for (ti,t) in enumerate(times)
            @printf(io, "    <DataSet timestep=\"%.10g\" group=\"\" part=\"0\" file=\"%s_%04d.vtu\"/>\n",
                    Float64(t), outprefix, ti)
        end
        println(io, "  </Collection>")
        println(io, "</VTKFile>")
    end
    return nothing
end



ResultFile = "Results/Result.jld2"
FaultFile = "Input_Discretized.jld2"

ResultData = load(ResultFile)
HistoryTime = Vector(ResultData["History_Time"][:,1]) ./ (365*24*3600)
SlipRate =  ResultData["History_V"]

FaultCenter = load(FaultFile, "FaultCenter")
FaultStrikeAngle = load(FaultFile, "FaultStrikeAngle")
FaultDipAngle = load(FaultFile, "FaultDipAngle")
FaultLengthStrike = load(FaultFile, "FaultLengthStrike")
FaultLengthDip = load(FaultFile, "FaultLengthDip")

Angle = hcat(FaultStrikeAngle, FaultDipAngle)
Length = hcat(FaultLengthStrike, FaultLengthDip)

outdir = "VTK/Fault_Recangular_vtu/"
isdir(outdir) || mkpath(outdir)
if isdir(outdir)
    for f in readdir(outdir)
        rm(joinpath(outdir, f); force=true)
    end
end

export_fault_patches_vtu_ENU(FaultCenter, Angle, Length, SlipRate;
    outprefix= outdir * "fault", times=HistoryTime, dip_variant=:perp_to_strike)
   
