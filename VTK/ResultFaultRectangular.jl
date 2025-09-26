using StaticArrays
using WriteVTK
using Printf

deg2radd(θ) = (π/180) * θ

# --- same orientation helpers as before ---
function strike_dip_frame_ENZ(strike_deg::Real, dip_deg::Real; dip_variant::Symbol=:perp_to_strike)
    ϕ = deg2radd(strike_deg)
    δ = deg2radd(dip_deg)
    ŝ = SVector(cos(ϕ), sin(ϕ), 0.0)                # along strike (horizontal)
    if dip_variant === :perp_to_strike
        ĥ = SVector(sin(ϕ), -cos(ϕ), 0.0)           # right of strike
        d̂ = SVector(ĥ[1]*cos(δ), ĥ[2]*cos(δ), sin(δ))  # downwards (+z in ENZ)
    elseif dip_variant === :x_to_z
        d̂ = SVector(cos(δ), 0.0, sin(δ))
    else
        error("Unknown dip_variant = $dip_variant")
    end
    return ŝ, d̂
end

function quad_corners_ENZ(center::NTuple{3,Float64}, Ls::Float64, Ld::Float64, ŝ, d̂)
    cx, cy, cz = center
    a = 0.5 * Ls; b = 0.5 * Ld
    p1 = (cx - a*ŝ[1] - b*d̂[1], cy - a*ŝ[2] - b*d̂[2], cz - a*ŝ[3] - b*d̂[3])
    p2 = (cx + a*ŝ[1] - b*d̂[1], cy + a*ŝ[2] - b*d̂[2], cz + a*ŝ[3] - b*d̂[3])
    p3 = (cx + a*ŝ[1] + b*d̂[1], cy + a*ŝ[2] + b*d̂[2], cz + a*ŝ[3] + b*d̂[3])
    p4 = (cx - a*ŝ[1] + b*d̂[1], cy - a*ŝ[2] + b*d̂[2], cz - a*ŝ[3] + b*d̂[3])
    return (p1, p2, p3, p4)
end

# --- FIXED VTU exporter (points as 3×Np matrix; 1-based connectivity) ---
function export_fault_patches_vtu(FaultCenter::AbstractMatrix{<:Real},
                                  Angle::AbstractMatrix{<:Real},
                                  Length::AbstractMatrix{<:Real},
                                  Slip::AbstractMatrix{<:Real};
                                  outprefix::AbstractString="fault",
                                  times::AbstractVector{<:Real}=collect(1:size(Slip,1)),
                                  dip_variant::Symbol=:perp_to_strike)

    N = size(FaultCenter,1)
    T = size(Slip,1)
    @assert size(FaultCenter,2) == 3
    @assert size(Angle,2) == 2
    @assert size(Length,2) == 2
    @assert size(Angle,1) == N == size(Length,1) == size(Slip,2)
    @assert length(times) == T

    # geometry once
    quads = Vector{NTuple{4,NTuple{3,Float64}}}(undef, N)
    for i in 1:N
        cx, cy, cz = FaultCenter[i,1], FaultCenter[i,2], FaultCenter[i,3]
        strike, dip = Angle[i,1], Angle[i,2]
        Ls, Ld      = Length[i,1], Length[i,2]
        ŝ, d̂ = strike_dip_frame_ENZ(strike, dip; dip_variant=dip_variant)
        quads[i] = quad_corners_ENZ((Float64(cx), Float64(cy), Float64(cz)),
                                    Float64(Ls), Float64(Ld), ŝ, d̂)
    end

    # points as a 3×(4N) matrix (columns are XYZ)
    pts = Array{Float64}(undef, 3, 4*N)
    cells = Vector{MeshCell}(undef, N)
    for i in 1:N
        p1,p2,p3,p4 = quads[i]
        base = 4*(i-1)
        pts[:, base+1] .= (p1[1], p1[2], p1[3])
        pts[:, base+2] .= (p2[1], p2[2], p2[3])
        pts[:, base+3] .= (p3[1], p3[2], p3[3])
        pts[:, base+4] .= (p4[1], p4[2], p4[3])
        # 1-based connectivity for WriteVTK:
        cells[i] = MeshCell(VTKCellTypes.VTK_QUAD, (base+1, base+2, base+3, base+4))
    end

    # write each time step
    for (ti,t) in enumerate(times)
        fn = @sprintf("%s_%04d", outprefix, ti)    # WriteVTK will append .vtu
        vtk = vtk_grid(fn, pts, cells; compress=false)
        vtk["slip", VTKCellData()] = vec(Slip[ti, :])  # one scalar per cell
        vtk_save(vtk)
    end

    # time collection (PVD pointing to VTUs)
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
HistoryTime = Vector(ResultData["History_Time"][:,1]) 
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

export_fault_patches_vtu(FaultCenter, Angle, Length, SlipRate;
    outprefix= outdir * "fault", times=HistoryTime, dip_variant=:perp_to_strike)
   

# In ParaView:
# 1) File → Open → select "fault_patches.pvd", click Apply.
# 2) In the Color Map Editor, set Coloring to "slip".
# 3) Use the VCR controls (bottom) to play through the time steps.
