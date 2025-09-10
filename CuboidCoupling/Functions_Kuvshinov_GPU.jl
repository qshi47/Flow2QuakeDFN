using CUDA
using StaticArrays


@inline function fFunc(x::T, y::T, z::T, R::T) where {T<:AbstractFloat}
    tiny = eps(T)
    denom = z*R
    denom = ifelse(denom == zero(T), (denom >= zero(T) ? tiny : -tiny), denom)
    ratio = (x*y) / denom

    a1 = abs(R + y); a1 = ifelse(a1 > tiny, a1, tiny)
    a2 = abs(R + x); a2 = ifelse(a2 > tiny, a2, tiny)

    return z*atan(ratio) - x*log(a1) - y*log(a2)
end

const Sσ = SVector{8, Int8}(-1, 1, 1, -1, 1, -1, -1, 1)

# --------------------------
# Pack Vector{Matrix{T}} -> Array{T,3} of size (8,3,N)
# --------------------------
function pack_vertices(V::Vector{<:AbstractMatrix{T}}) where {T}
    N = length(V)
    @assert N > 0 "Cuboids_Vertices vector is empty"
    out = Array{T,3}(undef, 8, 3, N)
    @inbounds for j in 1:N
        M = V[j]
        @assert size(M,1) == 8 && size(M,2) == 3 "Each cuboid matrix must be (8,3), got $(size(M)) at j=$j"
        @views out[:,:,j] .= M
    end
    return out
end

# --------------------------
# One thread computes Disp(3) & Sigma(6) for a single (i,j)
# --------------------------
function kernel_porodispstress!(
    outDisp,                 # CuDeviceArray{T,3}   (3, FaultCount, Cuboids_Count)
    outSig,                  # CuDeviceArray{T,3}   (6, FaultCount, Cuboids_Count)
    FaultCenter,             # CuDeviceMatrix{T}    (FaultCount, 3)
    Cuboids_Vertices,        # CuDeviceArray{T,3}   (8, 3, Cuboids_Count)
    Cm_per_cuboid,           # CuDeviceVector{T}    (Cuboids_Count)
    shear::T                 # scalar T
) where {T<:AbstractFloat}

    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    j = (blockIdx().y-1) * blockDim().y + threadIdx().y

    FC = size(FaultCenter, 1)
    CC = size(Cuboids_Vertices, 3)
    if i > FC || j > CC
        return
    end

    @inbounds begin
        x0 = FaultCenter[i,1];  y0 = FaultCenter[i,2];  z0 = FaultCenter[i,3]

        cm = Cm_per_cuboid[j]
        πT = T(π)
        pref_disp = cm / (4 * πT)
        pref_sig  = cm * shear / (2 * πT)
        tiny = eps(T)

        # accumulators
        d1 = zero(T); d2 = zero(T); d3 = zero(T)
        s11 = zero(T); s22 = zero(T); s33 = zero(T)
        s12 = zero(T); s13 = zero(T); s23 = zero(T)

        for jj in 1:8
            vx = Cuboids_Vertices[jj,1,j]
            vy = Cuboids_Vertices[jj,2,j]
            vz = Cuboids_Vertices[jj,3,j]

            xbar = vx - x0
            ybar = vy - y0
            ζm   = vz - z0
            rm   = sqrt(xbar*xbar + ybar*ybar + ζm*ζm)

            sgn = T(Sσ[jj])

            # --- Displacements (your fFunc) ---
            d1 +=  pref_disp * sgn * fFunc(ybar, ζm, xbar, rm)
            d2 +=  pref_disp * sgn * fFunc(xbar, ζm, ybar, rm)
            d3 += -pref_disp * sgn * fFunc(xbar, ybar, ζm,  rm)

            # --- Stresses (guards to avoid NaNs/Inf) ---
            den11 = ybar*ζm
            den22 = xbar*ζm
            den33 = xbar*ybar

            den11 = ifelse(den11 == 0, (den11 >= 0 ? tiny : -tiny), den11)
            den22 = ifelse(den22 == 0, (den22 >= 0 ? tiny : -tiny), den22)
            den33 = ifelse(den33 == 0, (den33 >= 0 ? tiny : -tiny), den33)

            s11 +=  pref_sig * sgn * atan( (xbar*rm) / den11 )
            s22 +=  pref_sig * sgn * atan( (ybar*rm) / den22 )
            s33 +=  pref_sig * sgn * atan( (ζm  *rm) / den33 )

            a12 = rm - ζm;  a12 = ifelse(a12 > tiny, a12, tiny)            # log(rm - ζm)
            a13 = rm + ybar; a13 = ifelse(a13 > tiny, a13, tiny)            # log(rm + ybar)
            a23 = rm + xbar; a23 = ifelse(a23 > tiny, a23, tiny)            # log(rm + xbar)

            s12 += -pref_sig * sgn * log(a12)
            s13 +=  pref_sig * sgn * log(a13)
            s23 +=  pref_sig * sgn * log(a23)
        end

        # write back
        outDisp[1,i,j] = d1; outDisp[2,i,j] = d2; outDisp[3,i,j] = d3
        outSig[1,i,j]  = s11; outSig[2,i,j] = s22; outSig[3,i,j] = s33
        outSig[4,i,j]  = s12; outSig[5,i,j] = s13; outSig[6,i,j] = s23
    end
    return
end

# --------------------------
# Host-side driver: accepts Vector{Matrix{T}} for Cuboids_Vertices
# --------------------------
function gpu_porodispstress!(
    PoroDisp_GF_patch_cuboid::AbstractArray{T,3},   # (3, FaultCount, Cuboids_Count)
    PoroStress_GF_patch_cuboid::AbstractArray{T,3}, # (6, FaultCount, Cuboids_Count)
    FaultCenter::AbstractMatrix{T},                 # (FaultCount, 3)
    Cuboids_Vertices_vec::Vector{<:AbstractMatrix{T}},  # length = Cuboids_Count, each (8,3)
    Compressibility,                                 # scalar or Vector{T} (Cuboids_Count)
    ShearModulus::T
) where {T<:AbstractFloat}

    FC = size(FaultCenter,1)
    CC = length(Cuboids_Vertices_vec)
    @assert size(PoroDisp_GF_patch_cuboid)  == (3, FC, CC)
    @assert size(PoroStress_GF_patch_cuboid) == (6, FC, CC)

    # Pack to contiguous (8,3,CC) once (reuse this if you call repeatedly)
    CV3 = pack_vertices(Cuboids_Vertices_vec)

    # Normalize Compressibility to a vector
    Cm_vec = isa(Compressibility, AbstractVector) ? Compressibility :
             fill(T(Compressibility), CC)

    # move data to GPU
    d_FC = cu(FaultCenter)
    d_CV = cu(CV3)
    d_Cm = cu(Cm_vec)

    d_Disp = CUDA.CuArray{T}(undef, 3, FC, CC)
    d_Sig  = CUDA.CuArray{T}(undef, 6, FC, CC)

    threads = (16, 16)
    blocks  = (cld(FC, threads[1]), cld(CC, threads[2]))
    @cuda threads=threads blocks=blocks kernel_porodispstress!(d_Disp, d_Sig, d_FC, d_CV, d_Cm, ShearModulus)
    synchronize()

    # copy back
    PoroDisp_GF_patch_cuboid   .= Array(d_Disp)
    PoroStress_GF_patch_cuboid .= Array(d_Sig)
    return nothing
end
