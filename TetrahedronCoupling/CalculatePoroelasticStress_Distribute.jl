using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
using Statistics


using ProgressBars
using HDF5

# for writing to VTK file for Paraview
using WriteVTK
using StaticArrays

using Distributed
addprocs(16)  

@everywhere begin
    include("./Functions_Kuvshinov_Tetrahedron_Faces.jl")
end
using SharedArrays


function If_FaultPatch_on_Reservoir(Criterion_for_Fault_within_Reservoir::Function, SingleFaultCenter, PorePressureChange, timeidx, idx_nearest_tetrahedron )
    if_faultpatch_on_reservoir = false
    porePressureChange_of_patch = 0.0

   if Criterion_for_Fault_within_Reservoir(SingleFaultCenter) == true
        if_faultpatch_on_reservoir = true
        porePressureChange_of_patch = PorePressureChange[timeidx, idx_nearest_tetrahedron] 
   end

   return if_faultpatch_on_reservoir, porePressureChange_of_patch

end

function Calculate_Nearest_CubeIdx_to_Fault(Tetrahedrons_Center, SingleFaultCenter)
    distances = sum((Tetrahedrons_Center .- SingleFaultCenter') .^ 2, dims=2)
    idx_nearest_tetrahedron =  argmin(distances)
    return idx_nearest_tetrahedron
end

function Calculate_Rotation_Matrix_for_Fault( FaultDipAngle, FaultStrikeAngle )
    # Rotate the fault center stress to reference frame to read shear and normal
    RotationMat_FromFault_Strike=
    [cosd(-FaultStrikeAngle) -sind(-FaultStrikeAngle)  0
    sind(-FaultStrikeAngle) cosd(-FaultStrikeAngle) 0
    0  0  1];

    RotationMat_FromFault_Dip=
    [1 0 0
    0 cosd(-FaultDipAngle) -sind(-FaultDipAngle)
    0 sind(-FaultDipAngle) cosd(-FaultDipAngle)]
    
    RotationMat_FromFault = RotationMat_FromFault_Dip*RotationMat_FromFault_Strike;
    
    return RotationMat_FromFault
end

function Rotate_Tensor_to_Stress_on_Fault(RotationMat_FromFault, sigEff_all, FaultRakeAngle)    
    Stress_Fault = RotationMat_FromFault * sigEff_all * RotationMat_FromFault'    
    
    D_Stress_Normal = - Stress_Fault[3,3] # tension is positive (negative to make it compression)
    D_Stress_Shear =  cosd(FaultRakeAngle) * Stress_Fault[1,3] + sind(FaultRakeAngle) * Stress_Fault[2,3]

    return D_Stress_Normal, D_Stress_Shear

end


# Input Files 
Input_Fault_File="Input_Discretized.jld2" 
Input_Tetrahedrons_File = "TetrahedronCoupling/Input_Tetrahedron.h5"
Input_PorePressure_File = "TetrahedronCoupling/Input_PorePressure.txt"
Output_ExternalStress_File = "Input_ExternalStressChange.jld2"

Input_Tetrahedron_Faces_File = Input_Tetrahedrons_File

# Criterion for Fault within Reservoir
function Criterion_for_Fault_within_Reservoir(SingleFaultCenter)
    return (SingleFaultCenter[3] >= 2800) & (SingleFaultCenter[3] <= 3000 )
end

# function main(Criterion_for_Fault_within_Reservoir::Function, Input_Fault_File, Input_Tetrahedron_Faces_File, Input_PorePressure_File, Output_ExternalStress_File)

# Load Fault parameters
println("---- Loading Fault and Tetrahedron Faces  ----")
FaultCenter= load(Input_Fault_File, "FaultCenter")
ShearModulus= load(Input_Fault_File, "ShearModulus")
PoissonRatio= load(Input_Fault_File, "PoissonRatio")
FaultStrikeAngle= load(Input_Fault_File, "FaultStrikeAngle")
FaultDipAngle= load(Input_Fault_File, "FaultDipAngle")
FaultRakeAngle= load(Input_Fault_File, "FaultRakeAngle")
FaultCount= load(Input_Fault_File, "FaultCount") 
Switch_StrikeSlip_or_ReverseNormal = load(Input_Fault_File, "Switch_StrikeSlip_or_ReverseNormal")
println("Fault count is:  ", FaultCount)

# Load Reservoir Tetrahedrons
TetrahedronsFile = h5open( Input_Tetrahedron_Faces_File, "r" )
Tetrahedrons = read( TetrahedronsFile["vertices_coord"] )
Tetrahedrons = permutedims(Tetrahedrons, (4, 3, 2, 1))
close(TetrahedronsFile)

TetrahedronCount = size(Tetrahedrons)[1]
Tetrahedrons_faces = collect(eachslice(Tetrahedrons, dims=1))
Tetrahedrons_Center = zeros(TetrahedronCount, 3)
for i = 1:TetrahedronCount
    Tetrahedrons_Center[i, :] = mean(Tetrahedrons_faces[i], dims= [1,2])[:]
end
println("Tetrahedron count is:  ", TetrahedronCount)

# Load external time array and pore pressure change 
# from CSV, skipping header
ExternalStress_TimeArray = readdlm(Input_PorePressure_File, ',', skipstart=1)[:,1]
PorePressureChange = readdlm(Input_PorePressure_File, ',', skipstart=1)[:,2:end]

# Reservoir Elastic Properties
PwaveModulus = 2 * ShearModulus * (1 - PoissonRatio) / (1 - 2 * PoissonRatio)
Compressibility = 1/PwaveModulus

# External Poro-elastic Stress Change on Fault 
# (initial stresses are assigned in QuickParameterChange.jl) 
TimeArrayCount=length(ExternalStress_TimeArray)

ExternalStress_Normal_Poro = zeros(TimeArrayCount,FaultCount)
ExternalStress_Shear_Poro  = zeros(TimeArrayCount,FaultCount)
PorePressure_Fault     = zeros(TimeArrayCount,FaultCount)

PoroStress_11 = zeros(TimeArrayCount,FaultCount)
PoroStress_22 = zeros(TimeArrayCount,FaultCount)
PoroStress_33 = zeros(TimeArrayCount,FaultCount)
PoroStress_12 = zeros(TimeArrayCount,FaultCount)
PoroStress_13 = zeros(TimeArrayCount,FaultCount)
PoroStress_23 = zeros(TimeArrayCount,FaultCount)

PoroDisp_1 = zeros(TimeArrayCount,FaultCount)
PoroDisp_2 = zeros(TimeArrayCount,FaultCount)
PoroDisp_3 = zeros(TimeArrayCount,FaultCount)

PoroDisp_GF_patch_tetrahedron   = SharedArray{Float64}(3, FaultCount, TetrahedronCount)
PoroStress_GF_patch_tetrahedron = SharedArray{Float64}(6, FaultCount, TetrahedronCount)


println("---- Calculating Green Functions of Poro Elastic Stresses on Fault ----")
println("---- ($(FaultCount) * $(TetrahedronCount) pairs) on $(nprocs()-1) workers ----")


@everywhere function _compute_stress!(A, FaultCenter, Tetrahedrons_faces, PwaveModulus, ShearModulus,
                                    i_start::Int, i_stop::Int, TetrahedronCount::Int)
    @inbounds @views for i in i_start:i_stop
        for j in 1:TetrahedronCount
            A[:, i, j] .= Calculate_PoroStress_SingleTetrahedron_FullSpace_GF(
                FaultCenter[i, :], Tetrahedrons_faces[j], PwaveModulus, ShearModulus)
        end
    end
    return nothing
end

# split the outer loop (i) across workers
chunks = collect(Iterators.partition(1:FaultCount, ceil(Int, FaultCount / (nprocs()-1))))

@sync begin
    for (w, ch) in zip(workers(), chunks)
        @async remotecall_wait( _compute_stress!, w,
            PoroStress_GF_patch_tetrahedron, FaultCenter, Tetrahedrons_faces,
            PwaveModulus, ShearModulus, first(ch), last(ch), TetrahedronCount)
    end
end


# Preparation work, avoid repeating computation
println("---- Preparing for Rotatating Stress Tensors ----")
Idx_Nearest_Cube_All = Vector{CartesianIndex}(undef, FaultCount) 
RotationMat_FromFault_All = zeros(FaultCount, 3, 3) 

for i = tqdm(1:FaultCount, unit="fault")
    Idx_Nearest_Cube_All[i] = Calculate_Nearest_CubeIdx_to_Fault(Tetrahedrons_Center, FaultCenter[i,:])
    RotationMat_FromFault_All[i, :, :] = Calculate_Rotation_Matrix_for_Fault(FaultDipAngle[i], FaultStrikeAngle[i])
end

# Rotate Poro Elastic Tensors to Stresses on Fault
println("---- Rotating Poro Elastic Tensors to Stress on Fault ----")
for (TimeIdx, Time) in tqdm( enumerate(ExternalStress_TimeArray), unit = " timestep" )
    for i =1:FaultCount
        # Convention: compression positive + Z downward positive (Left-hand system)
        PoroDisp_time_patch   = PoroDisp_GF_patch_tetrahedron[:,i,:]   * PorePressureChange[TimeIdx, :]
        PoroStress_time_patch = PoroStress_GF_patch_tetrahedron[:,i,:] * PorePressureChange[TimeIdx, :]

        PoroDisp_1[TimeIdx,i]   = PoroDisp_time_patch[1] 
        PoroDisp_2[TimeIdx,i]   = PoroDisp_time_patch[2] 
        PoroDisp_3[TimeIdx,i]   = PoroDisp_time_patch[3] 
        
        PoroStress_11[TimeIdx,i] = PoroStress_time_patch[1] 
        PoroStress_22[TimeIdx,i] = PoroStress_time_patch[2] 
        PoroStress_33[TimeIdx,i] = PoroStress_time_patch[3] 
        PoroStress_12[TimeIdx,i] = PoroStress_time_patch[4] 
        PoroStress_13[TimeIdx,i] = PoroStress_time_patch[5] 
        PoroStress_23[TimeIdx,i] = PoroStress_time_patch[6] 
        
        # Check if fault contact the reservoir 
        if_faultpatch_contact_reservoir, PorePressureChange_of_patch = If_FaultPatch_on_Reservoir(Criterion_for_Fault_within_Reservoir, FaultCenter[i,:], PorePressureChange, TimeIdx, Idx_Nearest_Cube_All[i])
        if if_faultpatch_contact_reservoir == true
            PorePressure_Fault[TimeIdx,i] = PorePressureChange_of_patch
        else
            PorePressure_Fault[TimeIdx,i] = 0.0
        end

        # Calculate Effective Stress (Change to Tensional stress positive)
        SigEff_all = zeros(3,3)
        SigEff_all = 
        [-PoroStress_11[TimeIdx,i]     -PoroStress_12[TimeIdx,i]   -PoroStress_13[TimeIdx,i]
        -PoroStress_12[TimeIdx,i]     -PoroStress_22[TimeIdx,i]   -PoroStress_23[TimeIdx,i]
        -PoroStress_13[TimeIdx,i]     -PoroStress_23[TimeIdx,i]   -PoroStress_33[TimeIdx,i] ] +
        [PorePressure_Fault[TimeIdx,i]    0.0            0.0
        0.0            PorePressure_Fault[TimeIdx,i]     0.0
        0.0            0.0            PorePressure_Fault[TimeIdx,i]] 

        ExternalStress_Normal_Poro[TimeIdx,i], ExternalStress_Shear_Poro[TimeIdx,i] = 
        Rotate_Tensor_to_Stress_on_Fault(RotationMat_FromFault_All[i,:,:], SigEff_all, FaultRakeAngle[i])

    end
end

# Write to JLD2 File for QuakeDFN solver
if Switch_StrikeSlip_or_ReverseNormal == 1 
    println("---- Stress calculated in strike-slip orientation ----")
elseif Switch_StrikeSlip_or_ReverseNormal == 2 
    println("---- Stress calculated in Reverse-Normal orientation ----")
end

save(Output_ExternalStress_File,
"ExternalStress_TimeArray", ExternalStress_TimeArray, 
"ExternalStress_Normal", ExternalStress_Normal_Poro,
"ExternalStress_Shear", ExternalStress_Shear_Poro,
"Pressure",  PorePressure_Fault,
"PorePressureChange", PorePressureChange,
"PoroStress_11", PoroStress_11,
"PoroStress_22", PoroStress_22,
"PoroStress_33", PoroStress_33,
"PoroStress_12", PoroStress_12,
"PoroStress_13", PoroStress_13,
"PoroStress_23", PoroStress_23,
"PoroDisp_1",   PoroDisp_1,
"PoroDisp_2",   PoroDisp_2,
"PoroDisp_3",   PoroDisp_3)

println("---- JLD2 File Saved: ", Output_ExternalStress_File, " ----")

# end





# main(Criterion_for_Fault_within_Reservoir, Input_Fault_File, Input_Tetrahedrons_File, Input_PorePressure_File, Output_ExternalStress_File)