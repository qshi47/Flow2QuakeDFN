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

include("../Functions_Kuvshinov_Cubes.jl")

function IF_FaultPatch_on_Reservoir(SingleFaultCenter, idx_nearest_cube, timeidx)
    if_faultpatch_on_reservoir = false
    porePressureChange_of_patch = 0.0

   if (SingleFaultCenter[3] >= 3950) & (SingleFaultCenter[3] <= 4050 )
        if_faultpatch_on_reservoir = true
        porePressureChange_of_patch = NonUniformPorePressureChange[timeidx, idx_nearest_cube] 
   end

   return if_faultpatch_on_reservoir, porePressureChange_of_patch

end

function Calculate_Nearest_CubeIdx_to_Fault(SingleFaultCenter)
    distances = sum((Reservoir_BLOCKS_Center_Pos .- SingleFaultCenter') .^ 2, dims=2)
    idx_nearest_cube =  argmin(distances)
    return idx_nearest_cube
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

function RotateTensor_to_Stress_on_Fault(RotationMat_FromFault, sigEff_all, FaultDipAngle, FaultStrikeAngle, FaultLLRR, SSorRN)    
    Stress_Fault = RotationMat_FromFault * sigEff_all * RotationMat_FromFault'    
    
    D_Stress_Normal = - Stress_Fault[3,3] # tension is positive (negative to make it compression)

    if SSorRN == 1 
        D_Stress_Shear = - FaultLLRR * Stress_Fault[1,3]  # Right Latteral become negative after rotation
    
    elseif SSorRN == 2

        if FaultDipAngle <= 90
                D_Stress_Shear = - FaultLLRR * Stress_Fault[2,3]  # Nomal orientation is negative when dip angle is <90
        else 
                D_Stress_Shear = FaultLLRR * Stress_Fault[2,3]  # Nomal orientation is positive when dip angle is <90
        end
    end

    return D_Stress_Normal, D_Stress_Shear

end

# Input Files 
LoadingInputFileName="Input_Discretized.jld2" 
LoadingReservoirCubesCoord = "Input_Example_Cubes.h5"
LoadingPorePressureChange = "Input_Example_PorePressureChange.csv"
OutputFileExternalStress = "Input_ExternalStressChange.jld2"

# Load Fault parameters
FaultCenter= load(LoadingInputFileName, "FaultCenter")
ShearModulus= load(LoadingInputFileName, "ShearModulus")
RockDensity= load(LoadingInputFileName, "RockDensity")
PoissonRatio= load(LoadingInputFileName, "PoissonRatio")
FaultLengthStrike= load(LoadingInputFileName, "FaultLengthStrike")
FaultLengthDip= load(LoadingInputFileName, "FaultLengthDip")
FaultStrikeAngle= load(LoadingInputFileName, "FaultStrikeAngle")
FaultDipAngle= load(LoadingInputFileName, "FaultDipAngle")
FaultLLRR= load(LoadingInputFileName, "FaultLLRR")

Fault_BulkIndex= load(LoadingInputFileName, "Fault_BulkIndex")
FaultLengthStrike_Bulk= load(LoadingInputFileName, "FaultLengthStrike_Bulk")
FaultLengthDip_Bulk= load(LoadingInputFileName, "FaultLengthDip_Bulk")
FaultCount= load(LoadingInputFileName, "FaultCount") 
LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount") 
Switch_StrikeSlip_or_ReverseNormal = load(LoadingInputFileName, "Switch_StrikeSlip_or_ReverseNormal")
println("Fault count number is:  ", FaultCount)

# Load Reservoir Positions
BlocksCount, Reservoir_BLOCKS_Center_Pos, Reservoir_BLOCKS_Vertices_Pos = Load_Reservoir_Cubes(LoadingReservoirCubesCoord)
println("Reservoir cube count number is:  ", BlocksCount)

# Load time array and pore pressure change 
# from CSV, skipping header
ExternalStress_TimeArray = readdlm(LoadingPorePressureChange, ',', skipstart=1)[:,1]
NonUniformPorePressureChange = readdlm(LoadingPorePressureChange, ',', skipstart=1)[:,2:end]

# Reservoir Elastic Properties
MyPossionRatio = 0.15
MyYoungModulus = 15e9
MyBulkModulus = MyYoungModulus/( 3*(1-2*MyPossionRatio) )
MyPwaveModulus = MyYoungModulus*(1-MyPossionRatio)/( (1+MyPossionRatio)*(1-2*MyPossionRatio) )
MyShearModulus = MyYoungModulus/(2*MyPossionRatio+2)
MyCompressibility = 1/MyPwaveModulus

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

PoroDisp_GF_patch_cube   = zeros(3, FaultCount, BlocksCount) # GF == Green Function
PoroStress_GF_patch_cube = zeros(6, FaultCount, BlocksCount)

println("---- Calculating Green Functions of Poro Elastic Stresses on Fault ----")
println("---- ", FaultCount, " * ", BlocksCount, " pairs  ----")

for i =  tqdm( 1:FaultCount, unit = "fault patch" )
    for j = 1:BlocksCount
        PoroDisp_GF_patch_cube[:,i,j], PoroStress_GF_patch_cube[:,i,j] = Calculate_PoroDispStress_SingleBlock_FullSpace_GF(FaultCenter[i,:], Reservoir_BLOCKS_Vertices_Pos[j], MyCompressibility, MyShearModulus, MyPossionRatio  )
    end
end


# Preparation work, avoid repeating computation
println("---- Preparing for Rotatating Stress Tensors ----")
Idx_Nearest_Cube_All = Vector{CartesianIndex}(undef, FaultCount) 
RotationMat_FromFault_All = zeros(FaultCount, 3, 3) 

for i = tqdm(1:FaultCount, unit="fault")
    Idx_Nearest_Cube_All[i] = Calculate_Nearest_CubeIdx_to_Fault(FaultCenter[i,:])
    RotationMat_FromFault_All[i, :, :] = Calculate_Rotation_Matrix_for_Fault(FaultDipAngle[i], FaultStrikeAngle[i])
end

# Rotate Poro Elastic Tensors to Stresses on Fault
println("---- Rotating Poro Elastic Tensors to Stress on Fault ----")
for (TimeIdx, Time) in tqdm( enumerate(ExternalStress_TimeArray), unit = " timestep" )
    for i =1:FaultCount
        # Convention: compression positive + Z downward positive (Left-hand system)
        PoroDisp_time_patch   = PoroDisp_GF_patch_cube[:,i,:]   * NonUniformPorePressureChange[TimeIdx]
        PoroStress_time_patch = PoroStress_GF_patch_cube[:,i,:] * NonUniformPorePressureChange[TimeIdx]

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
        if_faultpatch_contact_reservoir, PorePressureChange_of_patch = IF_FaultPatch_on_Reservoir(FaultCenter[i,:], Idx_Nearest_Cube_All[i], TimeIdx)
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
        RotateTensor_to_Stress_on_Fault(RotationMat_FromFault_All[i,:,:], SigEff_all,  
        FaultDipAngle[i], FaultStrikeAngle[i], FaultLLRR[i], Switch_StrikeSlip_or_ReverseNormal)

    end
end

# Write to JLD2 File for QuakeDFN solver
if Switch_StrikeSlip_or_ReverseNormal == 1 
    println("Stress calculated in strike-slip orientation")
elseif Switch_StrikeSlip_or_ReverseNormal == 2 
    println("Stress calculated in Reverse-Normal orientation")
end

save(OutputFileExternalStress,
"ExternalStress_TimeArray", ExternalStress_TimeArray, 
"ExternalStress_Normal", ExternalStress_Normal_Poro,
"ExternalStress_Shear", ExternalStress_Shear_Poro,
"Pressure",  PorePressure_Fault,
"NonUniformPorePressureChange", NonUniformPorePressureChange,
"PoroStress_11", PoroStress_11,
"PoroStress_22", PoroStress_22,
"PoroStress_33", PoroStress_33,
"PoroStress_12", PoroStress_12,
"PoroStress_13", PoroStress_13,
"PoroStress_23", PoroStress_23,
"PoroDisp_1",   PoroDisp_1,
"PoroDisp_2",   PoroDisp_2,
"PoroDisp_3",   PoroDisp_3)

println("JLD2 File Saved: ", OutputFileExternalStress  )