using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
# using PyPlot
using PyCall
using Statistics

# added by Qian
using ProgressBars
using SpecialFunctions: expinti
using HDF5

# pygui(false)

include("../Results/Functions_Plot.jl")
include("../Qian/Utilities_GRG.jl")


LoadingInputFileName="Input_Discretized.jld2" 

OutputFile="Input_ExternalStressChange.jld2"
OutputFile_Reservoir_Shape = "Input_Reservoir_Shape.jld2"
OutputFile_Reservoir_Center = "Input_Reservoir_Center.h5"
OutputFile_Reservoir_Vertices = "Input_Reservoir_Vertices.h5"


############################### Load Input Files ###############################
######++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++######
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
FaultCount= load(LoadingInputFileName, "FaultCount") # 400
LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount") # 400
Switch_StrikeSlip_or_ReverseNormal = load(LoadingInputFileName, "Switch_StrikeSlip_or_ReverseNormal")




TimeCount = 150

ExternalStress_TimeArray=zeros(TimeCount)
for i=1:TimeCount
    ExternalStress_TimeArray[i]=0.001*1.2^i # seconds 
end

MyPossionRatio = 0.25
MyYoungModulus = 1e9
MyBulkModulus = MyYoungModulus/( 3*(1-2*MyPossionRatio) )
MyPwaveMudulus = MyYoungModulus*(1-MyPossionRatio)/( (1+MyPossionRatio)*(1-2*MyPossionRatio) )
MyShearModulus = MyYoungModulus/(2*MyPossionRatio+2)


CubeLengthX = 200
CubeLengthY = 200
CubeThickness = 100
Reservoir_BLOCKS_Origin_Pos = Vector{Vector{Float64}}()
Reservoir_BLOCKS_Vertices_Pos = Vector{Matrix{Float64}}()


# >>>>> Set Reservoir Geometry, Compressibility and Pressure History >>>>>>
# Geometry1: Single Cube
"""
CubeOrigin = [ CubeLengthX/2, 0.0, 2000.0 ] # offset = 200, [ offset + CubeLengthX/2,0,2000]
CUBE = Get_Cube_8Vertices_Pos(CubeOrigin, CubeLengthX, CubeLengthY,CubeThickness)
Reservoir_BLOCKS_Origin_Pos = [CubeOrigin]
Reservoir_BLOCKS_Vertices_Pos = [CUBE] # previously called BLOCKS

MyConstantPressure = 1.0
Delta_P = MyConstantPressure*ones( TimeCount, length(Reservoir_BLOCKS_Origin_Pos) )
Cm_Delta_P = 1/MyPwaveMudulus*Delta_P


function IF_FaultPatch_on_Reservoir(SingleFaulteCenter time)
    If_faultpatch_on_reservoir = 0 # false
    PorePressure = 0.0

   if (SingleFaulteCenter[1]==0) & (-CubeLengthX/2 <= SingleFaulteCenter[2] <= +CubeLengthX/2) &
      (2000-CubeThickness/2 <= SingleFaulteCenter[3] <= 2000+CubeThickness/2)
        If_faultpatch_on_reservoir = 1
        PorePressure = MyConstantPressure
   end

   return If_faultpatch_on_reservoir, PorePressure

end
"""

# Geometry2: Plate 
# """
Reservoir_Center = [0, 0, 2000 - CubeThickness/2]
Num_for_EachSide = 8 # Original Num_for_EachSide*CubeLengthX  == 1600
Reservoir_BLOCKS_Origin_Pos = BuildReservoirGeometry_even( Reservoir_Center, CubeLengthX, CubeThickness, Num_for_EachSide )

for Cube_center in Reservoir_BLOCKS_Origin_Pos
        push!( Reservoir_BLOCKS_Vertices_Pos, Get_Cube_8Vertices_Pos( Cube_center, CubeLengthX, CubeLengthY, CubeThickness ) )
end 


MyDiffusivity = 10 # m/second 
MyP_star = 4e6 # Pa
InjectionOrigin = [Reservoir_Center[1], Reservoir_Center[2], 2000 ]

Delta_P = zeros( TimeCount, length(Reservoir_BLOCKS_Origin_Pos) )
DISP_to_injection_origin = zeros( length(Reservoir_BLOCKS_Origin_Pos) )


println("---- Calculate Pressure Change ----")
for (TimeIdx, Time) in tqdm( enumerate(ExternalStress_TimeArray), unit = " timestep"  )
    for (idxCube, Cube_center) in enumerate( Reservoir_BLOCKS_Origin_Pos )
        disp_to_injection_origin = norm(Cube_center - InjectionOrigin)
        DISP_to_injection_origin[idxCube] = disp_to_injection_origin
        Delta_P[TimeIdx, idxCube] = - expinti(  - disp_to_injection_origin^2 / (4*MyDiffusivity*Time)  ) * MyP_star
    end
end

Cm_Delta_P = 1/MyPwaveMudulus*Delta_P

function IF_FaultPatch_on_Reservoir(SingleFaulteCenter, time)
    If_faultpatch_on_reservoir = 0 # false
    PorePressure = 0.0

   if (SingleFaulteCenter[1]==0) & 
      (- CubeLengthX*Num_for_EachSide/2 <= SingleFaulteCenter[2] <= + CubeLengthX*Num_for_EachSide/2) &
      (2000 - CubeThickness/2 <= SingleFaulteCenter[3] <= 2000 + CubeThickness/2)
        
        If_faultpatch_on_reservoir = 1
        
        disp_to_injection_origin = norm( SingleFaulteCenter - InjectionOrigin )
        PorePressure = - expinti(  - disp_to_injection_origin^2 / (4*MyDiffusivity * time)  ) * MyP_star
   end

   return If_faultpatch_on_reservoir, PorePressure

end
# """

# <<<<<<<<<<<<<<<<<<<<<



TotalPlotFault=FaultCount
TimeArrayCount=length(ExternalStress_TimeArray)
ExternalStress_Normal = zeros(TimeArrayCount,TotalPlotFault)
ExternalStress_Shear = zeros(TimeArrayCount,TotalPlotFault)
Pressure             = zeros(TimeArrayCount,TotalPlotFault)

SIGMA_11 = zeros(TimeArrayCount,TotalPlotFault)
SIGMA_22 = zeros(TimeArrayCount,TotalPlotFault)
SIGMA_33 = zeros(TimeArrayCount,TotalPlotFault)
SIGMA_12 = zeros(TimeArrayCount,TotalPlotFault)
SIGMA_13 = zeros(TimeArrayCount,TotalPlotFault)
SIGMA_23 = zeros(TimeArrayCount,TotalPlotFault)

Disp_1 = zeros(TimeArrayCount,TotalPlotFault)
Disp_2 = zeros(TimeArrayCount,TotalPlotFault)
Disp_3 = zeros(TimeArrayCount,TotalPlotFault)



function CalculateStress_MultipleBlocks(time, SingleFaulteCenter, Blocks, cm_delta_P, 
                                        FaultDipAngle, FaultStrikeAngle, FaultLLRR, 
                                        SSorRN, ShearModulus, ν )
    Disp, Sigma = Calculate_DispStress_MultipleBlocks(SingleFaulteCenter, Blocks, cm_delta_P, ShearModulus, ν)
    
    sig_11 = Sigma[1]
    sig_22 = Sigma[2]
    sig_33 = Sigma[3]
    sig_12 = Sigma[4]
    sig_13 = Sigma[5]
    sig_23 = Sigma[6]

    if_faultpatch_on_reservoir, PorePressure = IF_FaultPatch_on_Reservoir(SingleFaulteCenter, time)
    
    # Calculate Effective Stress (Tensional stress is positive)
    sigEff_all=[-sig_11 -sig_12 -sig_13
                -sig_12 -sig_22 -sig_23
                -sig_13 -sig_23 -sig_33] + 
                [PorePressure   0.0            0.0
                 0.0            PorePressure   0.0
                 0.0            0.0            PorePressure] * if_faultpatch_on_reservoir

    # Rotate the fault center stress to reference frame to read shear and normal
    RotationMat_FromFault_Strike=
    [cosd(-FaultStrikeAngle) -sind(-FaultStrikeAngle)  0
    sind(-FaultStrikeAngle) cosd(-FaultStrikeAngle) 0
    0  0  1];

    RotationMat_FromFault_Dip=
    [1 0 0
    0 cosd(-FaultDipAngle) -sind(-FaultDipAngle)
    0 sind(-FaultDipAngle) cosd(-FaultDipAngle)]
    
    RotationMat_FromFault_All = RotationMat_FromFault_Dip*RotationMat_FromFault_Strike;
    
    Stress_Fault = RotationMat_FromFault_All * sigEff_all * RotationMat_FromFault_All'    
    
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

    return D_Stress_Normal, D_Stress_Shear, PorePressure, [sig_11, sig_22, sig_33, sig_12, sig_13, sig_23], Disp
end


println("---- Calculate Stresses on Fault ----")
for (TimeIdx, Time) in tqdm( enumerate(ExternalStress_TimeArray), unit = " timestep" )
    for i =1:TotalPlotFault
        ExternalStress_Normal[TimeIdx,i], ExternalStress_Shear[TimeIdx,i], Pressure[TimeIdx,i], SIGMA_time_patch, Disp_time_patch  = 
            CalculateStress_MultipleBlocks( Time, FaultCenter[i,:], Reservoir_BLOCKS_Vertices_Pos, Cm_Delta_P[TimeIdx,:],
                                            FaultDipAngle[i], FaultStrikeAngle[i], FaultLLRR[i], 
                                            Switch_StrikeSlip_or_ReverseNormal, MyShearModulus, MyPossionRatio )
        
        SIGMA_11[TimeIdx,i] = SIGMA_time_patch[1] 
        SIGMA_22[TimeIdx,i] = SIGMA_time_patch[2] 
        SIGMA_33[TimeIdx,i] = SIGMA_time_patch[3] 
        SIGMA_12[TimeIdx,i] = SIGMA_time_patch[4] 
        SIGMA_13[TimeIdx,i] = SIGMA_time_patch[5] 
        SIGMA_23[TimeIdx,i] = SIGMA_time_patch[6] 
        Disp_1[TimeIdx,i]   = Disp_time_patch[1] 
        Disp_2[TimeIdx,i]   = Disp_time_patch[2] 
        Disp_3[TimeIdx,i]   = Disp_time_patch[3] 

    end
end




function save_vector_of_matrices(filename::String, data::Vector{Matrix{Float64}})
    # Open an HDF5 file for writing
    h5open(filename, "w") do file
        # Create a group to store the vector of matrices
        g = create_group(file, "matrices")
        
        # Save each matrix in the vector
        for (i, matrix) in enumerate(data)
            # Create a dataset for each matrix
            g[  @sprintf("matrix_%04d", i-1) ] = matrix
        end
    end
end

function save_vector_of_vectors(filename::String, data::Vector{Vector{Float64}})
    # Open an HDF5 file for writing
    h5open(filename, "w") do file
        # Create a group to store the vector of matrices
        g = create_group(file, "vectors")
        
        # Save each matrix in the vector
        for (i, vector) in enumerate(data)
            # Create a dataset for each matrix
            g[ @sprintf("matrix_%04d", i-1)] = vector
        end
    end
end


save_vector_of_matrices(OutputFile_Reservoir_Vertices, Reservoir_BLOCKS_Vertices_Pos )
save_vector_of_vectors(OutputFile_Reservoir_Center, Reservoir_BLOCKS_Origin_Pos )


println("File Saved: ", OutputFile)
if Switch_StrikeSlip_or_ReverseNormal == 1 
        println("Stress calculated in strike-slip orientation")
elseif Switch_StrikeSlip_or_ReverseNormal == 2 
        println("Stress calculated in Reverse-Normal orientation")
end

save(OutputFile,
"ExternalStress_TimeArray", ExternalStress_TimeArray, 
"ExternalStress_Normal", ExternalStress_Normal, 
"ExternalStress_Shear", ExternalStress_Shear,
"Pressure", Pressure,
"Delta_P", Delta_P,
"Cm_Delta_P", Cm_Delta_P,
"SIGMA_11", SIGMA_11,
"SIGMA_22", SIGMA_22,
"SIGMA_33", SIGMA_33,
"SIGMA_12", SIGMA_12,
"SIGMA_13", SIGMA_13,
"SIGMA_23", SIGMA_23,
"Disp_1",   Disp_1,
"Disp_2",   Disp_2,
"Disp_3",   Disp_3)

save(OutputFile_Reservoir_Shape, 
"Reservoir_BLOCKS_Origin_Pos", Reservoir_BLOCKS_Origin_Pos,
"Reservoir_BLOCKS_Vertices_Pos", Reservoir_BLOCKS_Vertices_Pos)


println("[DONE!]")