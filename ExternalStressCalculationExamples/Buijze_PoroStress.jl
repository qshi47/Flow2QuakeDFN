"""
1. Load the geomrtry containing an inclining fault 
2. Transport the cubes vertices coord data from Python to Julia throught H5DF file
3. Read Buijze19 and import the initial pore PorePressure_Poro within and without the cubes
4. Set a constant pore PorePressure_Poro extraction rate and record the reduced pore PorePressure_Poro
5. Calculate the normal, shear stress changes on fault 
"""

using DelimitedFiles
using JLD2
using LinearAlgebra
using Printf
# using PyPlot
using PyCall
using Statistics

# New, added by Qian
using ProgressBars
using SpecialFunctions: expinti
using HDF5

# for writing to VTK file for Paraview
using WriteVTK
using StaticArrays


println("---- Pre-loading ----")
include("../Results/Functions_Plot.jl")
include("../Qian/Utilities_GRG.jl")
include("../Functions_Buijze_InitialStress.jl")

println("---- Loading ----")
LoadingBuijzeGeomery = "Input_Buijze19_Geometry.jld2"
LoadingInputFileName="Input_Discretized.jld2" 

ReservoirCubesCoordFileName = "Input_Buijze19_Cubes.h5"

OutputFile="Input_ExternalStressChange.jld2"
OutputVTKFileName = "Buijze19_50_FaultStress"

# OutputFile_Reservoir_Shape = "Input_Reservoir_Shape.jld2"
# OutputFile_Reservoir_Center = "Input_Reservoir_Center.h5"
# OutputFile_Reservoir_Vertices = "Input_Reservoir_Vertices.h5"
 

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
FaultCount= load(LoadingInputFileName, "FaultCount") 
LoadingFaultCount= load(LoadingInputFileName, "LoadingFaultCount") 
Switch_StrikeSlip_or_ReverseNormal = load(LoadingInputFileName, "Switch_StrikeSlip_or_ReverseNormal")

# Load Reservoir cubes data
CubesFile = h5open( ReservoirCubesCoordFileName, "r" )
BLOCKS = read( CubesFile["cubes"] )
BLOCKS = permutedims(BLOCKS, (3, 2, 1))
close(CubesFile)

BlocksCount = size(BLOCKS)[1]
Reservoir_BLOCKS_Vertices_Pos = collect(eachslice(BLOCKS, dims=1))



# Load Reservoir Basic Geometry Info
Offset = load(LoadingBuijzeGeomery, "Offset")
Phi =   load(LoadingBuijzeGeomery, "Phi")   # deg


# function Time_with_Timestep(time_idx)
#     # return in second
#     return 0.001*1.2^i
#     # 100 ~ 0.002 yr
#     # 110 ~ 0.016 yr
#     # 120 ~ 0.1 yr
#     # 150 ~ 24 yr
# end

function Time_with_Timestep(time_idx)
    # return in second
    Year_in_Second = 3600*24*365
    return Year_in_Second*time_idx/10
    # 100 ~ 10 yr
    # 150 ~ 15 yr
end


TimeCount =  200 # 150, 6
ExternalStress_TimeArray = zeros(TimeCount)
for i=1:TimeCount
    ExternalStress_TimeArray[i]= Time_with_Timestep(i) # seconds 
end



MyPossionRatio = 0.15 # 0.15
MyYoungModulus = 15e9
MyBulkModulus = MyYoungModulus/( 3*(1-2*MyPossionRatio) )
MyPwaveMudulus = MyYoungModulus*(1-MyPossionRatio)/( (1+MyPossionRatio)*(1-2*MyPossionRatio) )
MyShearModulus = MyYoungModulus/(2*MyPossionRatio+2)



# >>>> Calculate Uniform Pore Pressure Change >>>>
println("---- Calculate Uniform Pore Pressure Change ----")
UniformPorePressureChange = zeros(TimeCount)
EndExternalPorePressure = -40e6
EndExternalTime = Time_with_Timestep(TimeCount)
PorePressureExtractionRate = EndExternalPorePressure/EndExternalTime
println("PP change rate is: ", 365*24*3600*PorePressureExtractionRate*1e-6, "MPa/year")

# Linear pore pressure change
for (time_idx, time) in enumerate(ExternalStress_TimeArray)
    UniformPorePressureChange[time_idx] = PorePressureExtractionRate*time
end

Delta_P = UniformPorePressureChange .* ones( TimeCount, BlocksCount )
Cm_Delta_P = 1/MyPwaveMudulus*Delta_P

function IF_FaultPatch_on_Reservoir(SingleFaulteCenter, time)
    If_faultpatch_on_reservoir = 0 # false
    PorePressureChange = 0.0

   if (SingleFaulteCenter[3] >= 2800) & (SingleFaulteCenter[3] <= 3000 + Offset)
        If_faultpatch_on_reservoir = 1
        timeidx = findfirst( ExternalStress_TimeArray .== time)
        PorePressureChange = Delta_P[timeidx, 1]
   end

   return If_faultpatch_on_reservoir, PorePressureChange

end


# >>>> All Stresses on Fault >>>>
TotalPlotFault=FaultCount
TimeArrayCount=length(ExternalStress_TimeArray)

ExternalStress_Normal = zeros(TimeArrayCount,TotalPlotFault)
ExternalStress_Shear  = zeros(TimeArrayCount,TotalPlotFault)
PorePressure          = zeros(TimeArrayCount,TotalPlotFault)

# >>>> Initial Stress >>>>
println("---- Calculate Initial Stress on Fault ----")
ExternalStress_Normal_Initial = zeros(TimeArrayCount,TotalPlotFault)
ExternalStress_Shear_Initial  = zeros(TimeArrayCount,TotalPlotFault)
PorePressure_Initial          = zeros(TimeArrayCount,TotalPlotFault)

for i = 1:TotalPlotFault
    fault_patch_depth = FaultCenter[i,3]
    intialstress_effnormal, intialstress_shear, initial_porepressure = CalculateInitialStress_on_Fault(Phi, fault_patch_depth)
    ExternalStress_Normal_Initial[:,i] .= intialstress_effnormal
    ExternalStress_Shear_Initial[:,i]  .= intialstress_shear
    PorePressure_Initial[:,i]          .= initial_porepressure
end



# >>>> Poro-elastic Stress >>>>
println("---- Calculating Poro Elastic Stresses on Fault ----")
ExternalStress_Normal_Poro = zeros(TimeArrayCount,TotalPlotFault)
ExternalStress_Shear_Poro  = zeros(TimeArrayCount,TotalPlotFault)
PorePressure_Poro     = zeros(TimeArrayCount,TotalPlotFault)

PoroStress_11 = zeros(TimeArrayCount,TotalPlotFault)
PoroStress_22 = zeros(TimeArrayCount,TotalPlotFault)
PoroStress_33 = zeros(TimeArrayCount,TotalPlotFault)
PoroStress_12 = zeros(TimeArrayCount,TotalPlotFault)
PoroStress_13 = zeros(TimeArrayCount,TotalPlotFault)
PoroStress_23 = zeros(TimeArrayCount,TotalPlotFault)

PoroDisp_1 = zeros(TimeArrayCount,TotalPlotFault)
PoroDisp_2 = zeros(TimeArrayCount,TotalPlotFault)
PoroDisp_3 = zeros(TimeArrayCount,TotalPlotFault)


function CalculatePoroStress_on_Fault(sigEff_all, FaultDipAngle, FaultStrikeAngle, FaultLLRR, SSorRN)
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

    return D_Stress_Normal, D_Stress_Shear

end


for (TimeIdx, Time) in tqdm( enumerate(ExternalStress_TimeArray), unit = " timestep" )
    for i =1:TotalPlotFault
        # All stress convention: compress positive
        PoroDisp_time_patch, PoroStress_time_patch = Calculate_DispStress_MultipleBlocks(FaultCenter[i,:], Reservoir_BLOCKS_Vertices_Pos, Cm_Delta_P[TimeIdx,:], MyShearModulus, MyPossionRatio  )
        
        PoroDisp_1[TimeIdx,i]   = PoroDisp_time_patch[1] 
        PoroDisp_2[TimeIdx,i]   = PoroDisp_time_patch[2] 
        PoroDisp_3[TimeIdx,i]   = PoroDisp_time_patch[3] 
        
        PoroStress_11[TimeIdx,i] = PoroStress_time_patch[1] 
        PoroStress_22[TimeIdx,i] = PoroStress_time_patch[2] 
        PoroStress_33[TimeIdx,i] = PoroStress_time_patch[3] 
        PoroStress_12[TimeIdx,i] = PoroStress_time_patch[4] 
        PoroStress_13[TimeIdx,i] = PoroStress_time_patch[5] 
        PoroStress_23[TimeIdx,i] = PoroStress_time_patch[6] 
        
        # Judge if fault contact the reservoir 
        if_faultpatch_contact_reservoir, PorePressureChange = IF_FaultPatch_on_Reservoir(FaultCenter[i,:], Time)
        PorePressure_Poro[TimeIdx,i] = PorePressureChange * if_faultpatch_contact_reservoir

        # Calculate Effective Stress (Tensional stress is positive)
        SigEff_all = zeros(3,3)
        SigEff_all = 
        [-PoroStress_11[TimeIdx,i]     -PoroStress_12[TimeIdx,i]   -PoroStress_13[TimeIdx,i]
         -PoroStress_12[TimeIdx,i]     -PoroStress_22[TimeIdx,i]   -PoroStress_23[TimeIdx,i]
         -PoroStress_13[TimeIdx,i]     -PoroStress_23[TimeIdx,i]   -PoroStress_33[TimeIdx,i] ] +
        [PorePressure_Poro[TimeIdx,i]    0.0            0.0
        0.0            PorePressure_Poro[TimeIdx,i]     0.0
        0.0            0.0            PorePressure_Poro[TimeIdx,i]] 

        ExternalStress_Normal_Poro[TimeIdx,i], ExternalStress_Shear_Poro[TimeIdx,i] = CalculatePoroStress_on_Fault(SigEff_all,  
        FaultDipAngle[i], FaultStrikeAngle[i], FaultLLRR[i], Switch_StrikeSlip_or_ReverseNormal)

    end
end


ExternalStress_Normal = ExternalStress_Normal_Poro
ExternalStress_Shear  = ExternalStress_Shear_Poro
PorePressure          = PorePressure_Poro

println("---- Saving File: ", OutputFile, " ----" )
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

# save_vector_of_matrices(OutputFile_Reservoir_Vertices, Reservoir_BLOCKS_Vertices_Pos )
# save_vector_of_vectors(OutputFile_Reservoir_Center, Reservoir_BLOCKS_Origin_Pos )

if Switch_StrikeSlip_or_ReverseNormal == 1 
        println("Stress calculated in strike-slip orientation")
elseif Switch_StrikeSlip_or_ReverseNormal == 2 
        println("Stress calculated in Reverse-Normal orientation")
end

save(OutputFile,
"ExternalStress_TimeArray", ExternalStress_TimeArray, 
"ExternalStress_Normal", ExternalStress_Normal,
"ExternalStress_Normal_Initial", ExternalStress_Normal_Initial,  
"ExternalStress_Normal_Poro", ExternalStress_Normal_Poro, 
"ExternalStress_Shear", ExternalStress_Shear,
"ExternalStress_Shear_Initial", ExternalStress_Shear_Initial,
"ExternalStress_Shear_Poro", ExternalStress_Shear_Poro,
"Pressure",  PorePressure,
"PorePressure", PorePressure,
"PorePressure_Initial", PorePressure_Initial,
"PorePressure_Poro", PorePressure_Poro,
"Delta_P", Delta_P,
"Cm_Delta_P", Cm_Delta_P,
"PoroStress_11", PoroStress_11,
"PoroStress_22", PoroStress_22,
"PoroStress_33", PoroStress_33,
"PoroStress_12", PoroStress_12,
"PoroStress_13", PoroStress_13,
"PoroStress_23", PoroStress_23,
"PoroDisp_1",   PoroDisp_1,
"PoroDisp_2",   PoroDisp_2,
"PoroDisp_3",   PoroDisp_3)

# save(OutputFile_Reservoir_Shape, 
# "Reservoir_BLOCKS_Origin_Pos", Reservoir_BLOCKS_Origin_Pos,
# "Reservoir_BLOCKS_Vertices_Pos", Reservoir_BLOCKS_Vertices_Pos)




# >>>> Write to VTK File >>>>
println("---- Saving File: ", OutputVTKFileName, " ----" )
coordinates = [SVector(row...) for row in eachrow(FaultCenter)]
vtk = vtk_grid(OutputVTKFileName, coordinates )

# Add scalar data for each point
vtk["ExternalStress_Shear_Poro"] = ExternalStress_Shear_Poro
vtk["ExternalStress_Shear"] = ExternalStress_Shear
vtk["ExternalStress_Normal_Poro"] = ExternalStress_Normal_Poro
vtk["ExternalStress_Normal"] = ExternalStress_Normal
vtk["PorePressure_Poro"] = PorePressure_Poro
vtk["PorePressure"] = PorePressure
vtk["PoroStress_11"] = PoroStress_11
vtk["PoroStress_33"] = PoroStress_33

# Save the file
vtk_save(vtk)






println("[DONE!]")