using DelimitedFiles
using JLD2
using HDF5
using LinearAlgebra
using ProgressBars
using SpecialFunctions: expinti
using CSV
using DataFrames

include("../Functions_Kuvshinov_Cubes.jl")

function Time_with_Timestep(time_idx)
    """
    # return in second
    # 100 ~ 10 yr
    # 200 ~ 20 yr
    """
    Year_in_Second = 3600*24*365
    return Year_in_Second*time_idx/10
end

function Porepressure_Injection(InjectionOrigin, P_star, Diffusivity, ExternalStress_TimeArray, Reservoir_BLOCKS_Center_Pos)
    time_count, blocks_count = size(ExternalStress_TimeArray, 1), size(Reservoir_BLOCKS_Center_Pos, 1)
    NonUniformPorePressureChange = zeros( time_count, blocks_count )
    Distances_to_injection_origin = zeros( blocks_count )
    
    for idxCube = 1:blocks_count
        Cube_center = Reservoir_BLOCKS_Center_Pos[idxCube,:]
        Distances_to_injection_origin[idxCube] = norm(Cube_center - InjectionOrigin)
    end

    for (idxTime, Time) in tqdm( enumerate(ExternalStress_TimeArray), unit = " timestep"  )
        for idxCube = 1:blocks_count
            disp_to_injection_origin = Distances_to_injection_origin[idxCube]
            NonUniformPorePressureChange[idxTime, idxCube] = - expinti(  - disp_to_injection_origin^2 / (4*Diffusivity*Time)  ) * P_star
        end
    end
    return NonUniformPorePressureChange
end 

function Save_time_pressure_to_csv(filename::String, time_array::AbstractVector, pressure::AbstractMatrix)
    output_matrix = hcat(time_array, pressure)
    blocksCount = size(pressure, 2)

    header = ["time(sec)"]
    for i in 1:blocksCount
        push!(header, "cube$(i)")
    end
    df = DataFrame(output_matrix, Symbol.(header))
    CSV.write(filename, df)
end


function main()
    LoadingReservoirCubesCoord = "Input_Example_Cubes.h5"
    PorePressureChangeOutputFile = "Input_Example_PorePressureChange.csv"

    # Load Reservoir Centers
    BlocksCount, Reservoir_BLOCKS_Center_Pos, Reservoir_BLOCKS_Vertices_Pos = Load_Reservoir_Cubes(LoadingReservoirCubesCoord)
    println("Reservoir cube count number is:  ", BlocksCount)

    # Time for External Stresses
    TimeCount =  200
    ExternalStress_TimeArray = zeros(TimeCount)
    for i = 0:TimeCount-1
        ExternalStress_TimeArray[i+1]= Time_with_Timestep(i) # seconds 
    end
    println("Maximum time is ", ExternalStress_TimeArray[end]/(3600*24*365), " years")


    # Non-Uniform Pore Pressure Change
    println("---- Calculate Non-Uniform Pore Pressure Change ----")
    MyDiffusivity = 1e-2 # m/second  # 1e-2: 10yr ~ 2000m
    MyP_star = 4e6 # Pa
    InjectionOrigin = [0, 0, 4000]

    NonUniformPorePressureChange = Porepressure_Injection(InjectionOrigin, MyP_star, MyDiffusivity, ExternalStress_TimeArray, Reservoir_BLOCKS_Center_Pos)


    # Save Non-Uniform Pore Pressure Change to a csv file
    Save_time_pressure_to_csv(PorePressureChangeOutputFile, ExternalStress_TimeArray, NonUniformPorePressureChange)

    println("CSV OutputFile Saved: ", PorePressureChangeOutputFile)
end

main()

