using DelimitedFiles
using CSV
using DataFrames


include("../Functions_Kuvshinov_Cuboids.jl")



function Save_time_pressure_to_txt(filename::String, time_array::AbstractVector, pressure::AbstractMatrix)
    open(filename, "w") do io
        write(io, "time(sec), ")
        writedlm(io, time_array', ", ")
        for i in 1:size(pressure, 2)
            write(io, "cuboid$(i), ")
            writedlm(io, pressure[:, i]', ", ")
        end
    end
end

InputCuboidsFile = "CuboidCoupling/Input_Cuboids.txt"
OutputPorepressureFile = "CuboidCoupling/Input_PorePressure.txt"



function main(InputCuboidsFile, OutputPorepressureFile)
    # Load Reservoir Centers
    Cuboids_Count, Cuboids_Center, Cuboids_Length = Load_Reservoir_Cuboids(InputCuboidsFile)
    println("---- Cuboids count is:  ", Cuboids_Count, " ----")

    # Time for External Stresses
    Time_Count =  2
    ExternalStress_TimeArray = zeros(Time_Count)
    ExternalStress_TimeArray[1]  = 0
    ExternalStress_TimeArray[2] = 10*365*24*3600 # 10 years in seconds

    # Pore Pressure Change
    PorePressureChange = zeros(Time_Count, Cuboids_Count)
    PorePressureChange[1,:] .= 0.0
    PorePressureChange[2,:] .= -30e6 # -30 MPa

    # Save Pore Pressure Change to a txt file
    Save_time_pressure_to_txt(OutputPorepressureFile, ExternalStress_TimeArray, PorePressureChange)

    println("---- TXT OutputFile Saved: ", OutputPorepressureFile, " ----")
end



main(InputCuboidsFile, OutputPorepressureFile)

