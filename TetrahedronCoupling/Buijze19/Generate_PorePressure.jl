using DelimitedFiles


function Save_time_pressure_to_txt(filename::String, time_array::AbstractVector, pressure::AbstractMatrix)
    output_matrix = hcat(time_array, pressure)
    blocksCount = size(pressure, 2)

    # Construct header
    header = ["time(sec)"]
    for i in 1:blocksCount
        push!(header, "cuboid$(i)")
    end

    # Open file and write header + data
    open(filename, "w") do io
        write(io, join(header, ", "))           # write header line
        write(io, "\n")                        # new line after header
        writedlm(io, output_matrix, ", ")          # write data rows
    end
end

function main(InputTetrahedronFile, OutputPorepressureFile)
    # Load Reservoir Tetrahedrons
    data = readdlm(InputTetrahedronFile, ',', skipstart=1)
    Tetrahedrons_Count = Int(size(data, 1) / 4)
    println("---- Tetrahedrons count is:  ", Tetrahedrons_Count, " ----")

    # Time for External Stresses
    Time_Count =  2
    ExternalStress_TimeArray = zeros(Time_Count)
    ExternalStress_TimeArray[1]  = 0
    ExternalStress_TimeArray[2] = 10*365*24*3600 # 10 years in seconds

    # Pore Pressure Change
    PorePressureChange = zeros(Time_Count, Tetrahedrons_Count)
    PorePressureChange[1,:] .= 0.0
    PorePressureChange[2,:] .= -30e6 # -30 MPa

    # Save Pore Pressure Change to a txt file
    Save_time_pressure_to_txt(OutputPorepressureFile, ExternalStress_TimeArray, PorePressureChange)

    println("---- TXT OutputFile Saved: ", OutputPorepressureFile, " ----")
end


InputTetrahedronFile = "TetrahedronCoupling/Input_Tetrahedron.txt"
OutputPorepressureFile = "TetrahedronCoupling/Input_PorePressure.txt"

main(InputTetrahedronFile, OutputPorepressureFile)

