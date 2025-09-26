using JLD2
using WriteVTK
using StaticArrays
using ProgressBars

function write_pvd(pvd_filename, vts_filenames, result_time)
    """
    Write a `.pvd` file to link multiple `.vts` files for time series data.

    Arguments:
    - pvd_filename: Name of the output PVD file.
    - vts_filenames: Array of VTS file names.
    - result_time: Vector of time values.
    """
    open(pvd_filename, "w") do io
        println(io, """<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">""")
        println(io, "  <Collection>")
        for (i, time) in enumerate(result_time)
            println(io, """    <DataSet timestep="$time" file="$(vts_filenames[i])"/>""")
        end
        println(io, "  </Collection>")
        println(io, "</VTKFile>")
    end
end




function delete_all_files(folder_path::String)
    # Check if the folder exists
    # List all files in the folder
    files = readdir(folder_path)    
    # Loop through each file and delete it
    for file in files
        file_path = joinpath(folder_path, file)
        if isfile(file_path)
            rm(file_path)
        end
    end    
    println("All files in '$folder_path' have been deleted.")
end

function create_new_folder(folder_path::String)
    if !isdir(folder_path)  # Check if the folder already exists
        mkdir(folder_path)  # Create the folder
        println("Folder created: $folder_path")
    else
        println("Folder already exists: $folder_path")
    end
end

function write_time_series_vtk(base_filename, coordinates, result_data, result_time)
    """
    Write a series of `.vts` files for time-dependent data and a `.pvd` file to link them.

    Arguments:
    - base_filename: Base name for the output files (e.g., "output").
    - coordinates: Array of coordinates for grid points.
    - result_data: Dictionary containing scalar data arrays for each time step.
    - result_time: Vector of time values.
    """
    Year_in_Second = 365*24*3600 # seconds

    nt = length(result_time)  # Number of time steps
    vts_filenames = []

    @info "Start saving vtk files ----"
    for t in tqdm( 1:nt, unit = "time step" )
        vtk_filename = "$(base_filename)step_$(t).vts"
        push!(vts_filenames, vtk_filename)

        vtk = vtk_grid(vtk_filename, coordinates)

        # Add scalar data for this time step
        vtk["ResultTime"] = result_time[t]/Year_in_Second
        vtk["ResultV"] = result_data["ResultV"][t, :]
        # vtk["ResultDisp"] = result_data["ResultDisp"][t, :]
        # vtk["Result_NormalStress"] = result_data["Result_NormalStress"][t, :]
    
        vtk_save(vtk)
    end

    # Write the PVD file to link all time steps
    write_pvd(base_filename * ".pvd", vts_filenames, result_time)
end

function write_result_to_vtk_files(OutputFolder, LoadingInputFileName, DataFileName)
    if isdir(OutputFolder)
            delete_all_files(OutputFolder)
        else
            create_new_folder(OutputFolder)
    end

    FaultCenter= load(LoadingInputFileName, "FaultCenter")
    ResultTime, ResultV, ResultDisp, Result_NormalStress, History_Theta =
    load(DataFileName,"History_Time", "History_V", "History_Disp", "History_NormalStress","History_Theta")


    # Coordinates (array of SVector)
    FaultCenter[:,3] = - FaultCenter[:,3]
    coordinates = [SVector(row...) for row in eachrow(FaultCenter)]

    # Prepare result data
    result_data = Dict(
        "ResultTime" => ResultTime,
        "ResultV" => ResultV )
        # "ResultDisp" => ResultDisp,
        # "Result_NormalStress" => Result_NormalStress
        

    # Write the VTK time series and PVD file
    write_time_series_vtk(OutputFolder, coordinates, result_data, ResultTime)
end


Folder = "VTK/" # "/home/qshi2/tecto-data1/"
FileString = "Fault_Point_vtu" 

LoadingInputFileName = "Input_Discretized.jld2" 
DataFileName = "Results/Result.jld2"

write_result_to_vtk_files( Folder * FileString * "/", LoadingInputFileName, DataFileName)