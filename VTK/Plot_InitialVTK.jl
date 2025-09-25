using JLD2
using WriteVTK
using StaticArrays

Input_Fault_FileName = "Input_Discretized.jld2"
Input_Stress_Filename = "Input_ExternalStressChange.jld2"
Output_VTK_Filename = "./TetrahedronCoupling/ExternalStress"

FaultCenter= load(Input_Fault_FileName, "FaultCenter")
ExternalStress_Normal = load(Input_Stress_Filename, "ExternalStress_Normal")
ExternalStress_Shear = load(Input_Stress_Filename, "ExternalStress_Shear")

FaultCenter[:,3] = - FaultCenter[:,3]
coordinates = [SVector(row...) for row in eachrow(FaultCenter)]

vtk = vtk_grid(Output_VTK_Filename, coordinates)

vtk["sigma"] = ExternalStress_Normal[end,:]
vtk["tau"] = ExternalStress_Shear[end,:]

vtk_save(vtk)