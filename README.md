# Flow2QuakeDFN

## Introduction
This project focuses on simulating earthquake rupture in response to the change in pore pressure of reservoir. The physics underlying this process involves the coupling between pore pressure, poro-elastic stress, and rate-and-state fault dynamics. We'll use the case `Example` to demonstrate the whole workflow, where fluid is injected to a flat reservoir and induce earthquake rupture on a flat dipping fault offsetting to the injection point. The below movie showcase the above physical process and will be produced at the end of the workflow.

<video controls src="injection_rupture.mp4" title="Title"></video>


https://github.com/user-attachments/assets/70097d5e-351f-428e-918d-bce94bcdc325

## Dependencies
- Python 3.10+
- Packages: numpy, matplotlib, pyvista, tqdm, h5py
- Julia 1.11+
- Packages: Pkg, PyPlot, PyCall, Conda, DelimitedFiles, JLD2, LinearAlgebra, Printf, SpecialFunctions, StaticArrays, LowRankApprox, Distributed, Statistics, Clustering, ProgressBars, HDF5, CSV, DataFrames, WriteVTK

## Workflow
### Reservoir Mesh 
Create a computational mesh for the reservoir. In this example, the reservoir consists of 40*40 **cuboids** meshes.

```
python ReservoirMeshExamples/Example_Generate_Cubes.py
```


You can check the positions of all the reservoir cubes via the *Input_Example_Cubes.txt* file:

![alt text](image_cubes.png)

### Reservoir Pore Pressure
Assign the time and changes in pore pressure within each reservoir cuboid. 

```
julia PorePressureExamples/Example_Injection.jl
``` 

You can check the pore pressure changes at each time in each reservoir cube in the *Input_Example_PorePressureChange.csv* file.

![alt text](image_porepressure.png)


### Fault Geometry
Build bulk fault geometry.

```
julia InputGeometryExamples/Example_BuildGeometry_Single_Normal.jl
```

Discretize the fault and generate the fault meshes.
```
julia RUN_BUILD_INPUT.jl
```

### Fault Initial Stress
The initial shear&normal stresses on each fault patch are set in *QuickParameterAdjust.jl* file. This file will be implemented automatically when running the QuakeDFN simulation `RUN_QUAKEDFN.jl`.

### Fault Poro-elastic Stress Change
Compute the resulting poroelastic stress changes on each fault patch at each (prescribed) time step. 

```
julia ExternalStressCalculationExamples/Example_PoroStress.jl
```

*Plot_InitialStress.ipynb* file provides the necessary codes to visualize the computed external shear and effective normal stress.
![alt text](image_tau.png)
![alt text](image_sigma.png)


### QuakeDFN Simulation
```
julia RUN_QUAKEDFN.jl
```

### Post-processing (Paraview)
```
julia Plot_ResultWrite2VTK.jl
```
