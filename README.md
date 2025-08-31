# Flow2QuakeDFN
> [!WARNING]
> **Early Stage Code.** This is not a formal release. We're still working on a refined version. 
> 
> **Outdated Quake-DFN Version.** This is based on [Quake-DFN](https://github.com/limkjae/Quake-DFN) V1.0, not the most recent version of [Quake-DFN](https://github.com/limkjae/Quake-DFN).

## Introduction
This project simulates earthquake rupture driven by changes in reservoir pore pressure. The physics underlying involves the coupling between reservoir pore pressure, poro-elastic stress, and rate-and-state fault dynamics. We use the case `Example` to demonstrate the workflow: fluid is injected into a flat reservoir and induces rupture on a planar, dipping fault offset from the injection point. The movie below illustrates the coupled processes and is produced at the end of the workflow.

<video controls src="injection_rupture.mp4" title="Title"></video>

https://github.com/user-attachments/assets/70097d5e-351f-428e-918d-bce94bcdc325

## Dependencies
- **Python** 3.10+
- Packages: numpy, matplotlib, pyvista, tqdm, h5py
- **Julia** 1.11+
- Packages: Pkg, PyPlot, PyCall, Conda, DelimitedFiles, JLD2, LinearAlgebra, Printf, SpecialFunctions, StaticArrays, LowRankApprox, Distributed, Statistics, Clustering, PolyLog, ProgressBars, HDF5, CSV, DataFrames, WriteVTK

A detailed description of how to install the Julia packages is in [QuakeDFN_UserGuide.pdf](QuakeDFN_UserGuide.pdf)

## Workflow
### Reservoir Mesh 
Create a computational mesh for the reservoir. In this example, the reservoir consists of a 40 × 40 grid of **cuboid** elements.
```
python ReservoirMeshExamples/Example_Generate_Cubes.py
```


You can check the positions of all the reservoir cubes via the *Input_Example_Cubes.txt* file:

![alt text](image_cubes.png)

### Reservoir Pore Pressure
Assign the time and pore pressure change for each reservoir cuboid. 

```
julia PorePressureExamples/Example_Injection.jl
``` 

You can check the pore pressure changes at each time in each reservoir cube in *Input_Example_PorePressureChange.csv*.

![alt text](image_porepressure.png)


### Fault Geometry
Build the bulk fault geometry.

```
julia InputGeometryExamples/Example_BuildGeometry_Single_Normal.jl
```

Discretize the fault and generate the fault meshes.
```
julia RUN_BUILD_INPUT.jl
```

### Fault Initial Stress
Set the initial shear&normal stresses on each fault patch in *QuickParameterAdjust.jl*. This file is executed automatically when running the QuakeDFN simulation `RUN_QUAKEDFN.jl`.

### Fault Poro-elastic Stress Change
Compute the resulting poroelastic stress changes on each fault patch at each (prescribed) time step. 


```
julia ExternalStressCalculationExamples/Example_PoroStress.jl
```

This generates the *Input_ExternalStressChange.jld2* file, which contains the computed fault external shear and effective normal stress. *Plot_InitialStress.ipynb* file provides the necessary code to visualize these external stresses. 

![alt text](image_tau.png)
![alt text](image_sigma.png)


### QuakeDFN Simulation
```
julia RUN_QUAKEDFN.jl
```

### Post-processing (Paraview)
Write results to VTK files for visualization in ParaView.
```
julia Plot_ResultWrite2VTK.jl
```
