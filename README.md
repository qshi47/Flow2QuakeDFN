# Flow2QuakeDFN
> [!WARNING]
> **Early Stage Code.** This is not a formal release. We're still working on a refined version. 
> 

## Introduction
Flow2QuakeDFN, based on the latest version (1.2.0) of [Quake-DFN](https://github.com/limkjae/Quake-DFN), simulates earthquake rupture driven by changes in reservoir pore pressure. The underlying physics involves the coupling between reservoir pore pressure, reservoir poro-thermo-elastic stress, and fault rate-and-state rupture dynamics. This code can generate both injection-induced  ([video 1](#injection-rupture)) and depletion-induced ([video 2](#depletion-rupture)) earthquake rupture. 


### Injection-induced earthquake
<a id="injection-rupture"></a>
<video controls src="https://github.com/user-attachments/assets/70097d5e-351f-428e-918d-bce94bcdc325" title="Injection rupture video"></video>


### Depletion-induced earthquake
<a id="depletion-rupture"></a>
<video controls src="https://github.com/user-attachments/assets/8d35f1ac-b337-4a46-9f0b-f034e05bdcc9" title="Depletion rupture video"></video>







## Dependencies
- **Julia:** 1.11+
- **Packages:** Pkg, PyPlot, PyCall, Conda, DelimitedFiles, JLD2, LinearAlgebra, Printf, SpecialFunctions, StaticArrays, LowRankApprox, Distributed, Statistics, Clustering, PolyLog, ProgressBars, HDF5, CSV, DataFrames, WriteVTK

A detailed tutorial on how to install Julia and Julia packages is in [QuakeDFN_UserGuide_V1.2.pdf](QuakeDFN_UserGuide_V1.2.pdf)

## Workflow

We use the example case `Single`, a simpler depletion-induced case, to demonstrate the workflow. We first generate multiple cuboids that consist the reservoir. A pore pressure generator then gives the spatial heterogeneous and temporal evolving pore pressure in each cuboid. After assigning the fault geometry, elastic modulus and rupture mode, a poro-thermo-elastic stress calculator computes the corresponding external shear and normal stresses on each fault patch. Finally, we set the initial stress and initial rate-and-state parameters on fault before running the quasi-dynamic Quake-DFN simulation.  


### Reservoir Cuboids
Create the cuboids for the reservoir. In this example, the reservoir consists of a single cuboid.
```
julia CuboidCoupling/Single/Generate_Cuboids.jl
```

This generates the file containing the position and size of each cuboid `CuboidCoupling/Input_Cuboids.txt`.

<img width="459" height="83" alt="image_cuboids" src="https://github.com/user-attachments/assets/e7d5ec23-11da-4803-8f01-b80bbedf4570" />




### Reservoir Pore Pressure
Assign the temporal change in pore pressure and temperature  for each reservoir cuboid and the corresponding external time. 

```
julia CuboidCoupling/Single/Generate_PorePressure_Temperature.jl
``` 

The external time and pore pressure changes are stored in  `CuboidCoupling/Input_PorePressure.txt`.
Similarly, external temperature change is stored in `CuboidCoupling/Input_Temperature.txt`.


> [!NOTE]
> You can use your own mesh-generator and pore pressure solver code to produce `CuboidCoupling/Input_Cuboids.txt`, `CuboidCoupling/Input_PorePressure.txt` and `CuboidCoupling/Input_Temperature.txt`, as long as their formats remain the same as provided.



### Fault Geometry
Build the bulk fault geometry, including fault geometry, elastic modulus and rupture mode. In this example, the fault is a 70° dipping normal fault that cuts the reservoir with no offset.

```
julia CuboidCoupling/Single/Build_Geometry_Single_Normal.jl
```

Discretize the fault, pre-calculate and store the static stress interaction between each fault patch.
```
julia RUN_BUILD_INPUT.jl
```


### Fault Poro-elastic Stress Change



```
julia CuboidCoupling/CalculatePoroElasticStress.jl
```

This file reads `Input_Cuboids.txt`, `Input_PorePressure.txt`,`Input_Temperature.txt` and `Input_Discretized.jld2`, calculates the resulting poro-thermo-elastic external stress changes on each fault patch at each (prescribed) time, and ends up with generating `Input_ExternalStressChange.jld2`, which contains the computed external fault shear and effective normal stress and is loaded when running the Quake-DFN simulation.

### Fault Initial Condition
The initial fault conditions can be set in either `CuboidCoupling/Single/Build_Geometry_Single_Normal.jl` or `QuickParameterAdjust.jl`. We recommend systematically set all initial fault stress and rate-and-state parameters in `QuickParameterAdjust.jl`, since it overwrites the initial condition written in the former file. This file is executed automatically when running the Quake-DFN simulation.

<img width="746" height="474" alt="image_initial" src="https://github.com/user-attachments/assets/de017542-32ba-4f07-aba5-bbaa4c067c10" />


### Quake-DFN Simulation
```
julia RUN_QUAKEDFN.jl
``````
