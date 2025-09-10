# Flow2QuakeDFN
> [!WARNING]
> **Early Stage Code.** This is not a formal release. We're still working on a refined version. 
> 

## Introduction
Flow2QuakeDFN, based on the latest version (1.2.0) of [Quake-DFN](https://github.com/limkjae/Quake-DFN), simulates earthquake rupture driven by changes in reservoir pore pressure. The underlying physics involves the coupling between reservoir pore pressure, reservoir poro-elastic stress, and fault rate-and-state rupture dynamics. This code can generate both injection-induced  ([video 1](#injection-rupture)) and depletion-induced ([video 2](#depletion-rupture)) earthquake rupture. 


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

We use `Example`, a simpler depletion-induced case, to demonstrate the workflow. We first generate multiple cuboids that consist the reservoir. A pore pressure generator then gives the spatial heterogeneous and temporal evolving pore pressure in each cuboid. After assigning the fault geometry, elastic modulus and rupture mode, a poro-elastic stress calculator computes the corresponding external shear and normal stresses on each fault patch. Finally, we set the initial stress and initial rate-and-state parameters on fault before running the quasi-dynamic Quake-DFN simulation.  


### Reservoir Cuboids
Create the cuboids for the reservoir. In this example, the reservoir consists of a single cuboid.
```
julia CuboidCoupling/Example/Generate_Cuboids.jl
```

This generates the file containing the position and size of each cuboid `CuboidCoupling/Input_Cuboids.txt`.

<img width="459" height="83" alt="image_cuboids" src="https://github.com/user-attachments/assets/e7d5ec23-11da-4803-8f01-b80bbedf4570" />




### Reservoir Pore Pressure
Assign the temporal pore pressure change in each reservoir cuboid and the corresponding external time. 

```
julia CuboidCoupling/Example/Generate_PorePressure.jl
``` 

The external time and pore pressure changes are stored in  `CuboidCoupling/Input_PorePressureChange.txt`.
<img width="256" height="106" alt="image_porepressure" src="https://github.com/user-attachments/assets/2cc88289-4bca-4f1e-ad0b-39d667a011d6" />

> [!NOTE]
> You can use your own mesh-generator and pore pressure solver code to produce `CuboidCoupling/Input_Cuboids.txt` and `CuboidCoupling/Input_PorePressureChange.txt`, as long as their formats remain the same as provided.



### Fault Geometry
Build the bulk fault geometry, including fault geometry, elastic modulus and rupture mode. In this example, the fault is a 70° dipping normal fault that cuts the reservoir with no offset.

```
julia CuboidCoupling/Example/Build_Geometry_Single_Normal.jl
```

Discretize the fault, pre-calculate and store the static stress interaction between each fault patch.
```
julia RUN_BUILD_INPUT.jl
```


### Fault Poro-elastic Stress Change



```
julia CuboidCoupling/CalculatePoroElasticStress.jl
```

This file reads `CuboidCoupling/Input_Cuboids.txt`, `CuboidCoupling/Input_PorePressureChange.txt` and `Input_Discretized.jld2`, calculates the resulting poro-elastic external stress changes on each fault patch at each (prescribed) time, and  generates `Input_ExternalStressChange.jld2`, which contains the computed external fault shear and effective normal stress and is loaded when running the Quake-DFN simulation. In addition, you need to manually specify the criterion for which fault patches experience pore pressure and therefore have effective normal stress $\sigma' = \sigma - P$. 

<img width="529" height="179" alt="image_poro" src="https://github.com/user-attachments/assets/572efd52-e8c7-4e9a-b71c-4e40931ad4c5" />


### Fault Initial Condition
The initial fault conditions can be set in either `CuboidCoupling/Example/Build_Geometry_Single_Normal.jl` or `QuickParameterAdjust.jl`. We recommend systematically set all initial fault stress and rate-and-state parameters in `QuickParameterAdjust.jl`, since it overwrites the initial condition written in the former file. This file is executed automatically when running the Quake-DFN simulation.

<img width="746" height="474" alt="image_initial" src="https://github.com/user-attachments/assets/de017542-32ba-4f07-aba5-bbaa4c067c10" />


### Quake-DFN Simulation
```
julia RUN_QUAKEDFN.jl
``````
