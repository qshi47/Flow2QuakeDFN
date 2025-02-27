# Flow2QuakeDFN
Flow2QuakeDFN aims at simulating the fault ruptures in realistic stresses region.

[TOC]

## Step1: Initialize fault 
    julia ./InputGeometryExamples/Buijze19_BuildGeometry_Single_Normal.jl

## Step 2: Discritize fault & Calculate Okada GF functionbetween fault patches
    julia ./RUN_BUILD_INPUT.jl

## Step 3: Generate cubes of the reservoir.
    julia ./Generate_Mesh_Buijze.py

## Step 4: Calculate external poroelastic stresses on each fault patch
    julia ./ExternalStressCalculationExamples/Buijze_PoroStress.jl

## Step 5: Run QuakeDFN
    julia ./RUN_QUAKEDFN.jl

## Step 6: Post-processing: plot result slip/friction evolution 
    Plot_ResultSlip.ipynb
