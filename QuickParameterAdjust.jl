using DelimitedFiles
using JLD2
using PolyLog

include("./Functions_Buijze_InitialStress.jl")

function ParameterAdj(LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, 
    Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, Fault_V_Const,
    FaultStrikeAngle, FaultDipAngle, FaultCenter, Fault_BulkIndex, FaultLLRR, MinimumNormalStress)
    
    # >>>> Loading parameters >>>>>
    LoadingBuijzeGeomery = "Input_Buijze19_Geometry.jld2"
    Phi = load(LoadingBuijzeGeomery, "Phi") # deg


    FaultMass_Original = copy(FaultMass)
    Fault_a_Original = copy(Fault_a)
    Fault_b_Original = copy(Fault_b)
    Fault_Dc_Original = copy(Fault_Dc)
    Fault_Theta_i_Original = copy(Fault_Theta_i)
    Fault_V_i_Original = copy(Fault_V_i)
    Fault_Friction_i_Original = copy(Fault_Friction_i)
    Fault_NormalStress_Original = copy(Fault_NormalStress)
    Fault_V_Const_Original = copy(Fault_V_Const)
    FaultCenter_Original = copy(FaultCenter)
    FaultCount = length(Fault_a)
    MinimumNormalStress_Original = copy(MinimumNormalStress)  
    Count=0; 
    FaultIndex_Adjusted=0


    ######################################################################################################
    ############################    Alleviate the fault tip stress change ################################

    # ShearAllow = 1.2
    # NormalAllow = 1.2
    # FaultCount=length(FaultMass)
    # for i=1:FaultCount-LoadingFaultCount
    #     for j=1:FaultCount-LoadingFaultCount
    #             if i!=j
    #             StressTransferRatio = abs((StiffnessMatrixShear[j,i]- 0.6 * StiffnessMatrixNormal[j,i])/StiffnessMatrixShear[i,i])
    #             StressTransferRatioShear = abs((StiffnessMatrixShear[j,i])/StiffnessMatrixShear[i,i])
    #             StressTransferRatioNormal = abs((StiffnessMatrixNormal[j,i])/StiffnessMatrixShear[i,i])

    #                 if StressTransferRatioShear>ShearAllow
    #                     Count=Count+1
    #                     StiffnessMatrixShear[j,i]=StiffnessMatrixShear[i,i]*ShearAllow
    #                     FaultIndex_Adjusted=[FaultIndex_Adjusted;i]
    #                 end
                    
    #                 if StressTransferRatioNormal>NormalAllow
    #                     Count=Count+1
    #                     StiffnessMatrixNormal[j,i]=StiffnessMatrixShear[i,i]*NormalAllow
    #                     FaultIndex_Adjusted=[FaultIndex_Adjusted;i]
    #                 end
    #         end
    #     end
    # end
    # FaultIndex_Adjusted=FaultIndex_Adjusted[2:end]
    # FaultIndex_Adjusted=unique(FaultIndex_Adjusted)
    # println(FaultIndex_Adjusted)

    if FaultIndex_Adjusted == 0
        println("Adjusted Stiffness Count: 0")
    else
        println("Adjusted Stiffness Count: ", length(FaultIndex_Adjusted))
    end

    ##########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#############
    ######################################################################################################


    ######################################################################################################
    ##########################  Calculation of initial state from stress orientation #####################
    
    # MaxStressOrientation = 85. # between 0-180 degree
    # StressRatioMaxOverMin = 0.5
    # MinFrictionAllowed = 0.1 # smaller than this friction is not allowed

    # StressGradAtMaxOrientation = 6000.0
    # SurfaceStressAtMaxOrientation = 2e6
    # Fault_Theta_i .= 1e10
    # Fault_V_i .= 0.0
    # Friction_0 = ones(FaultCount) * 0.30
    # V0=1e-9;

    # Fault_Friction_i, Fault_NormalStress, Fault_V_i, Fault_Theta_i = 
    #             StressDependentFrictionParameters(MaxStressOrientation, StressRatioMaxOverMin, MinFrictionAllowed,
    #             StressGradAtMaxOrientation, SurfaceStressAtMaxOrientation,
    #             FaultStrikeAngle, FaultDipAngle, Fault_V_i, Fault_Theta_i, Fault_Friction_i, FaultLLRR,
    #             Fault_a, Fault_b, Fault_Dc, Fault_NormalStress, Friction_0, FaultCenter)

    ##########^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#############
    ######################################################################################################
    
    

    ######################################################################################################
    ######################################### Direct Adjust ##############################################
    # for i in eachindex(Fault_Dc)
    #     if Fault_BulkIndex[i]==1
    #         Fault_a[i]=0.003
    #     end
    # end
    
    
    # >>>> Set Initial Normal and Shear Stress >>>>
    
    for i=1:FaultCount
        fault_patch_depth = FaultCenter[i,3]
        intialstress_normal, intialstress_shear, initial_porepressure = CalculateInitialStress_on_Fault(Phi, fault_patch_depth)
        Fault_NormalStress[i] = intialstress_normal
        Fault_Friction_i[i] = intialstress_shear/intialstress_normal
    end

    Fault_fr = 0.35
    Fault_fs = 0.60
    Fault_Wsw = 0.00375 # unit: mm^-1
    Fault_Dsw = (Fault_fs-Fault_fr)/Fault_Wsw*1e-3 # unit: m
    Fault_Vr = 1e-1
    Fault_Vp = 1e-6 
    Fault_V0 = 1e-9
    
    Fault_Gc = 0.5*(Fault_fs-Fault_fr)*Fault_Dsw 
    Fault_Dc = (Fault_Gc)./(Fault_b_Original*( -reli2(1-exp((Fault_fs-Fault_fr)/Fault_b_Original[1])) ))
    Fault_f0 = Fault_fr.-(Fault_a_Original-Fault_b_Original).*log(Fault_Vr/Fault_V0)
    Fault_Theta_i = Fault_Dc/Fault_V0.*exp.( (Fault_fs.-Fault_f0.-Fault_a_Original*log(Fault_Vp/Fault_V0))./Fault_b_Original )
    Fault_V_i = Fault_V0*exp.( (Fault_Friction_i.-Fault_f0.-Fault_b_Original.*log.(Fault_Theta_i*Fault_V0./Fault_Dc))./Fault_a_Original )
    
    println( "Size of Dc:  ", size(Fault_Dc),"  ", Fault_Dc[1] )
    println( "Size of theta_i:  ", size(Fault_Theta_i),"  ", Fault_Theta_i[1] )
    println( "Size of V_i:  ", size(Fault_V_i),"  ", maximum(Fault_V_i) )
    println( "Size of f0:  ", size(Fault_f0), "  ", Fault_f0[1] )
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    
    
    # Fault_Dc .= 1e-4 # 1e-3 
    # Fault_a .= 0.05
    # Fault_b .= 0.003
    # Fault_NormalStress .= 20e6
    # Fault_V_i .= 1e-10
    # Fault_Theta_i .= 1e3 # 1.5027018579405773e9
    # MinimumNormalStress= 1e6
    # FaultMass .= 1e6


    
    ###^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^###
    #####################################################################################################


    if FaultMass_Original != FaultMass
        println("- Fault Mass Adjusted")
    end
    
    if Fault_a_Original != Fault_a
        println("- Fault a Adjusted")
    end
    if Fault_b_Original != Fault_b
        println("- Fault b Adjusted")
    end
    if Fault_Dc_Original != Fault_Dc
        println("- Fault Dc Adjusted")
    end
    if Fault_Theta_i_Original != Fault_Theta_i
        println("- Fault Theta_i Adjusted")
    end
    if Fault_V_i_Original != Fault_V_i
        println("- Fault Fault_V_i Adjusted")
    end
    if Fault_Friction_i_Original != Fault_Friction_i
        println("- Fault Friction_i Adjusted")
    end
    if Fault_NormalStress_Original != Fault_NormalStress
        println("- Fault Normal stress Adjusted")
    end
    if Fault_V_Const_Original != Fault_V_Const
        println("- Fault V_Const Adjusted")
    end
    if FaultCenter_Original != FaultCenter
        println("- Fault Center Adjusted")
    end

    if MinimumNormalStress_Original != MinimumNormalStress
        println("- Minimum NormalStress Adjusted")
    end

    return LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, Fault_Theta_i, Fault_V_i, 
    Fault_Friction_i, Fault_NormalStress, Fault_V_Const, FaultCenter, FaultIndex_Adjusted, MinimumNormalStress


end



function StressDependentFrictionParameters(MaxStressOrientation, StressRatioMaxOverMin, MinFrictionAllowed,
    StressGradAtMaxOrientation, SurfaceStressAtMaxOrientation,
    FaultStrikeAngle, FaultDipAngle, Fault_V_i, Fault_Theta_i, Fault_Friction_i, FaultLLRR,
    Fault_a, Fault_b, Fault_Dc, Fault_NormalStress, Friction_0, FaultCenter)
    V0 = 1e-9
    NormalStressParameter = (1+StressRatioMaxOverMin)/2 .+ (1-StressRatioMaxOverMin)/2 .* cosd.(2 .* (FaultStrikeAngle .- 90.0 .- MaxStressOrientation))
    ShearStressParameter = -(1-StressRatioMaxOverMin)/2 .* sind.(2 .* (FaultStrikeAngle .- 90.0 .- MaxStressOrientation))
    Fault_Friction_i .= abs.(ShearStressParameter ./ NormalStressParameter)
    
    if FaultLLRR != sign.(ShearStressParameter)
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("!!! Warning: Stress Orientation inconsistant with Slip Direction !!! ")
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    end
        
    for i in eachindex(Fault_Friction_i)
        if Fault_Friction_i[i] < MinFrictionAllowed
            Fault_Friction_i[i] = MinFrictionAllowed
        end
        Fault_NormalStress[i] = (StressGradAtMaxOrientation * FaultCenter[i,3] + SurfaceStressAtMaxOrientation) * NormalStressParameter[i] 
    end
    
    if iszero(Fault_V_i)
        Fault_V_i = V0 .* exp.( (Fault_Friction_i .- Friction_0 .- Fault_b .* log.(Fault_Theta_i .* V0./Fault_Dc)) ./ Fault_a)
    end
    
    
    if iszero(Fault_Theta_i)
        Fault_Theta_i = Fault_Dc ./ V0 .* exp.( (Fault_Friction_i .- Friction_0 .- Fault_a .* log.(Fault_V_i ./ V0)) ./ Fault_b)
    end
    
    return Fault_Friction_i, Fault_NormalStress, Fault_V_i, Fault_Theta_i
end


