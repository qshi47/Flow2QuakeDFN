using DelimitedFiles
using JLD2
using PolyLog


function SetInitialStress_Example!(FaultCount, fault_normalStress, fault_friction_i)
    for i=1:FaultCount
        fault_normalStress[i] = 50e6
        fault_friction_i[i] = 0.6
    end

    return fault_normalStress, fault_friction_i
end

function SetInitialRSParameters_Example!(Fault_a_Original, Fault_b_Original, Fault_Friction_i)
    faultcount = size(Fault_a_Original)
    Fault_V_i = zeros(faultcount)
    Fault_f0 = zeros(faultcount)
    Fault_Dc = zeros(faultcount)
    Fault_Theta_i = zeros(faultcount)

    Fault_V0 = 1e-9
    Fault_V_i .= 1e-10
    Fault_Dc .= 1e-4 
    Fault_Theta_i .= 1e3
    
    """ 
    # Fault_f0 .= 0.4
    # Fault_V_i = Fault_V0*exp.( (Fault_Friction_i.-Fault_f0.-Fault_b_Original.*log.(Fault_Theta_i*Fault_V0./Fault_Dc))./Fault_a_Original )
    """
    
    println( "Size of V_i:  ", size(Fault_V_i),"  ", maximum(Fault_V_i) )
    return Fault_Dc, Fault_Theta_i, Fault_V_i

end

function ParameterAdj(LoadingFaultCount, FaultMass, Fault_a, Fault_b, Fault_Dc, 
    Fault_Theta_i, Fault_V_i, Fault_Friction_i, Fault_NormalStress, Fault_V_Const,
    FaultStrikeAngle, FaultDipAngle, FaultCenter, Fault_BulkIndex, FaultLLRR, MinimumNormalStress)
    

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

    #####################################################################################################
    ###^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^###
    
    # >>>> Set Initial Normal and Shear Stress >>>>
    Fault_NormalStress, Fault_Friction_i  = SetInitialStress_Example!(FaultCount, Fault_NormalStress, Fault_Friction_i)

    # >>>> Set Initial R&S Parameters >>>>
    Fault_Dc, Fault_Theta_i, Fault_V_i = SetInitialRSParameters_Example!(Fault_a_Original, Fault_b_Original, Fault_Friction_i)

    
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


