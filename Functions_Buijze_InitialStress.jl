using Interpolations

# function Set_InitialStress_IntialPoroPressure()
Rho_rock = 2400
Rho_water = 1150
Rho_gas = 200
grav = 9.8
K0 = 0.75 # (Buijze et al., 2019): K0=0.75
UpperBottomDepth = 1200

function CalculateInitialPoroPressure(z)
    if z >= 3000
        return 35e6 + Rho_water*grav*(z-3000)
    elseif z >= 2800 # z will never touch the PP boundary
        return  35e6 - Rho_gas*grav*(3000-z)
    else 
        return Rho_water*grav*z
    end
end


function CalculateInitialSigmaV(z)
    return 23.6e6 + Rho_rock*grav*(z-UpperBottomDepth)
end


function CalculateInitialSigmaH(z)
    return K0*CalculateInitialSigmaV(z)
end



InitialStressZ = range(2000., 4000., step =1)
InitialPorePressure = zeros(size(InitialStressZ))
InitialSigmaH = zeros(size(InitialStressZ))
InitialSigmaV = zeros(size(InitialStressZ))


for (idx,depth) in enumerate( InitialStressZ )
    InitialPorePressure[idx] = CalculateInitialPoroPressure(depth)
    InitialSigmaV[idx] = CalculateInitialSigmaV(depth)
    InitialSigmaH[idx] = CalculateInitialSigmaH(depth)
end

Interp_SigmaV = linear_interpolation(InitialStressZ, InitialSigmaV )
Interp_SigmaH = linear_interpolation(InitialStressZ, InitialSigmaH )
Interp_PorePressure = linear_interpolation(InitialStressZ, InitialPorePressure )

function CalculateInitialStress_on_Fault(phi, depth)
    """
    sigma_v > sigma_h
    """
    sigma_v = Interp_SigmaV(depth)
    sigma_h = Interp_SigmaH(depth)
    pore_pressure = Interp_PorePressure(depth)

    mohr_radius = (sigma_v - sigma_h)/2
    stress_normal = sigma_h + mohr_radius*(1 - cosd(180 - 2*phi))
    stress_shear = mohr_radius*sind(180 - 2*phi)
    stress_eff_normal = stress_normal - pore_pressure
    return stress_eff_normal, stress_shear, pore_pressure
end



    # Depth_test = 3000
    # println(  CalculateInitialStress_on_Fault(Phi, Depth_test) )
