function Calculate_Correspondense_Matrix(block_center_pos, pressure_pos) 
    """
        Calculate_Correspondense_Matrix(block_center_pos, pressure_pos) 
    
        Return correspondense matrix {C} 
    C[i,j] means the ith cube's pressure change caused by the jth pressure point
    """
    tolerance = 500  #dim change to meters

    correspondance_matrix = zeros( (size(block_center_pos)[1], size(pressure_pos)[1]) )

    for ii in range(1, size(block_center_pos)[1], step = 1)
        #For each reservoir's block, we search the closest point of pressure_pos
        disp_to_each_pressure = sqrt.(  (pressure_pos[:,1] .- block_center_pos[ii,1]).^2 + (pressure_pos[:,2] .- block_center_pos[ii,2]).^2  )
        indx = argmin(disp_to_each_pressure) 
        val  = disp_to_each_pressure[indx] # minimum(disp_to_each_pressure)
        if val < tolerance
            correspondance_matrix[ii, indx] = 1.0
        end
    end

    return correspondance_matrix

end


function fFunc(x,y,z,R)
    return z*atan((x*y)/(z*R)) - x*log(abs(R + y)) - y*log(abs(R + x))
end


function Calculate_PoroDispStress_MultipleBlocks_GF( ObsPoint, Blocks, Blocks_Cm, 
    ShearModulus, PossionRatio )
    """
    输入坐标是Z轴向下的左手系的XYZ坐标
    所有Li2021中计算公式，假设成立于左手Z下系中，输出应力压缩为正
    输出位移遵循原始坐标系，应力仍然保持压缩为正
    ObsPoint: (1,3)
    Blocks: can be N (8,3)
    """
    
    pi = 3.1415

    Disp = [0.0, 0.0, 0.0]

    sig_11 = 0.0
    sig_12 = 0.0
    sig_13 = 0.0
    sig_22 = 0.0
    sig_23 = 0.0
    sig_33 = 0.0

    
    for (idx_block, block) in enumerate(Blocks)
        Prefactor_disp = Blocks_Cm[idx_block]/(4*pi)
        Prefactor_sig  = Blocks_Cm[idx_block]*ShearModulus/(2*pi)
        for jj in range( start=1, stop=size(block)[1] ) # block:(8,3) Mat
            vertx = block[jj,:] # (1,3)
            xbar  = vertx[1] -  ObsPoint[1]
            ybar  = vertx[2] -  ObsPoint[2]
            ζp = vertx[3] + ObsPoint[3]
            ζm = vertx[3] - ObsPoint[3]
            rp = sqrt(   xbar^2 + ybar^2 + ζp^2    )
            rm = sqrt(   xbar^2 + ybar^2 + ζm^2    )

            # minus sign comes from that Julia convention starting from 1, not 0! 
            Sσ = - [-1, 1, 1, -1, 1, -1, -1, 1]

            # -- Displacements --
            # same as (B. li et al., 2021)
            # Z downward positive 
            Disp[1] +=  Prefactor_disp*Sσ[jj]*(fFunc(ybar,ζm,xbar,rm) 
                                + (3. - 4*PossionRatio)*fFunc(ybar, ζp,  xbar, rp) 
                                + 2*ObsPoint[3]*log(abs(rp + ybar)))

            Disp[2] +=  Prefactor_disp*Sσ[jj]*(fFunc(xbar,ζm,ybar,rm)
                                + (3. - 4*PossionRatio)*fFunc(xbar, ζp,  ybar, rp)
                                + 2*ObsPoint[3]*log(abs(rp + xbar)))

            Disp[3] += - Prefactor_disp*Sσ[jj]*(fFunc(xbar,ybar,ζm,rm)
                                + (3. - 4*PossionRatio)*fFunc(xbar,ybar, ζp, rp) 
                                - 2 * ObsPoint[3] * atan((ζp*rp)/(xbar*ybar)) )

            # -- Stresses --
            # same as Flow2Quake and (B. li et al., 2021)
            # compression positive + Z downward positive
            
            sig_11 += Prefactor_sig*Sσ[jj]*(
                        - atan((xbar*rm)/(ybar*ζm))
                        - (3.0-4*PossionRatio)* atan((xbar*rp)/(ybar*ζp))
                        + 4*PossionRatio* atan((ζp*rp)/(xbar*ybar))
                        + 2*ObsPoint[3]*xbar*ybar/( xbar^2 +  ζp^2 )/rp)

            # σYY
            sig_22 += Prefactor_sig*Sσ[jj]*(
                        -  atan((ybar*rm)/(xbar*ζm))
                        - (3.0-4*PossionRatio)* atan((ybar*rp)/(xbar*ζp))
                        + 4*PossionRatio* atan((ζp*rp)/(xbar*ybar))
                        + 2*ObsPoint[3]*xbar*ybar/( ybar^2 + ζp^2 )/rp)

            # σZZ
            sig_33 +=  Prefactor_sig*Sσ[jj]*(
                            -  atan((ζm*rm)/(xbar*ybar))
                            +  atan((ζp*rp)/(xbar*ybar))
                            - 2*ObsPoint[3]*xbar*ybar/rp*(1/(xbar^2 +  ζp^2) + 1/(ybar^2 +  ζp^2)) )

            # σXY
            sig_12 += - Prefactor_sig*Sσ[jj]*( log(abs(rm + ζm))
                            + (3.0-4*PossionRatio)* log(abs(rp+ζp))
                            + 2*ObsPoint[3]/rp   )

            # σXZ 
            sig_13 +=  Prefactor_sig*Sσ[jj]*( log( abs((rm+ybar)/(rp+ybar)))
                            + 2*ObsPoint[3]*ybar*ζp/( xbar^2+ ζp^2)/rp)
            # σYZ
            sig_23 +=  Prefactor_sig*Sσ[jj]*( log( abs((rm+xbar)/(rp+xbar)))
                            + 2*ObsPoint[3]*xbar*ζp/(ybar^2+ ζp^2)/rp)
        end

    end
    
    # ----对比COMSOL结果，手动符号修正----
    Disp[3] *= -1
    sig_13 *= -1
    sig_23 *= -1
    # ------------------------------------
    Sigma = [sig_11, sig_22, sig_33, sig_12, sig_13, sig_23]
    
    return Disp, Sigma 
end


function Calculate_DispStress_MultipleBlocks( ObsPoint, Blocks, Blocks_Cm_P, 
    ShearModulus, PossionRatio )
    """
    输入坐标是Z轴向下的左手系的XYZ坐标
    所有Li2021中计算公式，假设成立于左手Z下系中，输出应力压缩为正
    输出位移遵循原始坐标系，应力仍然保持压缩为正
    ObsPoint: (1,3)
    Blocks: can be N (8,3)
    """
    
    pi = 3.1415

    Disp = [0.0, 0.0, 0.0]

    sig_11 = 0.0
    sig_12 = 0.0
    sig_13 = 0.0
    sig_22 = 0.0
    sig_23 = 0.0
    sig_33 = 0.0

    
    for (idx_block, block) in enumerate(Blocks)
        Prefactor_disp = Blocks_Cm_P[idx_block]/(4*pi)
        Prefactor_sig  = Blocks_Cm_P[idx_block]*ShearModulus/(2*pi)
        for jj in range( start=1, stop=size(block)[1] ) # block:(8,3) Mat
            vertx = block[jj,:] # (1,3)
            xbar  = vertx[1] -  ObsPoint[1]
            ybar  = vertx[2] -  ObsPoint[2]
            ζp = vertx[3] + ObsPoint[3]
            ζm = vertx[3] - ObsPoint[3]
            rp = sqrt(   xbar^2 + ybar^2 + ζp^2    )
            rm = sqrt(   xbar^2 + ybar^2 + ζm^2    )

            # minus sign comes from that Julia convention starting from 1, not 0! 
            Sσ = - [-1, 1, 1, -1, 1, -1, -1, 1]

            # -- Displacements --
            # same as (B. li et al., 2021)
            # Z downward positive 
            Disp[1] +=  Prefactor_disp*Sσ[jj]*(fFunc(ybar,ζm,xbar,rm) 
                                + (3. - 4*PossionRatio)*fFunc(ybar, ζp,  xbar, rp) 
                                + 2*ObsPoint[3]*log(abs(rp + ybar)))

            Disp[2] +=  Prefactor_disp*Sσ[jj]*(fFunc(xbar,ζm,ybar,rm)
                                + (3. - 4*PossionRatio)*fFunc(xbar, ζp,  ybar, rp)
                                + 2*ObsPoint[3]*log(abs(rp + xbar)))

            Disp[3] += - Prefactor_disp*Sσ[jj]*(fFunc(xbar,ybar,ζm,rm)
                                + (3. - 4*PossionRatio)*fFunc(xbar,ybar, ζp, rp) 
                                - 2 * ObsPoint[3] * atan((ζp*rp)/(xbar*ybar)) )

            # -- Stresses --
            # same as Flow2Quake and (B. li et al., 2021)
            # compression positive + Z downward positive
            
            sig_11 += Prefactor_sig*Sσ[jj]*(
                        - atan((xbar*rm)/(ybar*ζm))
                        - (3.0-4*PossionRatio)* atan((xbar*rp)/(ybar*ζp))
                        + 4*PossionRatio* atan((ζp*rp)/(xbar*ybar))
                        + 2*ObsPoint[3]*xbar*ybar/( xbar^2 +  ζp^2 )/rp)

            # σYY
            sig_22 += Prefactor_sig*Sσ[jj]*(
                        -  atan((ybar*rm)/(xbar*ζm))
                        - (3.0-4*PossionRatio)* atan((ybar*rp)/(xbar*ζp))
                        + 4*PossionRatio* atan((ζp*rp)/(xbar*ybar))
                        + 2*ObsPoint[3]*xbar*ybar/( ybar^2 + ζp^2 )/rp)

            # σZZ
            sig_33 +=  Prefactor_sig*Sσ[jj]*(
                            -  atan((ζm*rm)/(xbar*ybar))
                            +  atan((ζp*rp)/(xbar*ybar))
                            - 2*ObsPoint[3]*xbar*ybar/rp*(1/(xbar^2 +  ζp^2) + 1/(ybar^2 +  ζp^2)) )

            # σXY
            sig_12 += - Prefactor_sig*Sσ[jj]*( log(abs(rm + ζm))
                            + (3.0-4*PossionRatio)* log(abs(rp+ζp))
                            + 2*ObsPoint[3]/rp   )

            # σXZ 
            sig_13 +=  Prefactor_sig*Sσ[jj]*( log( abs((rm+ybar)/(rp+ybar)))
                            + 2*ObsPoint[3]*ybar*ζp/( xbar^2+ ζp^2)/rp)
            # σYZ
            sig_23 +=  Prefactor_sig*Sσ[jj]*( log( abs((rm+xbar)/(rp+xbar)))
                            + 2*ObsPoint[3]*xbar*ζp/(ybar^2+ ζp^2)/rp)
        end

    end
    
    # ----对比COMSOL结果，手动符号修正----
    Disp[3] *= -1
    sig_13 *= -1
    sig_23 *= -1
    # ------------------------------------
    Sigma = [sig_11, sig_22, sig_33, sig_12, sig_13, sig_23]
    
    return Disp, Sigma 
end


function Get_Cube_8Vertices_Pos(CubeOrigin, CubeLengthX, CubeLengthY, CubeThickness)
    Block = zeros(8,3)
    Block[1,:] = [CubeOrigin[1]-CubeLengthX/2,CubeOrigin[2]-CubeLengthY/2,CubeOrigin[3]-CubeThickness/2]
    Block[2,:] = [CubeOrigin[1]+CubeLengthX/2,CubeOrigin[2]-CubeLengthY/2,CubeOrigin[3]-CubeThickness/2]
    Block[3,:] = [CubeOrigin[1]-CubeLengthX/2,CubeOrigin[2]+CubeLengthY/2,CubeOrigin[3]-CubeThickness/2]
    Block[4,:] = [CubeOrigin[1]+CubeLengthX/2,CubeOrigin[2]+CubeLengthY/2,CubeOrigin[3]-CubeThickness/2]
    Block[5,:] = [CubeOrigin[1]-CubeLengthX/2,CubeOrigin[2]-CubeLengthY/2,CubeOrigin[3]+CubeThickness/2]
    Block[6,:] = [CubeOrigin[1]+CubeLengthX/2,CubeOrigin[2]-CubeLengthY/2,CubeOrigin[3]+CubeThickness/2]
    Block[7,:] = [CubeOrigin[1]-CubeLengthX/2,CubeOrigin[2]+CubeLengthY/2,CubeOrigin[3]+CubeThickness/2]
    Block[8,:] = [CubeOrigin[1]+CubeLengthX/2,CubeOrigin[2]+CubeLengthY/2,CubeOrigin[3]+CubeThickness/2]
    return Block
end


function BuildReservoirGeometry_even( CenterPoint, Width, Thickness, Num_for_EachSide = 4 )

    # Definition of "corner": Xmin, Ymin
    X_corner = CenterPoint[1] - ( 0.5 + (Num_for_EachSide/2 -1 ) )*Width
    Y_corner = CenterPoint[2] - ( 0.5 + (Num_for_EachSide/2 -1 ) )*Width
    CubeCenter_Z = CenterPoint[3] + Thickness/2

    BLOCKS_Center_pos_XYZ = Vector{Vector{Float64}}()
    for idx_y in range(0, Num_for_EachSide-1, step=1) 
        CubeCenter_Y = Y_corner + idx_y * Width        
        for idx_x in range(0, Num_for_EachSide-1, step=1)
            CubeCenter_X = X_corner + idx_x * Width 
            push!(BLOCKS_Center_pos_XYZ, [CubeCenter_X, CubeCenter_Y, CubeCenter_Z])
        end
    end
    return BLOCKS_Center_pos_XYZ
end