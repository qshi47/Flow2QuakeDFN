include("Functions_Kuvshinov_Cubes.jl")

ObsPoint = [2000.1, 1000.1, 0]

CubeOrigin = [2000, 1000, 2900]
CubeLengthX = 4000
CubeLengthY = 2000
CubeThickness = 200

block = Get_Cube_8Vertices_Pos(CubeOrigin, CubeLengthX, CubeLengthY, CubeThickness)

ShearModulus = 1
block_Cm = 1
PossionRatio = 0.25

Disp, Sigma = Calculate_PoroDispStress_SingleBlock_FullSpace_GF( ObsPoint, block, block_Cm, ShearModulus, PossionRatio )

