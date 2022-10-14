using SymPy
using LinearAlgebra

include("../rotations.jl")

function omniwheel_kinematics(angle_x, angle_y, angle_z)
    v = [0; 0; -1]

    v1 = rot_3d(angle_x, angle_y, 0) * v
    v2 = rot_3d(angle_x, angle_y, angle_z) * v
    v3 = rot_3d(angle_x, angle_y, -angle_z) * v

    inv_kin = transpose([v1 v2 v3])
    for_kin = inv(inv_kin)

    return for_kin, inv_kin
end