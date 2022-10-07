using SymPy

include("../rotations.jl")

v = [0; 0; 1]

angle_x = pi / 4
angle_y = pi / 12
angle_z = 2 / 3 * pi      # 120 deg

vx = rot_3d(angle_x, angle_y, 0) * v
vy = rot_3d(angle_x, angle_y, angle_z) * v
vz = rot_3d(angle_x, angle_y, -angle_z) * v

kin = [vx vy vz]