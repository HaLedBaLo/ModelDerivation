using SymPy
using LinearAlgebra

include("../misc.jl")
include("../rotations.jl")
include("kinematic.jl")

## Symbolic variable declaration
# state 
@syms x, dx, ddx
@syms y, dy, ddy
@syms ex, dex, ddex
@syms ey, dey, ddey
@syms ez, dez, ddez
@syms tau1, tau2, tau3
# system parameters 
@syms t
@syms r_b, r_w
@syms m_b, m_t
@syms I_bi, I_tx, I_ty, I_tz, I_wi
@syms lx, ly, lz
@syms g
# time dependecies
x_ = SymFunction("x")(t)
y_ = SymFunction("y")(t)
ex_ = SymFunction("ex")(t)
ey_ = SymFunction("ey")(t)
ez_ = SymFunction("ez")(t)

## Vector assembly
# state
q = [x; y; ex; ey; ez]
dq = [dx; dy; dex; dey; dez]
ddq = [ddx; ddy; ddex; ddey; ddez]
# state time_dependent
q_ = [x_; y_; ex_; ey_; ez_]
dq_ = q_.diff(t)
ddq_ = dq_.diff(t)
# system parameter
l = [0; 0; lz]
I_b = diagm([I_bi; I_bi; I_bi])
I_t = diagm([I_tx; I_tx; I_tz])

for_kin, inv_kin = omniwheel_kinematics(
    sympy.pi / 4, sympy.S.Zero, sympy.pi * 2 / 3)

## positions
# ball position
p_b = [q_[1:2]; r_b]
# torso position
p_t = p_b + rot_3d(q_[3], q_[4], q_[5]) * l

## velocities
# ball velocity
lin_vel_b = p_b.diff(t)
rot_vel_b = lin_vel_b / r_b
# torso velocity
lin_vel_t = p_t.diff(t)
rot_vel_t = dq_[3:5]
# motor velocity
