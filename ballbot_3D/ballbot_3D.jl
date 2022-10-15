using SymPy

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
I_w = diagm([I_wi; I_wi; I_wi])
grav = [0; 0; g]

# kinematics that map the rotation of the wheels to the corresponding rotation
# in x, y, z between ball and torso
for_kin, inv_kin = omniwheel_kinematics(
    sympy.pi / 4, sympy.S.Zero, sympy.pi * 2 / 3)
for_kin *= r_w / r_b
inv_kin *= r_b / r_w

rot_global_to_torso = rot_3d(q_[3], q_[4], q_[5])

# Positions
p_ball = [q_[1:2]; r_b]
p_torso = p_ball + rot_global_to_torso * l

# Velocities
lin_vel_ball = p_ball.diff(t)
rot_vel_ball = lin_vel_ball / r_b

lin_vel_torso = p_torso.diff(t)
rot_vel_torso = dq_[3:5]

rot_vel_wheel = inv_kin * (rot_global_to_torso * rot_vel_ball - rot_vel_torso)

# Potential energies
V_ball = m_b * grav.T * p_ball
V_torso = m_t * grav.T * p_torso

# Kinetic energies
T_ball = lin_vel_ball.T * m_b / 2 * lin_vel_ball +
         rot_vel_ball.T * I_b / 2 * rot_vel_ball
T_torso = lin_vel_torso.T * m_t / 2 * lin_vel_torso +
          rot_vel_torso.T * I_t / 2 * rot_vel_torso
T_wheel = rot_vel_wheel.T * I_w / 2 * rot_vel_wheel

# Lagrangian
lagrangian = T_ball[] + T_torso[] + T_wheel[] - V_torso[] - V_ball[]
eqm = calculate_equationofmotion_from_lagrangian(lagrangian, q_)
eqm = remove_time_dependencies(eqm, [ddq_; dq_; q_], [ddq; dq; q])

# Mass Matrix
M, eqm_tmp = collect_variable(eqm, ddq, [5, 5])

# Corioli Matrix
C = sympy.zeros(5, 5);
for i = 1:5
    for j = 1:5
        collection = collect(eqm_tmp[i], dq[j]^2)
        if (collection.coeff(dq[j]^2) != 0)
            C[i, j] = collection.coeff(dq[j]^2) * dq[j]
            eqm_tmp[i] = expand(eqm_tmp[i] - C[i, j] * dq[j])
        end
    end
end

# Gravity vector
G = eqm_tmp

# Verify Matrix calculation
if (simplify.(eqm - M * ddq - C * dq - G) != zeros(5, 1))
    print("equation of motion did not resolve to zero")
end

# Motor torques
rot_vel_wheel = remove_time_dependencies(rot_vel_wheel, [dq_ q_], [dq q])
J_wheel = sympy.zeros(5, 3);
for i = 1:5
    for j = 1:3
        J_wheel[i, j] = rot_vel_wheel[j].diff(dq[i])
    end
end

# Export variables
export_function(M, "bb3d_M", [q; dq])
export_function(C, "bb3d_C", [q; dq])
export_function(G, "bb3d_G", [q; dq])
export_function(J_wheel, "bb3d_Q", [q; dq])
