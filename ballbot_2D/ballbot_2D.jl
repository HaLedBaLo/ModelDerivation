using SymPy

include("../misc.jl")

# time
@syms t
# state variables
@syms x, dx, ddx
@syms alpha, dalpha, ddalpha
@syms tau
# model constants
@syms r_b, r_w
@syms m_b, m_t
@syms I_b, I_t, I_w
@syms l
@syms g

# time dependecies
x_ = SymFunction("x")(t)
alpha_ = SymFunction("alpha")(t)

# state
q = [x; alpha]
dq = [dx; dalpha]
ddq = [ddx; ddalpha]

q_ = [x_; alpha_]
dq_ = q_.diff(t)
ddq_ = dq_.diff(t)

## positions
# ball position
p_b = [x_; 0]
# torso position
p_t = [x_ - l * sin(alpha_); l * cos(alpha_)]

## velocities 
# ball velocity
v_b = p_b.diff(t)
omega_b = v_b / r_b
# torso velocity
v_t = p_t.diff(t)
omega_t = diff.(alpha_, t)
# motor velocity in terms of state
omega_m = r_b / r_w * (omega_t - omega_b[1])

## Energies
# Potential energies
V_t = m_t * g * p_t[2]
# Kinetic energy
T_b = v_b.T * m_b / 2 * v_b + omega_b.T * I_b / 2 * omega_b
T_t = v_t.T * m_t / 2 * v_t + [omega_t * I_t / 2 * omega_t]
T_w = I_w / 2 * omega_m^2;

## Lagrangian 
L = T_b[] + T_t[] + T_w[] - V_t[]

## Lagrange's equations
eqm = calculate_equationofmotion_from_lagrangian(L, q_)
eqm = remove_time_dependencies(eqm, [ddq_; dq_; q_], [ddq; dq; q])

## Mass matrix
M, eqm_tmp = collect_variable(eqm, ddq, [2, 2])

## Corriolis matrix
C = sympy.zeros(2, 2);

for i = 1:2
    for j = 1:2
        collection = collect(eqm_tmp[i], dq[j]^2)
        if (collection.coeff(dq[j]^2) != 0)
            C[i, j] = collection.coeff(dq[j]^2) * dq[j]
            eqm_tmp[i] = expand(eqm_tmp[i] - C[i, j] * dq[j])
        end
    end
end

## Gravity vector
G = eqm_tmp

## Verify Matrix calculation
if (simplify.(eqm - M * ddq - C * dq - G) != zeros(2, 1))
    print("equation of motion did not resolve to zero")
end

## Generalized forces
omega_m = remove_time_dependencies(omega_m, dq_, dq)
J_m, rest = collect_variable(omega_m, dq, [1, 2])

Q = J_m.transpose()

## Export variables
export_function(M, "bb2d_M", [q; dq])
export_function(C, "bb2d_C", [q; dq])
export_function(G, "bb2d_G", [q; dq])
export_function(Q, "bb2d_Q")
