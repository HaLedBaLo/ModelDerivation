using SymPy

include("misc.jl")

# time
@syms t
# state variables
@syms x, dx, ddx 
@syms alpha, dalpha, ddalpha
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
omega_t = diff.(alpha_,t)
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
len_state = length(q_)
dL_dq = sympy.zeros(len_state,1)
dL_ddq = sympy.zeros(len_state,1)
for i = 1:len_state
    dL_dq[i] = L.diff(q_[i])
    dL_ddq[i] = L.diff(dq_[i])
end
eqm = dL_ddq.diff(t) - dL_dq
# substitute time dependent variables
eqm = remove_time_dependencies(eqm, ddq_, dq_, q_, ddq, dq, q)

## Mass matrix
M = sympy.zeros(2,2)

eqm_tmp = expand.(eqm)
for i = 1:2
    for j = 1:2
        collection = collect(eqm_tmp[i], ddq[j])
        if (collection.coeff(ddq[j]) != 0)
            M[i,j] = collection.coeff(ddq[j])
            eqm_tmp[i] = expand(eqm_tmp[i] - M[i,j]*ddq[j])
        end
    end
end

## Corriolis matrix
C = sympy.zeros(2,2);

for i = 1:2
    for j = 1:2
        collection = collect(eqm_tmp[i], dq[j]^2)
        if (collection.coeff(dq[j]^2) != 0)
            C[i,j] = collection.coeff(dq[j]^2) * dq[j]
            eqm_tmp[i] = expand(eqm_tmp[i] - C[i,j] * dq[j])
        end
    end
end
  
## Gravity vector
G = eqm_tmp

## Verify Matrix calculation
simplify.(eqm - M * ddq - C * dq - G)

## Generalized forces
# omega_m