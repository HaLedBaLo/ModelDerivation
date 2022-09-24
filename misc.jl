

function remove_time_dependencies(eqm, ddq_, dq_, q_, ddq, dq, q)
    len_state = length(q_)
    for i = 1:len_state
        eqm = eqm.subs(ddq_[i], ddq[i])
        eqm = eqm.subs(dq_[i], dq[i])
        eqm = eqm.subs(q_[i], q[i])
    end
    return eqm
end