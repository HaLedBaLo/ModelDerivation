function calculate_equationofmotion_from_lagrangian(lagrangian, states)
    len_state = length(states)
    dL_dq = sympy.zeros(len_state, 1)
    dL_ddq = sympy.zeros(len_state, 1)
    for i = 1:len_state
        dL_dq[i] = lagrangian.diff(q_[i])
        dL_ddq[i] = lagrangian.diff(dq_[i])
    end
    eqm = dL_ddq.diff(t) - dL_dq
    return eqm
end

function replace_variables(var, old, new)
    for i = 1:lastindex(old)
        var = var.subs(old[i], new[i])
    end
    return var
end

function collect_variable(input::Array, vars, dims)
    collection = sympy.zeros(dims[1], dims[2])
    rest = expand.(input)

    for i = 1:dims[1]
        for j = 1:dims[2]
            col = collect(rest[i], vars[j])
            coeff = col.coeff(vars[j])
            if (coeff != 0)
                collection[i, j] = coeff
                rest[i] -= expand(coeff * vars[j])
            end
        end
    end
    return collection, rest
end

function collect_variable(input::Sym, vars, dims)
    collection = sympy.zeros(dims[1], dims[2])
    rest = expand(input)
    for i = 1:maximum(dims)
        col = collect(rest, vars[i])
        coeff = col.coeff(vars[i])
        if (coeff != 0)
            collection[i] = coeff
            rest -= expand(coeff * vars[i])
        end
    end
    return collection, rest
end

function collect_corioli_matrix(input::Array, vars, dims)
    corioli = sympy.zeros(dims[1], dims[2])
    rest = expand.(input)
    for i = 1:dims[1]
        for j = 1:dims[2]
            col = collect(rest[i], vars[j]^2)
            coeff = col.coeff(vars[j]^2)
            if (coeff != 0)
                corioli[i, j] = coeff * vars[j]
                rest[i] -= expand(coeff * vars[j]^2)
            end
        end
        for j = 1:dims[2]
            col = collect(rest[i], vars[j])
            coeff = col.coeff(vars[j])
            if (coeff != 0)
                corioli[i, j] += coeff
                rest[i] -= expand(coeff * vars[j])
            end
        end
    end
    return corioli, rest
end

function check_variable_existence(var1, var2)
    exists = false
    for i = 1:lastindex(var1)
        if (var2 in var1[i].free_symbols)
            exists = true
            break
        end
    end
    return exists
end

function check_variable_existence(var1, var2::Vector)
    exists = false
    for i = 1:lastindex(var1)
        for j = 1:lastindex(var2)
            if (var2[j] in var1[i].free_symbols)
                exists = true
                break
            end
        end
    end
    return exists
end

function export_variable(variable, name::String)
    filename = string("export/", name, ".jl")
    filecontent = string(name, " = ", variable)
    if (!isdir("export"))
        mkdir("export")
    end
    io = open(filename, "w")
    write(io, filecontent)
    close(io)
end

function export_function(variable, name::String)
    filename = string("export/", name, ".jl")
    if (!isdir("export"))
        mkdir("export")
    end
    io = open(filename, "w")
    write(io, "function ")
    write(io, name)
    write(io, "()\n\treturn ")
    write(io, string(variable)[4:end])
    write(io, "\nend")
    close(io)
end

function export_function(variable, name::String, argin::Sym)
    filename = string("export/", name, ".jl")
    if (!isdir("export"))
        mkdir("export")
    end
    io = open(filename, "w")
    write(io, "function ")
    write(io, name)
    write(io, "(")
    write(io, string(argin))
    write(io, ")\n\treturn ")
    write(io, string(variable)[4:end])
    write(io, "\nend")
    close(io)
end

function export_function(variable, name::String, argin::Vector{Sym})
    filename = string("export/", name, ".jl")
    if (!isdir("export"))
        mkdir("export")
    end
    io = open(filename, "w")
    write(io, "function ")
    write(io, name)
    write(io, "(")
    n_argin = lastindex(argin)
    argin_found = falses(n_argin, 1)
    for i = 1:n_argin
        if (check_variable_existence(variable, argin[i]))
            argin_found[i] = true
        end
    end
    last_index = get_last_true_index(argin_found)
    if (last_index > 0)
        for i = 1:n_argin
            if argin_found[i]
                write(io, string(argin[i]))
                if (i < last_index)
                    write(io, ", ")
                end
            end
        end
    end
    write(io, ")\n\treturn ")
    write(io, string(variable)[4:end])
    write(io, "\nend")
    close(io)
end

function get_last_true_index(vector)
    index = 0
    for i = 1:lastindex(vector)
        if (vector[i])
            index = i
        end
    end
    return index
end
