using SymPy

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

    # get arguments of function and common subexpressions for faster execution
    arguments = free_symbols(vec(variable))
    subexpressions, variable = sympy.cse(variable)

    io = open(filename, "w")
    write(io, string("function ", name, "("))
    write_vector(arguments, io)
    write(io, ")\n")
    if (~isempty(subexpressions))
        write_tuple(subexpressions, io)
    end
    write(io, "\treturn ")
    write(io, string(variable[1])[4:end])
    write(io, "\nend")
    close(io)
end

function write_vector(vector, io)
    last_index = lastindex(vector)
    for i = 1:last_index
        if (i == last_index)
            write(io, string(vector[i]))
        else
            write(io, string(vector[i], ", "))
        end
    end
end

function write_tuple(tuple, io)
    for i = 1:lastindex(tuple)
        write(io, string("  ", tuple[i][1], " = ", tuple[i][2], "\n"))
    end
end