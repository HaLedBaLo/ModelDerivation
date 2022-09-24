

function remove_time_dependencies(var, old, new)
    for i = 1:lastindex(old)
        var = var.subs(old[i], new[i])
    end
    return var
end

function collect_variable(input::Array, variable, dims)
    collection = sympy.zeros(dims[1], dims[2])
    rest = expand.(input)

    for i = 1:dims[1]
        for j = 1:dims[2]
            col = collect(rest[i], variable[j])
            if (col.coeff(variable[j]) != 0)
                collection[i,j] = col.coeff(variable[j])
                rest[i] -= expand(collection[i,j] * variable[j])
            end
        end
    end
    return collection, rest
end

function collect_variable(input, variable, dims)
    collection = sympy.zeros(dims[1], dims[2])
    rest = expand(input)

    for i = 1:maximum(dims)
        col = collect(rest, variable[i])
        if (col.coeff(variable[i]) != 0)
            collection[i] = col.coeff(variable[i])
            rest -= expand(collection[i] * variable[i])
        end
    end
    return collection, rest
end