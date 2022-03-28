using LaTeXStrings

function print_nice(state::State)
    for field in fieldnames(typeof(state))
        println(field, ": ", getfield(state, field))
    end
    return nothing
end
export print_nice

function print_basis_state(state, fields)
    str = "\\left|"
    for (i, field) in enumerate(fields)
        val = getfield(state, field)
        if isinteger(val)
            val_str = string(val.num)
        else
            val_str = "\\frac{$(val.num)}{$(val.den)}"
        end
        str *= string(field) * "=" * val_str
        if i < length(fields)
            str *= string(", ")
        end
    end
    str *= "\\right\\rangle"
    return str
end

function Base.show(io::IO, m::MIME"text/latex", state::State)
    basis = state.basis
    coeffs = state.coeffs
    basis_type = typeof(basis[1])
    
    printed = false
    fields = fieldnames(basis_type)[2:end-1]
    
    for (i, coeff) in enumerate(coeffs)
        str = "\$\$"
        basis_state = basis[i]
        if norm(coeff) > 1e-3
            plus_sign = true
            state_str = ""
            real_val = string(round(real(coeff), digits=3))
            imag_val = string(round(imag(coeff), digits=3))
            if (abs(real(coeff)) > 1e-3) && (abs(imag(coeff)) < 1e-5)
                state_str *= real_val
                if real(coeff) < 0
                    plus_sign = false
                end
            elseif (abs(imag(coeff)) > 1e-3) && abs((real(coeff)) < 1e-5)
                state_str *= imag_val * "i"
                if imag(coeff) < 0
                    plus_sign = false
                end
            else
                state_str *= "(" * real_val * " + " * imag_val * "i" * ")"
            end
            if plus_sign && printed
                str *= " + "
            end
            str *= state_str
            str *= print_basis_state(basis_state, fields)
            str *= "\$\$"
            latex_str = latexstring(str)
            println(io, latex_str)
            printed = true
        end
    end
    return nothing
end

function print_matrix(A)
    println("Real part:")
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            A_ij = real(A[i,j])
            if abs(A_ij) < 1e-5
                A_ij = convert(Int64, floor(A_ij))
            else
                A_ij = round(A_ij, digits=2)
            end
            print(A_ij, " ")
        end
        println()
    end
    println()
    println("Imaginary part:")
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            A_ij = imag(A[i,j])
            if abs(A_ij) < 1e-5
                A_ij = convert(Int64, floor(A_ij))
            else
                A_ij = round(A_ij, digits=2)
            end
            print(A_ij, " ")
        end
        println()
    end
end
export print_matrix