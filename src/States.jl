import WignerSymbols: wigner3j, wigner6j
using CompositeStructs

abstract type BasisState end
export BasisState

Base.@kwdef struct QuantumState <: BasisState
    E::Float64 = 0.0
end
    
struct State
    basis::Vector{BasisState}
    coeffs::Vector{ComplexF64}
end
export State

function wigner9j_(j1, j2, j3, j4, j5, j6, j7, j8, j9)
    val = 0.0
    kmin = maximum([abs(j1 - j9), abs(j4 - j8), abs(j2 - j6)])
    kmax = minimum([abs(j1 + j9), abs(j4 + j8), abs(j2 + j6)])
    if kmax >= kmin
        val += sum(
            (-1)^(2k) * (2k + 1) * 
            wigner6j(j1, j4, j7, j8, j9, k) * 
            wigner6j(j2, j5, j8, j4, k, j6) *
            wigner6j(j3, j6, j9, k, j1, j2) for k in kmin:kmax)
    end
    return val
end

function wigner3j_(j1, j2, j3, m1, m2, m3)
    try 
        wigner3j(j1, j2, j3, m1, m2, m3) 
    catch 
        0.0
    end
end

function wigner6j_(j1, j2, j3, m1, m2, m3)
    try 
        wigner6j(j1, j2, j3, m1, m2, m3) 
    catch 
        0.0
    end
end

function order(basis::Vector{<:BasisState}, ordering_QN)
    all_QN = [getfield(state, ordering_QN) for state in basis]
    idxs_sorted = sortperm(all_QN)
    return basis[idxs_sorted]
end
export order