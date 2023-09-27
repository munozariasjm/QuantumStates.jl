using HalfIntegers

import Base.factorial
factorial(x::HalfInt64) = factorial(div(x.twice, 2))
export factorial

Δ(a,b,c) = factorial(a+b-c) * factorial(a-b+c) * factorial(-a+b+c) / factorial(a+b+c+1)
export Δ

@inline x(t,j1,j2,j3,m1,m2) =
    factorial(t) * factorial(t+j3-j2+m1) * factorial(t+j3-j1-m2) * factorial(j1+j2-j3-t) * factorial(j1-m1-t) * factorial(j2+m2-t)

function wigner3j(j1,j2,j3,m1,m2,m3)
    threeJ_value = 0.0
    
    if (-j1 <= m1 <= j1) && (-j2 <= m2 <= j2) && (-j3 <= m3 <= j3) && (m1 + m2 == -m3) && (abs(j1 - j2) <= j3 <= j1 + j2)
        min_t = max(0, j2-j3-m1, j1-j3+m2)
        max_t = min(j1+j2-j3, j1-m1, j2+m2)
        
        if min_t <= max_t
            threeJ_value += (
                (-1)^(j1-j2-m3) * sqrt(Δ(j1,j2,j3)) 
                * sqrt(factorial(j1+m1)) * sqrt(factorial(j1-m1)) * sqrt(factorial(j2+m2)) * sqrt(factorial(j2-m2)) * sqrt(factorial(j3+m3)) * sqrt(factorial(j3-m3))
                * sum((-1)^t / x(t,j1,j2,j3,m1,m2) for t in min_t:max_t)
            )
        end
    end
    return threeJ_value
end
export wigner3j

@inline f(t::HalfInteger,
    j1j2j3::HalfInteger,
    j1J2J3::HalfInteger,
    J1j2J3::HalfInteger,
    J1J2j3::HalfInteger,
    j1j2J1J2::HalfInteger,
    j2j3J2J3::HalfInteger,
    j3j1J3J1::HalfInteger
    ) = factorial(t-j1j2j3) * factorial(t-j1J2J3) * factorial(t-J1j2J3) * factorial(t-J1J2j3) * factorial(j1j2J1J2-t) * factorial(j2j3J2J3-t) * factorial(j3j1J3J1-t)

@inline f(t,
    j1j2j3,
    j1J2J3,
    J1j2J3,
    J1J2j3,
    j1j2J1J2,
    j2j3J2J3,
    j3j1J3J1
    ) = factorial(t-j1j2j3) * factorial(t-j1J2J3) * factorial(t-J1j2J3) * factorial(t-J1J2j3) * factorial(j1j2J1J2-t) * factorial(j2j3J2J3-t) * factorial(j3j1J3J1-t)

function wigner6j(j1,j2,j3,J1,J2,J3)
    sixJ_value = 0.0
    
    if (abs(j1-j2) <= j3 <= j1+j2) && (abs(j1-J2) <= J3 <= j1+J2) && (abs(J1-j2) <= J3 <= J1+j2) && (abs(J1-J2) <= j3 <= J1+J2)
        
        j1j2j3 = j1+j2+j3
        j1J2J3 = j1+J2+J3
        J1j2J3 = J1+j2+J3
        J1J2j3 = J1+J2+j3

        j1j2J1J2 = j1+j2+J1+J2
        j2j3J2J3 = j2+j3+J2+J3
        j3j1J3J1 = j3+j1+J3+J1

        min_t = max(j1j2j3, j1J2J3, J1j2J3, J1J2j3)
        max_t = min(j1j2J1J2, j2j3J2J3, j3j1J3J1)

        if min_t <= max_t
            sixJ_value += (
                sqrt(Δ(j1,j2,j3) * Δ(j1,J2,J3) * Δ(J1,j2,J3) * Δ(J1,J2,j3))
                * sum( (-1)^t * factorial(t+1) / f(t,j1j2j3,j1J2J3,J1j2J3,J1J2j3,j1j2J1J2,j2j3J2J3,j3j1J3J1) for t in min_t:max_t)
                )
        end
    end
    return sixJ_value
        
end
export wigner6j
    
function wigner9j(j1, j2, j3, j4, j5, j6, j7, j8, j9)::Float64
    val = 0.0
    kmin = max(abs(j1 - j9), abs(j4 - j8), abs(j2 - j6))
    kmax = min(abs(j1 + j9), abs(j4 + j8), abs(j2 + j6))
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
        WignerSymbols.wigner3j(Float64, j1, j2, j3, m1, m2, m3) 
    catch 
        0.0
    end
end

function wigner6j_(j1, j2, j3, m1, m2, m3)
    try 
        WignerSymbols.wigner6j(Float64, j1, j2, j3, m1, m2, m3) 
    catch 
        0.0
    end
end
export wigner6j_

# struct QuantumNumber
#     val::Int64
#     is_half::Bool
# end

# import Base: isequal, <=, +, -
# isequal(x::QuantumNumber, y::QuantumNumber) = (x.val == y.val) && (x.is_half == y.is_half)
# <=(x::QuantumNumber, y::QuantumNumber) = (x.val <= y.val)
# function +(x::QuantumNumber, y::QuantumNumber)
#     double = x.is_half && y.is_half
#     is_half = ~(x.is_half == y.is_half)
#     val = (1 + ~(x.is_half)) * x.val + (1 + ~(y.is_half)) * y.val
#     QuantumNumber(div(val, 2 - double) * (~is_half) + val * is_half, is_half)
# end
