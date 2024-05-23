Base.@kwdef struct VibrationalState <: BasisState
    E = 0.0
    v1::Int64
    v2::Int64
    v3::Int64
    constraints = ()
end
export VibrationalState

function I1(state, state′)
    n  = state.v1
    n′ = state′.v1
    return n == n′
end
export I1

function I2(state, state′)
    n  = state.v2
    n′ = state′.v2
    return n == n′
end
export I2

function I3(state, state′)
    n  = state.v3
    n′ = state′.v3
    return n == n′
end
export I3

function x1(state, state′)
    n  = state.v1
    n′ = state′.v1
    return (sqrt(n′ + 1) * (n == (n′ + 1)) + sqrt(n′) * (n == (n′ - 1))) * (state.v2 == state′.v2) * (state.v3 == state′.v3) #* sqrt( ħ / (2*m1*ω1) )
end
export x1

function x2(state, state′)
    n  = state.v2
    n′ = state′.v2
    return (sqrt(n′ + 1) * (n == (n′ + 1)) + sqrt(n′) * (n == (n′ - 1))) * (state.v1 == state′.v1) * (state.v3 == state′.v3) #* sqrt( ħ / (2*m2*ω2) )
end
export x2

function x3(state, state′)
    n  = state.v3
    n′ = state′.v3
    return (sqrt(n′ + 1) * (n == (n′ + 1)) + sqrt(n′) * (n == (n′ - 1))) * (state.v1 == state′.v1) * (state.v2 == state′.v2) #* sqrt( ħ / (2*m3*ω3) )
end
export x3

function p1(state, state′)
    n  = state.v1
    n′ = state′.v1
    return im * (sqrt(n′ + 1) * (n == (n′ + 1)) - sqrt(n′) * (n == (n′ - 1))) * (state.v2 == state′.v2) * (state.v3 == state′.v3) #* sqrt( ħ*m*ω / 2 )
end
export p1

function p2(state, state′)
    n  = state.v2
    n′ = state′.v2
    return im * (sqrt(n′ + 1) * (n == (n′ + 1)) - sqrt(n′) * (n == (n′ - 1))) * (state.v1 == state′.v1) * (state.v3 == state′.v3) #* sqrt( ħ*m*ω / 2 )
end
export p2

function p3(state, state′)
    n  = state.v3
    n′ = state′.v3
    return im * (sqrt(n′ + 1) * (n == (n′ + 1)) - sqrt(n′) * (n == (n′ - 1))) * (state.v1 == state′.v1) * (state.v2 == state′.v2) #* sqrt( ħ*m*ω / 2 )
end
export p3