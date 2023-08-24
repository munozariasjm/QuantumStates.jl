mutable struct Transition
    ground_state::State
    excited_state::State
    frequency::Float64
    tdm::ComplexF64
end
export Transition

function compute_transitions(H::Hamiltonian, p, threshold=1e-8)
    transitions = Transition[]
    
    for (i, basis_state) ∈ enumerate(H.basis)
        for (j, basis_state′) ∈ enumerate(H.basis)
            H.tdms[i,j,p] = TDM(basis_state, basis_state′, p)
        end
    end
            
    for (i, state) ∈ enumerate(H.states)
        for (j, state′) ∈ enumerate(H.states)
            if state′.E > state.E
                tdm = state.coeffs ⋅ (H.tdms * state′.coeffs)
                if norm(tdm) > threshold
                    transition = Transition(state, state′, state′.E - state.E, tdm)
                    push!(transitions, transition)
                    H.tdms
                end
            end
        end
    end
    return transitions
end

"""
    Compute the transitions between two sets of states, `states` and `states′`.
"""
function compute_transitions(states::Vector{<:State}, states′::Vector{<:State}, p; threshold=1e-8, compute_tdms=true)
    transitions = Transition[]
    basis = states[1].basis

    tdms = zeros(length(basis), length(basis))
    if compute_tdms
        for (i, basis_state) ∈ enumerate(basis)
            for (j, basis_state′) ∈ enumerate(basis)
                tdms[i,j] = TDM(basis_state, basis_state′, p)
            end
        end
    end

    for state ∈ states
        for state′ ∈ states′
            if state′.E > state.E
                tdm = state.coeffs ⋅ (tdms * state′.coeffs)
                f = state′.E - state.E
                if (norm(tdm) > threshold || ~compute_tdms) && abs(f) > 1
                    transition = Transition(state, state′, f, tdm)
                    push!(transitions, transition)
                end
            end
        end
    end
    return transitions
end
export compute_transitions

ground_state(transition::Transition) = transition.ground_state
excited_state(transition::Transition) = transition.excited_state
frequency(transition::Transition) = transition.frequency
tdm(transition::Transition) = transition.tdm
export ground_state
export excited_state

function transitions_table(transitions::Vector{Transition}, relabelling_states=nothing::Union{Nothing, Vector{<:State}}; threshold=1e-8)
    
    ground_states = ground_state.(transitions)
    excited_states = excited_state.(transitions)
    frequencies = frequency.(transitions)
    tdms = tdm.(transitions)
    
    if ~isnothing(relabelling_states)
        relabelled_ground_states = relabelling_states[getfield.(ground_states, :idx)]
        relabelled_excited_states = relabelling_states[getfield.(excited_states, :idx)]
    else
        relabelled_ground_states = ground_states
        relabelled_excited_states = excited_states
    end
    
    df1 = states_table(relabelled_ground_states, threshold=1e-1, dominant_state_only=true)
    df2 = states_table(relabelled_excited_states, threshold=1e-1, dominant_state_only=true)
    
    df_transitions = hcat(df1, df2, makeunique=true)
    df_transitions[!, :f] = frequencies
    df_transitions[!, :tdm] = real.(tdms)
    sort!(df_transitions, :f)
    
    return df_transitions
end
export transitions_table

function calculate_state_overlaps!(states, states′, overlaps)
    for i ∈ eachindex(states)
        for j ∈ eachindex(states′)
            state = states[i]
            state′ = states′[j]
            overlaps[i,j] = norm(state ⋅ state′)
        end
    end
    return nothing
end  