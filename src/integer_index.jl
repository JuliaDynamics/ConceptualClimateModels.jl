# TODO: this is necessary because of requiring the parameters to
# be changed in the optimization library. However, is it ever used anywhere...?

export integer_parameter_index, integer_state_index

"""
    integer_parameter_index(symbol, ds::DynamicalSystem) → i::Int

Convert the given symbolic variable representing a parameter to its integer
index in the parameter container ([`current_parameters`](@ref)).
Return `nothing` if `ds` doesn't reference a symbolic model
or the model does not have the given `symbol` as parameter.
"""
integer_parameter_index(s, ds::DynamicalSystem) = integer_parameter_index(s, referrenced_sciml_model(ds))
function integer_parameter_index(symbol, model)
    isnothing(model) && return nothing
    findfirst(isequal(symbol), parameters(model))
end

"""
    integer_state_index(symbol, ds::DynamicalSystem) → i::Int

Convert the given symbolic variable representing a state (dynamic) variable to its integer
index in the state vector ([`current_state`](@ref)).
Return `nothing` if `ds` doesn't reference a symbolic model
or the model does not have the given `symbol` as a state variable.
"""
integer_state_index(s, ds::DynamicalSystem) = integer_state_index(s, referrenced_sciml_model(ds))
function integer_state_index(symbol, model)
    isnothing(model) && return nothing
    findfirst(isequal(symbol), states(model))
end
