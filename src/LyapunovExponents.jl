"""
    A module to calculate lyapunov exponents.
"""
module LyapunovExponents

using ForwardDiff: _lyap_div!
using LinearAlgebra: AbstractMatrix
using DifferentialEquations 
using LinearAlgebra 
using DocStringExtensions 
using ForwardDiff

include("dynamics.jl") 
include("lyapunovs.jl")


end # module 
