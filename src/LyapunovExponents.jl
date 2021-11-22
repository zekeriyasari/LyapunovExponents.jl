"""
    A module to calculate lyapunov exponents.
"""
module LyapunovExponents

using DifferentialEquations 
using LinearAlgebra 
using DocStringExtensions 
using ForwardDiff

include("dynamics.jl") 
include("lyapunovs.jl")
include("network.jl")

end # module 
