"""
    A module to calculate lyapunov exponents.
"""
module LyapunovExponents

using LinearAlgebra: AbstractMatrix
using DifferentialEquations 
using LinearAlgebra 
using DocStringExtensions 
using ForwardDiff

include("dynamics.jl") 
include("lyapunovs.jl")
include("msf.jl")


end # module 
