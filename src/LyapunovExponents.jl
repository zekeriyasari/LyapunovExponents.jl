"""
    A module to calculate lyapunov exponents.
"""
module LyapunovExponents

using Base: entry_point_and_project_file
using LinearAlgebra: AbstractMatrix
using DifferentialEquations 
using LinearAlgebra 
using DocStringExtensions 
using ForwardDiff

include("dynamics.jl") 
include("lyapunovs.jl")
include("network.jl")

end # module 
