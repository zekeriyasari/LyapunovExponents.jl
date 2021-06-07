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

export Lorenz, Chua, Chen, Rossler, HRNeuron, lyapunovs

end # module 
