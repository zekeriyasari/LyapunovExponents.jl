module LyapunovExponents

# TODO: Define jabobian functions for other dynamics.

using DifferentialEquations 
using LinearAlgebra 
using DocStringExtensions 
using ForwardDiff

const SOLVER = Tsit5()

##### Define dynamics 

abstract type Dynamics end

Base.@kwdef struct Lorenz <: Dynamics 
    σ::Float64 = 10 
    β::Float64 = 8 / 3 
    ρ::Float64 = 28 
    t::Float64 = 0.
    x0::Vector{Float64} = rand(3) 
end 

function (ds::Lorenz)(dx::AbstractVector, x::AbstractVector, u, t)  # Diffeq 
    dx[1] = ds.σ * (x[2] - x[1]) 
    dx[2] = x[1] * (ds.ρ - x[3]) - x[2] 
    dx[3] = x[1] * x[2] - ds.β * x[3]    
end

function (ds::Lorenz)(J::AbstractMatrix, x::AbstractVector, u, t)  # Jacobian 
    J .= [ 
        -ds.σ           ds.σ    0; 
        ds.ρ - x[3]     -1      -x[1]; 
        x[2]            x[1]    -β
    ]
end


Base.@kwdef struct Rossler <: Dynamics 
    α::Float64 = 0.2 
    β::Float64 = 0.2 
    γ::Float64 = 9
    t::Float64 = 0.
    x0::Vector{Float64} = rand(3) 
end 

function (ds::Rossler)(dx, x, u, t)
    dx[1] = -x[2] - x[3] 
    dx[2] = x[1] + ds.α * x[2] 
    dx[3] = ds.β + (x[1] - ds.γ) * x[3]    
end

Base.@kwdef struct Chen <: Dynamics 
    a::Float64 = 35 
    b::Float64 = 28  
    c::Float64 = 8/3 
    t::Float64 = 0.
    x0::Vector{Float64} = rand(3) 
end 

function (ds::Chen)(dx, x, u, t)
    dx[1] = ds.a * (x[2] - x[1]) 
    dx[2] = x[1] * (ds.c - ds.a - x[3]) + ds.c * x[2] 
    dx[3] = x[1] * x[2] - ds.β * x[3] 
end

Base.@kwdef struct Diode 
    a::Float64 = -1.27 
    b::Float64 = -0.68
end 
function (diode::Diode)(x) 
    if x > 1 
        -ds.b * x - ds.a + ds.b 
    elseif abs(x) ≤ 1 
        -ds.a * x 
    else
        -ds.b * x + ds.a - ds.b    
    end  
end 


Base.@kwdef struct Chua <: Dynamics 
    diode::Diode = Diode() 
    α::Float64 = 10. 
    β::Float64 = 14.87 
    t::Float64 = 0.
    x0::Vector{Float64} = rand(3)
end 

function (ds::Chua)(dx, x, u, t) 
    dx[1] = ds.α * (x[2] - x[1] + ds.diode(x[1]))
    dx[2] = x[1] - x[2] + x[3] 
    dx[3] = -ds.β * x[2] 
end 

Base.@kwdef struct HRNeuron 
    r::Float64 = 0.006 
    s::Float64 = 4 
    I::Float64 = 3.2 
    t::Float64 = 0.
    x0::Vector{Float64} = rand(3)
end 
function (ds::HRNeuron)(dx, x, u, t) 
    dx[1] = x[2] + 3 * x[1]^2 - x[1]^3 - x[3] + ds.I 
    dx[2] = 1 - 5 * x[1]^2 - x[2] 
    dx[3] = -ds.r * ds.z + ds.r * ds.s * (x[1] + 1.6) 
end 

##### Calculate Lyapunov exponents 

function lyapunovs(ds::Dynamics, nsteps::Int, dt::Real=0.01)
    # Construct integrator 
    integ = initinteg(ds) 

    # Advance integrator 
    t0 = ds.t 
    d = length(ds.x0) 
    λ = zeros(d) 
    for k in 1 : nsteps 
        step!(integ, dt, true) 
        Δ = integ.u[:, 2 : end] 
        QR = qr(Δ) 
        integ.u[:, 2 : end] .= QR.Q 
        λ .+= log.(abs.(diag(QR.R)))
        integ.u[:, 2 : end] .= QR.Q
    end 
    λ ./ (integ.t - t0)
end 

function initinteg(ds::Dynamics) 
    d = size(ds.x0, 1)
    J0 = ds(zeros(d, d), ds.x0, nothing, t)     # Initial jacobian matrix 
    Δ0 = diagm(ones(d))
    Φ0 = [ds.x0 Δ0]
    tangentf = augment(ds, J0)
    tspan = (ds.t, Inf)
    prob = ODEProblem(tangentf, Φ0, tspan)
    init(prob, SOLVER)
end 

function augment(ds::Dynamics, J)
    function f(dΦ, Φ, J, t) 
        x  = @view Φ[:, 1] 
        dx = @view dΦ[:, 1] 
        Δ  = @view Φ[:, 2 : end] 
        dΔ = @view dΦ[:, 2 : end]
        ds(dx, x, u, t)     # Update ds 
        ds(J, x, u, t)      # Update J 
        dΔ .= J * Δ         # Update dΔ
    end 
end 

##### Exports 
export Lorenz, Chua, Rossler, HRNeuron, lyapunovs

end # module 
