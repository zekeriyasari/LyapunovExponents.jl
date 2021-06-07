# This file includes some well-known dynamics 

const SOLVER = Tsit5()

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
        x[2]            x[1]    -ds.β
    ]
end


Base.@kwdef struct Rossler <: Dynamics 
    α::Float64 = 0.2 
    β::Float64 = 0.2 
    γ::Float64 = 9
    t::Float64 = 0.
    x0::Vector{Float64} = rand(3) 
end 

function (ds::Rossler)(dx::AbstractVector, x::AbstractVector, u, t)
    dx[1] = -x[2] - x[3] 
    dx[2] = x[1] + ds.α * x[2] 
    dx[3] = ds.β + (x[1] - ds.γ) * x[3]    
end

function (ds::Rossler)(J::AbstractMatrix, x::AbstractVector, u, t)  # Jacobian 
    J .= [ 
        0       -1      -1; 
        1       ds.α    0; 
        x[3]    0       x[1] - ds.γ
    ]
end

Base.@kwdef struct Chen <: Dynamics 
    a::Float64 = 35 
    b::Float64 = 8/3   
    c::Float64 = 28
    t::Float64 = 0.
    x0::Vector{Float64} = rand(3) 
end 

function (ds::Chen)(dx::AbstractVector, x::AbstractVector, u, t)
    dx[1] = ds.a * (x[2] - x[1]) 
    dx[2] = x[1] * (ds.c - ds.a - x[3]) + ds.c * x[2] 
    dx[3] = x[1] * x[2] - ds.b * x[3] 
end

function (ds::Chen)(J::AbstractMatrix, x::AbstractVector, u, t)  # Jacobian 
    J .= [ 
        -ds.a                   ds.a       0; 
        ds.c - ds.a - x[3]      ds.c        -x[1]; 
        x[2]                    x[1]        -ds.b
    ]
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

function (ds::Chua)(dx::AbstractVector, x::AbstractVector, u, t) 
    dx[1] = ds.α * (x[2] - x[1] + ds.diode(x[1]))
    dx[2] = x[1] - x[2] + x[3] 
    dx[3] = -ds.β * x[2] 
end 

# TODO: Define jabobian functions for Chua dynamics.

Base.@kwdef struct HRNeuron <: Dynamics 
    r::Float64 = 0.006 
    s::Float64 = 4 
    I::Float64 = 3.2 
    t::Float64 = 0.
    x0::Vector{Float64} = rand(3)
end 

function (ds::HRNeuron)(dx::AbstractVector, x::AbstractVector, u, t) 
    dx[1] = x[2] + 3 * x[1]^2 - x[1]^3 - x[3] + ds.I 
    dx[2] = 1 - 5 * x[1]^2 - x[2] 
    dx[3] = -ds.r * x[3] + ds.r * ds.s * (x[1] + 1.6) 
end 

function (ds::HRNeuron)(J::AbstractMatrix, x::AbstractVector, u, t)  # Jacobian 
    J .= [ 
        6 * x[1] - 3 * x[1]^2   1           -1; 
        -10 * x[1]              -1          0; 
        ds.r * ds.s             0           -ds.r
    ]
end
