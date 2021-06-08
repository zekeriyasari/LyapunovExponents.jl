# This file includes methods for master stability function 

struct CoupledSystem{T1 <: Dynamics, T2<:Real, T3<:AbstractMatrix} <: Dynamics
    ds::T1 
    λ::T2 
    P::T3
    x0::Vector{Float64} 
    t::Real
end 
CoupledSystem(ds::Dynamics, λ::Real, P::AbstractMatrix) = CoupledSystem(ds, λ, P, ds.x0, ds.t)

function (cds::CoupledSystem)(dx, x, u, t)
    cds.ds(dx, x, u, t) 
    dx .-= λ * cds.P * dx
end 


msf(ds::Dynamics, λ::Real, P::AbstractMatrix) = CoupledSystem(ds, λ, P) |> lyapunovs
