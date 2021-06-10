using DifferentialEquations
using LinearAlgebra
using Plots 

function lorenz!(dx::AbstractVector, x::AbstractVector, u, t, σ=10, β=8/3, ρ=28) 
    dx[1] = σ * (x[2] - x[1]) 
    dx[2] = x[1] * (ρ - x[3]) - x[2] 
    dx[3] = x[1] * x[2] - β * x[3] 
end 

function lorenzjac!(J::AbstractMatrix, x::AbstractVector, u, t, σ=10, β=8/3, ρ=28) 
    J .= [
        -σ σ 0; 
        ρ - x[3] -1 -x[1]; 
        x[2] x[1] -β
    ]
end 

function tangentf(dΦ, Φ, J::AbstractMatrix, t) 
    x = @view Φ[:, 1] 
    dx = @view dΦ[:, 1] 
    Δ = @view Φ[:, 2 : end] 
    dΔ = @view dΦ[:, 2 : end] 
    lorenz!(dx, x, nothing, t)      # Update dx 
    lorenzjac!(J, x, nothing, t)    # Update J 
    dΔ .= J * Δ                     # Update dΔ
end 

function tangentf(dΦ, Φ, (J, λ, P)::Tuple, t) 
    Δ = @view Φ[:, 2 : end]
    dΔ = @view dΦ[:, 2 : end]
    tangentf(dΦ, Φ, J, t) 
    dΔ .-= λ * P * Δ
end 

function initinteg(f, fjac, msfparams...) 
    d = 3
    tspan = (0., Inf) 
    x0 = rand(d) 
    J0 = fjac(zeros(d, d), x0, nothing, 0.)
    Δ0 = diagm(ones(d))
    Φ0 = [x0 Δ0] 
    params = isempty(msfparams) ? J0 : (J0, msfparams...)
    prob = ODEProblem(tangentf, Φ0, tspan, params)
    init(prob, Tsit5())
end 

# Lyapunov calculation should start with an initialized with a tangent integrator. For other method calls such as calling
# lyapunov for a system defined with f, the corresponding tangent integrator must be prepared first.

function lyap(f::Function, fjac::Function; kwargs...)
    integ = initinteg(f, fjac) 
    _lyap(integ; kwargs...)
end

function msf(f::Function, fjac::Function, λ::Real, P::AbstractMatrix; kwargs...) 
    integ = initinteg(f, fjac, λ, P) 
    _lyap(integ; kwargs...) |> maximum
end 

function _lyap(integ; nsteps=Int(3e4), ntrsteps=Int(1e4), dt=0.001)
    d = size(integ.u, 1)
    for i in 1 : ntrsteps
        step!(integ, dt, true)
        Δ = integ.u[:, 2 : end]
        QR = qr(Δ) 
        integ.u[:, 2 : end] .= QR.Q
    end 
    Λ = zeros(d) 
    t0 = integ.t 
    for i in 1 : nsteps
        step!(integ, dt, true)
        Δ = integ.u[:, 2 : end]
        QR = qr(Δ) 
        Λ .+= log.(abs.(diag(QR.R)))
        integ.u[:, 2 : end] .= QR.Q
    end 
    T = integ.t - t0 
    Λ ./ T
end 


# -------------------------------- Example -------------------------------- #
  
function runlyapex() 
    f = lorenz!
    fjac = lorenzjac!
    dt = 0.001
    nsteps = Int(3e4) 
    ntrsteps = Int(1e4) 
    vals = lyap(f, fjac, nsteps=nsteps, ntrsteps=ntrsteps, dt=dt)
end 

function runmsfex(λ) 
    f = lorenz!
    fjac = lorenzjac!
    λ = range(1, 100, length=11)
    P = [
        1 0 0; 
        0 0 0; 
        0 0 0
        ]
    dt = 0.001
    nsteps = Int(3e4) 
    ntrsteps = Int(1e4) 
    Λ = map(λi -> msf(f, fjac, λi, P, nsteps=nsteps, ntrsteps=ntrsteps, dt=dt), λ)
    plot(λ, Λ)
end 

# runlyapex()
runmsfex()
