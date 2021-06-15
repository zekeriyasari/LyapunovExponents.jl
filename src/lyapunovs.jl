# This file includes methods for master stability function.

export lyapunovs, msf

"""
    $SIGNATURES 

Computes lyapunov exponents of `ds`.kwargs may include 

* nsteps::Int = 3e4 
    Number of steps to calculate lypaunov exponents 
* ntrsteps::Int = 0
    Number of transient steps for transient steps. Transient steps are taken before calculation of lypaunov exponents 
* dt::Real = 0.01
    Step size of integrator. 
"""
lyapunovs(ds::Dynamics; kwargs...) = _lyapunovs(tangentinteg(ds); kwargs...)

"""
    $SIGNATURES 

Master stability function (msf). Returns maximum lyapunov exponent. kwargs may include 

* nsteps::Int = 3e4 
    Number of steps to calculate lypaunov exponents 
* ntrsteps::Int = 0
    Number of transient steps for transient steps. Transient steps are taken before calculation of lypaunov exponents 
* dt::Real = 0.01
    Step size of integrator. 
"""
msf(ds::Dynamics, λ::Real, P::AbstractMatrix; kwargs...) = _lyapunovs(tangentinteg(ds, λ, P); kwargs...) |> maximum

function tangentinteg(ds::Dynamics) 
    function tf(dΦ, Φ, J, t) 
        x  = @view Φ[:, 1] 
        Δ  = @view Φ[:, 2 : end] 
        dx = @view dΦ[:, 1] 
        dΔ = @view dΦ[:, 2 : end]   
        ds(dx, x, nothing, t)   # Update dx 
        ds(J, x, nothing, t)    # Update J 
        dΔ .= J * Δ             # Update dΔ
    end 
    d = dimension(ds) 
    J = zeros(d, d) 
    ds(J, ds.x0, nothing, ds.t) 
    prob = ODEProblem(tf, [ds.x0 diagm(ones(d))], (ds.t, Inf), J)
    init(prob, Tsit5())
end 

function tangentinteg(ds::Dynamics, λ::Real, P::AbstractMatrix) 
    function tf(dΦ, Φ, (J, λ, P), t) 
        x  = @view Φ[:, 1] 
        Δ  = @view Φ[:, 2 : end] 
        dx = @view dΦ[:, 1] 
        dΔ = @view dΦ[:, 2 : end]   # Variations in the tangent space 
        ds(dx, x, nothing, t)   # Update dx 
        ds(J, x, nothing, t)    # Update J 
        dΔ .= (J - λ * P) * Δ   # Update dΔ
    end 
    d = dimension(ds) 
    J = zeros(d, d) 
    ds(J, ds.x0, nothing, ds.t) 
    prob = ODEProblem(tf, [ds.x0 diagm(ones(d))], (ds.t, Inf), (J, λ, P))
    init(prob, Tsit5())
end 

function _lyapunovs(integ; nsteps::Int=Int(3e4), ntrsteps::Int=0, dt::Real=0.01)
    # Take transient steps 
    for i in 1 : ntrsteps
        step!(integ, dt, true) 
        Δ = integ.u[:, 2 : end] 
        QR = qr(Δ) 
        integ.u[:, 2 : end] .= QR.Q
    end 

    # Calculate lyapunov exponents 
    t0 = integ.t
    Λ = zeros(size(integ.u, 1)) 
    for i in 1 : nsteps 
        step!(integ, dt, true) 
        Δ = integ.u[:, 2 : end] 
        QR = qr(Δ) 
        Λ .+= log.(abs.(diag(QR.R)))
        integ.u[:, 2 : end] .= QR.Q
    end 
    Λ ./ (integ.t - t0)
end
