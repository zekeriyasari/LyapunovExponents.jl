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

Master stability function (msf). Returns maximum lyapunov exponent of the network defined as 
```math 
\\dot{x}_i = f(x)_i + \\sum_{j = 1}^n \\epsilon_{ij} P x_j \\quad i = 1, 2, \\ldots, n
```
where ``x_i \\in \\mathbb{R}^d`` is the state vector of the ``i^{th}`` node, ``f`` is the vector-valued function defining the individual
node dynamics. ``\\epsilon_{ij} \\geq 0 `` if there is a connection between the nodes ``i`` and ``j``.
``\\epsilon_{ij} = 0``, if there is no connection between the nodes ``i`` and ``j``. ``E = [\\epsilon_{ij}]`` determines
he network topology. ``P = diag(p_1, p_2, \\ldots, p_d)`` determines the variables by which the nodes are conected to
each other. The network can be expressed more compaclty as
```math 
\\dot{X} = F(X) + (E \\otimes P) X 
```
where ``X = [x_1, x_2, \\ldots, x_n], F(X) = [f(x_1), f_(x_2), \\ldots, f_(x_n)], E = [\\epsilon_{ij}], P = diag(p_1,
p_2, \\ldots, p_n)``
    
`kwargs` are 
* nsteps::Int = 3e4 
    Number of steps to calculate lypaunov exponents 
* ntrsteps::Int = 0
    Number of transient steps for transient steps. Transient steps are taken before calculation of lypaunov exponents 
* dt::Real = 0.01
    Step size of integrator. 

!!! note 
    Note that ``\\epsilon_{ij} \\geq 0, j \\neq i`, ``E`` has nonpositive eigenvalues. Let ``\\eta`` denote the
    eigenvalue of ``E``, than we have, 
    ```math 
        \\nu_n \\leq \\nu_{n - 1} \\leq \\ldots \\leq \\nu_2 < \\nu_1 = 0
    ```
"""
function msf(ds::Dynamics, η::Real, P::AbstractMatrix; kwargs...)
    η ≤ 0 || error("Expected nonpositive eigenvalue η, got $η instead.")
    _lyapunovs(tangentinteg(ds, η, P); kwargs...) |> maximum
end 

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

function tangentinteg(ds::Dynamics, η::Real, P::AbstractMatrix) 
    function tf(dΦ, Φ, (J, η, P), t) 
        x  = @view Φ[:, 1] 
        Δ  = @view Φ[:, 2 : end] 
        dx = @view dΦ[:, 1] 
        dΔ = @view dΦ[:, 2 : end]   # Variations in the tangent space 
        ds(dx, x, nothing, t)   # Update dx 
        ds(J, x, nothing, t)    # Update J 
        dΔ .= (J + η * P) * Δ   # Update dΔ
    end 
    d = dimension(ds) 
    J = zeros(d, d) 
    ds(J, ds.x0, nothing, ds.t) 
    prob = ODEProblem(tf, [ds.x0 diagm(ones(d))], (ds.t, Inf), (J, η, P))
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
    η = zeros(size(integ.u, 1)) 
    for i in 1 : nsteps 
        step!(integ, dt, true) 
        Δ = integ.u[:, 2 : end] 
        QR = qr(Δ) 
        η .+= log.(abs.(diag(QR.R)))
        integ.u[:, 2 : end] .= QR.Q
    end 
    η ./ (integ.t - t0)
end


