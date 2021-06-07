# This file include methods to calculate lyapunov exponents 

"""
    $SIGNATURES 

Computes the lyapunov exponents of `ds` for `nsteps` with a step size of `dt` seconds. 
"""
function lyapunovs(ds::Dynamics; nsteps::Int=Int(1e5), ntrsteps::Int=0, dt::Real=0.001)
    integ = initinteg(ds) 
    ntrsteps == 0 || transientsteps!(integ, ntrsteps, dt)
    t0 = integ.t        # Initial time 
    d = length(ds.x0)   # State space dimension 
    λ = zeros(d)        # Lyapunov exponents
    for k in 1 : nsteps 
        # While advacing the integrator, we perform QR decomposition to avoid numeical instabilities. 
        step!(integ, dt, true) 
        Δ = integ.u[:, 2 : end] 
        QR = qr(Δ) 
        integ.u[:, 2 : end] .= QR.Q 
        λ .+= log.(abs.(diag(QR.R)))
        integ.u[:, 2 : end] .= QR.Q
    end 
    T = integ.t - t0 
    λ ./ T
end 

# Takes transient steps to calculate the lyapunov exponents on the attractor.
function transientsteps!(integ, ntrsteps, dt)
    for k in 1 : ntrsteps
        step!(integ, dt, true)     
        Δ = integ.u[:, 2 : end] 
        QR = qr(Δ) 
        integ.u[:, 2 : end] .= QR.Q
    end 
end 

function augment(ds::Dynamics, J::AbstractMatrix)
    #= 
        Returns augmented dynamics given as 
        ̇x = f(x) 
        ̇Δ = Df(s)⋅Δ
        where Δ is the variation in the tangent space. 
        We collect the variables x and Δ with the following notation Φ = [x | Δ]
    =#
    function f(dΦ, Φ, J, t) 
        x  = @view Φ[:, 1] 
        dx = @view dΦ[:, 1] 
        Δ  = @view Φ[:, 2 : end] 
        dΔ = @view dΦ[:, 2 : end]
        ds(dx, x, nothing, t)     # Update ds 
        ds(J, x, nothing, t)      # Update J 
        dΔ .= J * Δ               # Update dΔ
    end 
end 

# Initialized augmented integrator of augmented system. 
function initinteg(ds::Dynamics) 
    d = dimension(ds) 
    J0 = ds(zeros(d, d), ds.x0, nothing, ds.t)  # Initial jacobian matrix 
    Δ0 = diagm(ones(d))                         # Initial variation in the tangent space 
    Φ0 = [ds.x0 Δ0]
    tangentf = augment(ds, J0)                  # Tangent space dynamics 
    # Since a large number of steps may be required for the calculation of lyapunov exponents, the final time is adjusted to
    # Inf here. 
    tspan = (ds.t, Inf)
    prob = ODEProblem(tangentf, Φ0, tspan, J0)
    init(prob, SOLVER)
end 
