using LyapunovExponents
using Plots

ds = Lorenz()
P = [
    0 0 0;
    0 0 0;
    0 0 1
]

nsteps = Int(3e4)       # Number of steps
ntrsteps = Int(1e4)     # Number of transient steps
dt = 0.001              # Step size of the integration

Ψ(η) = msf(ds, η, P, nsteps=nsteps, ntrsteps=ntrsteps, dt=dt)

η = range(0, 100, length=51)
Λ = Ψ.(η)

theme(:default)
plot(η, Λ, marker=(:circle, 3))
hline!([0], linestyle=:dash)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

