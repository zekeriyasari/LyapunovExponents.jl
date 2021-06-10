# This file plots master stability function plots 

using LyapunovExponents 
using Plots 

ds = Lorenz() 
P = [
    0 0 0; 
    0 0 0; 
    0 0 1
]
nsteps = Int(3e4)
ntrsteps = Int(1e4) 
dt = 0.001 
λ = range(0, 100, length=51) 
Λ = map(λi -> msf(ds, λi, P, nsteps=nsteps, ntrsteps=ntrsteps, dt=dt), λ) 
@info "Done calculation"
plot(λ, Λ, marker=(:circle, 3))
hline!([0], linestyle=:dash)
