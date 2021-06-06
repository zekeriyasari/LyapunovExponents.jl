
using LyapunovExponents
nsteps = 10000
dt = 0.01
for ds in [Lorenz(), Chen(), Rossler(), HRNeuron()]
    println("$ds ", " => ", lyapunovs(ds, nsteps, dt))
end 
