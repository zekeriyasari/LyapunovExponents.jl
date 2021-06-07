# This file prints lyapunov expoenents of some of well-known systems. 

using LyapunovExponents

nsteps = 10000  # Number of steps 
dt = 0.01       # Step size

for ds in [Lorenz(), Chen(), Rossler(), HRNeuron()]
    println("$ds ", " => ", lyapunovs(ds, nsteps, dt))
end 
