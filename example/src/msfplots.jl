

# # Master Stability Function (MSF)

# ### Network 
# Consider the network given as 
# $$ 
# \dot{x}_i = f(x_i) - \epsilon \sum_{j = 1}^n \xi_{ij} P x_j \quad i = 1, 2, \ldots, n 
# $$ 
# where 
# * $ x_i \in \mathbb{R}^d $ 
# * $ f : \mathbb{R}^d \mapsto \mathbb{R}^d $ 
# * $ Ξ = [\xi_{ij}] $ 
# * $ 0 = \lambda_1 < \lambda_2 \leq \lambda_3 \leq \ldots \leq \lambda_n $ are the eigenvalues of $\Xi$. 

# ### Full Synchronization 
# The system achieves full synchronzation if 
# $$ 
# \lim_{t \mapsto \infty} |x_i(t) - x_j(t)| \quad \forall i 
# $$ 

# ### Master Stability Function 
# The full synchronization of the network given above can be associated with 
# $$ 
# \dot{\zeta} = (Df(s) - \epsilon \lambda_k P) \zeta \quad k = 1, 2, \ldots, n 
# $$  
# where $s$ is the special solution of the node dynamics $\dot{s} = f(s)$, $Df(s)$ is the jacobian of $f$ calculated at
# $s$. For the system to achieve full synchronization, maximum lyapunov exponents calculated usign the auxilary system
# must be negative. Let $\eta = \epsilon \lambda$. Then, the master stability fucntion $\Lambda = \Psi(\eta)$ where $\Lambda$ is the maximum lyapunnov
# expoent corresponding to $\eta$. In the following, we plot MSF. 
#-
# Load packages
using LyapunovExponents 
using Plots 

#  Parameters 
ds = Lorenz() 
P = [
    0 0 0; 
    0 0 0; 
    0 0 1
]

# MSF calculation parameters 
nsteps = Int(3e4)       # Number of steps 
ntrsteps = Int(1e4)     # Number of transient steps 
dt = 0.001              # Step size of the integration

# Define $Ψ$ 
Ψ(η) = msf(ds, η, P, nsteps=nsteps, ntrsteps=ntrsteps, dt=dt)

# Calculate MSF values 
η = range(0, 100, length=51) 
Λ = Ψ.(η)

# Plots 
theme(:default) 
plot(η, Λ, marker=(:circle, 3))
hline!([0], linestyle=:dash)

# #### References 
# * Huang, L., Chen, Q., Lai, Y. C., & Pecora, L. M. (2009). Generic behavior of master-stability functions in coupled nonlinear dynamical systems. Physical Review E, 80(3), 036204.
