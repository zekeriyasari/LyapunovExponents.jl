```@meta
EditURL = "<unknown>/example/src/msfplots.jl"
```

# Master Stability Function (MSF)

### Network
Consider the network given as
$$
\dot{x}_i = f(x_i) - ϵ \sum_{j = 1}^n ξ_{ij} P x_j \quad i = 1, 2, \ldots, n
$$
where
* $ x_i ∈ \mathbb{R}^d $
* $ f : \mathbb{R}^d \mapsto \mathbb{R}^d $
* $ Ξ = [\xi_{ij}] $
* $ 0 = λ_1 < λ_2 ≤ λ_3 ≤ … ≤ λ_n $ are the eigenvalues of $Ξ$.

### Full Synchronization
The system achieves full synchronzation if
$$
\lim_{t \mapsto \infty} |x_i(t) - x_j(t)| \quad ∀ i
$$

### Master Stability Function
The full synchronization of the network given above can be associated with
$$
\dot{ζ} = (Df(s) - ϵ λ_k P) ζ \quad k = 1, 2, \ldots, n
$$
where $s$ is the special solution of the node dynamics $\dot{s} = f(s)$, $Df(s)$ is the jacobian of $f$ calculated at
$s$. For the system to achieve full synchronization, maximum lyapunov exponents calculated usign the auxilary system
must be negative. Let $η = ϵ λ$. Then, the master stability fucntion $Λ = Ψ(η)$ where $Λ$ is the maximum lyapunnov
expoent corresponding to $η$. In the following, we plot MSF.

Load packages

```@example msfplots
using LyapunovExponents
using Plots
```

 Parameters

```@example msfplots
ds = Lorenz()
P = [
    0 0 0;
    0 0 0;
    0 0 1
]
```

MSF calculation parameters

```@example msfplots
nsteps = Int(3e4)       # Number of steps
ntrsteps = Int(1e4)     # Number of transient steps
dt = 0.001              # Step size of the integration
```

Define $Ψ$

```@example msfplots
Ψ(η) = msf(ds, η, P, nsteps=nsteps, ntrsteps=ntrsteps, dt=dt)
```

Calculate MSF values

```@example msfplots
η = range(0, 100, length=51)
Λ = Ψ.(η)
```

Plots

```@example msfplots
plot(λ, Λ, marker=(:circle, 3))
hline!([0], linestyle=:dash)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

