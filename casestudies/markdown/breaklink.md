```@meta
EditURL = "<unknown>/example/src/breaklink.jl"
```

# Mutating The Topology
In this notebook, we will consider mutating the topology of a network while simulating the network. For this purpose,
we will consider a very simple network consisting of a just two nodes whose dynamics are Lorenz dynamics defined by
$$
f([x, y, x]) = \left\{
\begin{aligned}
\dot{x} &= \sigma (y - x) \\
\dot{y} &= x (\rho - z) - y \\
\dot{z} & = x y - \beta x \\
\end{aligned}
\right.
$$
The network can be represented more compactly as
$$
\dot{X} = F(X) + \epsilon (\Xi \otimes P) X
$$
where $ \Xi = [\xi_{ij}] $ is the outer coupling matrix that is determined by the network topology. For this simple
network we have
$$
\Xi = \begin{bmatrix}
-1 &  1 \\
 1 & -1 \\
\end{bmatrix}
In this use-case, we will control the coupling between the nodes.

Load the packages

```@example breaklink
using DifferentialEquations
using Plots
using LightGraphs, GraphPlot
using LinearAlgebra
```

Construct the network ODE problem

```@example breaklink
function lorenz(dx, x, u, t, σ=10, β=8/3, ρ=35)
    dx[1] = σ * (x[2] - x[1])
    dx[2] = x[1] * (ρ - x[3]) - x[2]
    dx[3] = x[1] * x[2] - β * x[3]
end

mutable struct Net2{T1}
    f::T1
    ϵ::Float64
    Ξ::Matrix{Float64}
    P::Matrix{Float64}
end

function (net::Net2)(dx, x, u, t)
    for i in Iterators.partition(1 : length(x), size(net.P, 1))
        net.f(view(dx, i), view(x, i), nothing, t)
    end
    dx .+= net.ϵ[1] * kron(net.Ξ, net.P) * x
end

ϵl = 0.     # Low level of coupling. The nodes are decoupled.
ϵh = 50.    # High level of coupling. The nodes are coupled .
ϵ = ϵh
graph = path_graph(2)
Ξ = -collect(laplacian_matrix(graph)) .|> float
P = diagm([1, 0, 0]) .|> float
n = size(Ξ, 1)
d = size(P, 1)
net = Net2(lorenz, ϵ, Ξ, P)
prob = ODEProblem(net, rand(n * d), (0., 100.));
nothing #hide
```

For the control of the connection, we will construct to callbacks. We will start the simulation with the nodes
connected. Then, with the first callback, we will decoupled the nodes at $ t = 25 $. Then, the with the second
callback, we will couple the nodes again. When the nodes are coupled, they synchronize to each other, and when they
are decoupled the nodes desynchronize. Let us contruct the callbacks.

```@example breaklink
tc = 25.    # Decoupling callback time.
cond_couple(x, t, integ) = t == tc
act_couple!(integ) = (integ.f.f.ϵ = ϵl)
td = 40.    # Coupling callback time.
cond_decouple(x, t, integ) = t == td
act_decouple!(integ) = (integ.f.f.ϵ = ϵh)
cbs = CallbackSet(
    DiscreteCallback(cond_couple, act_couple!),
    DiscreteCallback(cond_decouple, act_decouple!)
);
nothing #hide
```

Then, we solve the problem

```@example breaklink
sol = solve(prob, callback=cbs, tstops=[tc, td])
```

Now we can plot the time waveform of the nodes.

```@example breaklink
theme(:default)
t, x = sol.t, sol.u
plot(t, abs.(getindex.(x, 1) - getindex.(x, 4)), label="e12")
vline!([tc],linewidth = 2, linestyle=:dash, label="decoupling")
vline!([td],linewidth = 2, linestyle=:dash, label="recoupling")
xlabel!("t")
```

From the plot, we can see that the nodes synchronize when start the simulation. Then at $ t = 25 $, the coupling
between the nodes are disconnected. After this the nodes desynchronize. The, at $ t = 40 $ the nodes are connected
again, and the synchronization error goes to zero shortly .

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

