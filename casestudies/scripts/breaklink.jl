using DifferentialEquations
using Plots
using LightGraphs, GraphPlot
using LinearAlgebra

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

sol = solve(prob, callback=cbs, tstops=[tc, td])

theme(:default)
t, x = sol.t, sol.u
plot(t, abs.(getindex.(x, 1) - getindex.(x, 4)), label="e12")
vline!([tc],linewidth = 2, linestyle=:dash, label="decoupling")
vline!([td],linewidth = 2, linestyle=:dash, label="recoupling")
xlabel!("t")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

