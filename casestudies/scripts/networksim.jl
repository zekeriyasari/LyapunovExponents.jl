using LightGraphs
using GraphPlot

n = 5   # Number of nodes in the network
gr = cycle_graph(n)
add_edge!(gr, 1, 3)
add_edge!(gr, 1, 4)
add_edge!(gr, 2, 5)
gplot(gr, nodelabel=1:nv(gr), layout=circular_layout)

using LyapunovExponents
using LinearAlgebra
using Plots

ds = Lorenz()

P = diagm([1, 0, 0])
d = size(P, 1)  # Node dimension

Ψ(η) = msf(ds, η, P)

η = range(0, 100, length=21)
Λ = Ψ.(η)
plot(η, Λ, marker=(:circle, 3))
hline!([0.], color="red")

ϵ = 100
Ξ = collect(laplacian_matrix(gr))
λ = eigvals(Ξ)
Λnet = Ψ.(ϵ * λ[2 : end])

net = Network(Lorenz, -ϵ * Ξ, P)
sol = solvenet(net, (0., 100.), saveat=0.001)
t, x = sol.t, sol.u;

x1 = getindex.(x, 1)
plt = plot(layout = (n - 1, 1))
for i in 2 : n
    xi = getindex.(x, (i - 1) * d + 1)
    plot!(t, abs.(x1 - xi), label = "e1$i", subplot=i - 1)
end
display(plt)

Q = [
    1 0 0 0 0;
    0 0 0 0 1;
    0 0 0 1 0;
    0 0 1 0 0;
    0 1 0 0 0
]

M = eigvecs(Q)
G = inv(M) * Ξ * M
B = G[1 : 2, 1 : 2]
D = G[3 : 5, 3 : 5]
λB = eigvals(B)
λD = eigvals(D)
@show λB
@show λD;

@show Ψ.(ϵ * λB)
@show Ψ.(ϵ * λD)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

